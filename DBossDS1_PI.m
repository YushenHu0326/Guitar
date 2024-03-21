classdef DBossDS1_PI < audioPlugin
    %DBossDS1: Simplified Model of the Boss-DS1 guitar distortion pedal
    
    properties
        dist = 0.001, tc = 0.4, lvl = 1, os = false, fs;
    end
    
    properties (Constant, Hidden)
        PluginInterface = audioPluginInterface('PluginName', 'DBossDS1', ...
        audioPluginParameter('dist', 'DisplayName', 'Distortion gain',...
        'Mapping', {'lin', 0.001, 1}, 'Style', 'rotaryknob', 'Layout', [1 2]),...
        audioPluginParameter('tc', 'DisplayName', 'Tone control',...
        'Mapping', {'lin', 0, 1}, 'Style', 'rotaryknob', 'Layout', [1 1]),...
        audioPluginParameter('lvl', 'DisplayName', 'Level',...
        'Mapping', {'pow', 2, 0, 2}, 'Style', 'rotaryknob', 'Layout', [3 1]),...
         audioPluginParameter('os',...
                'DisplayName','Oversampling', ...
                'Mapping', {'enum', 'None', '2x'}, ...
                'Style', 'vrocker', 'Layout', [3 2]), ...
        audioPluginGridLayout('RowHeight', [160 20 100 20], ...
            'ColumnWidth', [100 100], 'Padding', [10 30 10 20]))
    end
    
    % Define private properties: filter coefficients and filter state for
    % real time processing
    properties (Access = private)
        lowpassState = zeros(2);
        highpassState = zeros(2);
        BJTstate = zeros(2);
        OpAmpState = zeros(2);
        lowpassAdjState = zeros(2);
        blendState = zeros(2);
        bBJT = zeros(1,3);
        aBJT = zeros(1,3);
        bOpAmp = zeros(1,3);
        aOpAmp = zeros(1,3);
        bLow = zeros(1,3);
        bHigh = zeros(1,3);
        aLow = zeros(1,3);
        aHigh = zeros(1,3);
        bLowBnd = zeros(1,3);
        aLowBnd = zeros(1,3);
        twoDy = zeros(1024, 2);

    end
    
    methods
        function plugin = DBossDS1_PI  
            % Constructor; init coeffs
            fs = getSampleRate(plugin);
            % BJT tranzistor amplification stage coefficients
            [plugin.bBJT, plugin.aBJT] = calcBJTcoef(plugin, fs);
            % Op-Amp stage coefficients
            [plugin.bOpAmp, plugin.aOpAmp] = calcOpAmpCoef(plugin, plugin.dist, fs);
            % Filter coefficients for the two lowpass and one highpass
            % filters
            [plugin.bLow, plugin.aLow, plugin.bHigh, plugin.aHigh, plugin.bLowBnd, plugin.aLowBnd] = calcFiltCoef(plugin, fs);
        end
        
        function y = process(plugin, x)
            % Audio processing function
            
            % BJT transistor gain stage
            [y, plugin.BJTstate] = filter(plugin.bBJT, plugin.aBJT, x, plugin.BJTstate);
            
            % Filter to reduce harshness; blend in 20%
            [blendLp, plugin.blendState] = filter(plugin.bLow, plugin.aLow, y, plugin.blendState);
            
            y = 0.2*blendLp + y;
            
            % Op-Amp gain stage
            [y, plugin.OpAmpState] = filter(plugin.bOpAmp, plugin.aOpAmp, y, plugin.OpAmpState);
  
            % Resampling           
            if plugin.os == true
                % Upsample signal and select first two columns containing
                % the stereo channels
                x_os = resample(y, 2, 1);
                x_os = x_os(:,[1,2]);
                % Apply LPF filter 
                %[y, plugin.resState] = filter(plugin.bRes,plugin.aRes, x_us, plugin.resState);
                % Apply clipping stage; n = 2.5
                x_clip = (x_os ./ (1 + abs(x_os).^2.5).^(1./2.5))*0.5 + sign(x_os).*(1-exp(-abs(x_os)))*0.5;            
                % Downsample the signal and select the first two columns
                y = resample(x_clip, 1, 2);
                y = y(:,[1,2]);
            else
                   % Clipping stage; n = 2.5 -> hard-coded; dictates curve
                   % of the exponential function
                    y = (y ./ (1 + abs(y).^2.5).^(1./2.5))*0.5 + sign(y).*(1-exp(-abs(y)))*0.5; 
            end
  
            % Tone control stage
            [lowPassRC, plugin.lowpassState] = filter(plugin.bLow, plugin.aLow, y, plugin.lowpassState);
            [highPassRC, plugin.highpassState] = filter(plugin.bHigh, plugin.aHigh, y, plugin.highpassState);

            % Blend between the two filters 
            y = highPassRC*plugin.tc + lowPassRC*(1-plugin.tc);
            
            % Lowpass filter to smooth the higher frequencies a little bit more
            [y, plugin.lowpassAdjState] = filter(plugin.bLowBnd, plugin.aLowBnd, y, plugin.lowpassAdjState);
            
            % Apply level gain to the signal and tanh to limit
            y = 1.2 * y * plugin.lvl;
            y = tanh(y);            
        end
        
        function set.dist(plugin, val)
           plugin.dist = val;
           [plugin.bOpAmp, plugin.aOpAmp] = calcOpAmpCoef(plugin, val, getSampleRate(plugin));
        end
        
        function reset(plugin)
            sr = getSampleRate(plugin);
            % Recompute coefficients dependent on sampling rate
            [plugin.bBJT, plugin.aBJT] = calcBJTcoef(plugin, sr);
            [plugin.bOpAmp, plugin.aOpAmp] = calcOpAmpCoef(plugin, plugin.dist, sr);
            [plugin.bLow, plugin.aLow, plugin.bHigh, plugin.aHigh, plugin.bLowBnd, plugin.aLowBnd] = calcFiltCoef(plugin, sr);
            
            plugin.lowpassState = zeros(2);
            plugin.highpassState = zeros(2);
            plugin.BJTstate = zeros(2);
            plugin.OpAmpState = zeros(2);
            plugin.blendState = zeros(2);
            plugin.lowpassAdjState = zeros(2);
        end
       
    end
    
    methods (Access = private)
            
        function [bBJT, aBJT] = calcBJTcoef(~, fs)
            % BJT gain stage
            
            % discrete time coefficient
            coeff = pi/fs/2;

            % Discretet-time denominator coefficients
            B0 = 1;
            B1 = -2;
            B2 = 1;

            % Discrete-time numerator coefficients
            A0 = 7200.*coeff.^2 + 1206.*coeff + 1;
            A1 = 14400.*coeff.^2 - 2;
            A2 = 7200.*coeff.^2 - 1206.*coeff + 1;

            % Stage needs 36dB of bandbass gain
            % Equation for gain 36dB = 20*log10(x)
            amp = 10.^(36/20);
            bBJT = amp.*[B0, B1, B2];
            aBJT = [A0, A1, A2];
        end
        
        function [bOpAmp, aOpAmp] = calcOpAmpCoef(~, dist, fs)
            % Op-Amp gain stage

            % Resistors and capacitors values
            Rt = 100000 * dist + 0.001;
            Rb = 100000*(1-dist) + 4700;
            Cz = 0.00000047;
            Cc = 0.000000000100;

            % Constant for the bilinear transform
            c = 2*fs;

            % Continuous-time coefficients
            % a0 = b0
            b0 = 1 / (Rt*Cc*Rb*Cz);
            a1 = 1/(Rb*Cz) + 1/(Rt*Cc);
            b1 = a1 + 1/(Rb*Cc);
            a1 = 1/(Rb*Cz) + 1/(Rt*Cc);

            % Discrete-time coefficients
            B0 = b0 + b1*c + c.^2;
            B1 = 2*b0 - 2*c.^2;
            B2 = b0 - b1*c + c.^2;

            A0 = b0 + a1*c + c.^2;
            A1 = B1;
            A2 = b0 - a1*c + c.^2;

            bOpAmp = [B0, B1, B2];
            aOpAmp = [A0, A1, A2]; 
        end
        
        function [bLow, aLow, bHigh, aHigh, bLowBnd, aLowBnd] = calcFiltCoef(~, fs)
            % Filter coefficients
            
            % Lowpass Adj coeff
            w0LpA    = 2*pi*2300/fs;
            alphaLpA = sin(w0LpA)/sqrt(2);
            cosw0LpA = cos(w0LpA);
            normLpA  = 1/(1+alphaLpA);
            
            % Lowpass coeff
            w0Lp    = 2*pi*234/fs;
            alphaLp = sin(w0Lp)/sqrt(2);
            cosw0Lp = cos(w0Lp);
            normLp  = 1/(1+alphaLp);
            
            % Highpass coeff
            w0Hp   = 2*pi*1063/fs;
            alphaHp = sin(w0Hp)/sqrt(2);
            cosw0Hp = cos(w0Hp);
            normHp  = 1/(1+alphaHp);
            
            % HP
            bHigh = (1 + cosw0Hp)*normHp * [.5 -1 .5];
            aHigh = [1 -2*cosw0Hp*normHp (1 - alphaHp)*normHp];
            % LP
            bLow = (1 - cosw0Lp)*normLp * [.5 1 .5];
            aLow = [1 -2*cosw0Lp*normLp (1 - alphaLp)*normLp];
            % LP Adj
            bLowBnd = (1 - cosw0LpA)*normLpA * [.5 1 .5];
            aLowBnd = [1 -2*cosw0LpA*normLpA (1 - alphaLpA)*normLpA];
            
        end
        
    end
        
end

