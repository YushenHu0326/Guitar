function wah(audio)
    % Input: audio file
    [ audio_in, fs] = audioread(audio);
    damping = 0.05;
    width = 1000;
    min_cutoff = 250;
    max_cutoff = 5000;
    center_freq = width/fs;
    cutoff_freq=min_cutoff:center_freq:max_cutoff;
    while(length(cutoff_freq) < length(audio_in) )
        cutoff_freq = [ cutoff_freq (max_cutoff:-center_freq:min_cutoff) ];
        cutoff_freq = [ cutoff_freq (min_cutoff:center_freq:max_cutoff) ];
    end
    cutoff_freq = cutoff_freq(1:length(audio_in));
    % control the center frequency
    F1 = 2*sin((pi*cutoff_freq(1))/fs);
    Q1 = 2*damping;
    % Create and Zero Vectors to Match Length of Audio Input File
    highpass=zeros(size(audio_in));
    bandpass=zeros(size(audio_in));
    lowpass=zeros(size(audio_in));
    highpass(1) = audio_in(1);
    bandpass(1) = F1*highpass(1);
    lowpass(1) = F1*bandpass(1);
    for n=2:length(audio_in)
        highpass(n) = audio_in(n) - lowpass(n-1) - Q1*bandpass(n-1);
        bandpass(n) = F1*highpass(n) + bandpass(n-1);
        lowpass(n) = F1*bandpass(n) + lowpass(n-1);
        F1 = 2*sin((pi*cutoff_freq(n))/fs);
    end
    % Normalize and play back
    normed = bandpass./max(max(abs(bandpass)));
    audiowrite('wah wahed.wav', normed, fs);
    sound (normed, fs);
end
