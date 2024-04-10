% TO DO LIST
% MAKE GUITAR SOUND -check
% PLAY SINGLE NOTE -check
% RELEASE & LET RING -check
% TREMOLO -check
% BENDING -check
% VIBRATO -check
% TEMPO SYSTEM -check
% CHORD -check
% PALM MUTE -check
% HAMMER ON & PULL OFF
% AH -check
% PH -check
% STOCHASTIC PROCESS
% HUMBUCKER PICKUP

function guitar()

function [x]=gdist(a,x)
    k = 2*a/(1-a);
    x = (1+k)*(x)./(1+k*abs(x));
end

function play_note(t,duration,s,f)
    if(f>0 && f<25)
        TIMESTAMP(s,ceil(t/2*(60/BPM)/dt))=f;
        TIMESTAMP(s,ceil((t+duration)/2*60/BPM/dt)-2)=-1;
        for i=ceil(t/2*(60/BPM)/dt):ceil((t+duration)/2*(60/BPM)/dt)
            fp(s,i)=ceil(fret(f)*J);
        end
    elseif(f==0)
        TIMESTAMP(s,ceil(t/2*(60/BPM)/dt))=25;
        TIMESTAMP(s,ceil((t+duration)/2*60/BPM/dt)-2)=-1;
    elseif(f==-1)
        TIMESTAMP(s,ceil(t/2*(60/BPM)/dt))=59;
        TIMESTAMP(s,ceil((t+duration)/2*60/BPM/dt)-2)=-1;
    end
end

function play_slide(t,duration,s,f,f0)
    if(f>0 && f<25)
        TIMESTAMP(s,ceil(t/2*(60/BPM)/dt))=f;
        TIMESTAMP(s,ceil((t+duration)/2*60/BPM/dt)-2)=-1;
        for i=ceil(t/2*(60/BPM)/dt):ceil((t+duration)/2*(60/BPM)/dt)
            fp(s,i)=ceil(fret(ceil(f+(f0-f)*i/(ceil((t+duration)/2*(60/BPM)/dt)-ceil(t/2*(60/BPM)/dt))))*J);
        end
    end
end

function play_bottle_slide(t,duration,s,f,f0)
    if(f>0 && f<25)
        TIMESTAMP(s,ceil(t/2*(60/BPM)/dt))=f;
        TIMESTAMP(s,ceil((t+duration)/2*60/BPM/dt)-2)=-1;
        for i=ceil(t/2*(60/BPM)/dt):ceil((t+duration)/2*(60/BPM)/dt)
            fp(s,i)=ceil((fret(f)+(fret(f0)-fret(f))*i/(ceil((t+duration)/2*(60/BPM)/dt)-ceil(t/2*(60/BPM)/dt)))*J);
        end
    end
end

function play_tremolo(t,duration,s,f)
    if(f>0 && f<25)
        TIMESTAMP(s,ceil(t/2*(60/BPM)/dt))=f;
        TIMESTAMP(s,ceil((t+duration*0.25)/2*(60/BPM)/dt))=f;
        TIMESTAMP(s,ceil((t+duration*0.5)/2*(60/BPM)/dt))=f;
        TIMESTAMP(s,ceil((t+duration*0.75)/2*(60/BPM)/dt))=f;
        TIMESTAMP(s,ceil((t+duration)/2*60/BPM/dt)-2)=-1;
        for i=ceil(t/2*(60/BPM)/dt):ceil((t+duration)/2*(60/BPM)/dt)
            fp(s,i)=ceil(fret(f)*J);
        end
    elseif(f==0)
        TIMESTAMP(s,ceil(t/2*(60/BPM)/dt))=25;
        TIMESTAMP(s,ceil((t+duration*0.25)/2*(60/BPM)/dt))=25;
        TIMESTAMP(s,ceil((t+duration*0.5)/2*(60/BPM)/dt))=25;
        TIMESTAMP(s,ceil((t+duration*0.75)/2*(60/BPM)/dt))=25;
        TIMESTAMP(s,ceil((t+duration)/2*60/BPM/dt)-2)=-1;
    elseif(f==-1)
        TIMESTAMP(s,ceil(t/2*(60/BPM)/dt))=59;
        TIMESTAMP(s,ceil((t+duration*0.25)/2*(60/BPM)/dt))=59;
        TIMESTAMP(s,ceil((t+duration*0.5)/2*(60/BPM)/dt))=59;
        TIMESTAMP(s,ceil((t+duration*0.75)/2*(60/BPM)/dt))=59;
        TIMESTAMP(s,ceil((t+duration)/2*60/BPM/dt)-2)=-1;
    end
end

function play_chord(t,duration,s,f)
    interval=0;
    for i=1:size(s,2)
        if(f(i)>0 && f(i)<25)
            TIMESTAMP(s(i),ceil(t/2*(60/BPM)/dt)+interval)=f(i);
            TIMESTAMP(s(i),ceil((t+duration(i))/2*(60/BPM)/dt)-2)=-1;
            interval=interval+1;
            for iii=ceil(t/2*(60/BPM)/dt):ceil((t+duration(i))/2*(60/BPM)/dt)
                fp(s(i),iii)=ceil(fret(f(i))*J);
            end
        elseif(f(i)==-1)
            TIMESTAMP(s(i),ceil(t/2*(60/BPM)/dt)+interval)=59;
            TIMESTAMP(s(i),ceil((t+duration(i))/2*(60/BPM)/dt)-2)=-1;
            interval=interval+1;
        end
    end
end

function right_palm_mute(t,duration,s,f)
    if(f>0 && f<25)
        TIMESTAMP(s,ceil(t/2*(60/BPM)/dt))=f+30;
        TIMESTAMP(s,ceil((t+duration)/2*60/BPM/dt)-2)=-1;
        for i=ceil(t/2*(60/BPM)/dt):ceil((t+duration)/2*(60/BPM)/dt)
            fp(s,i)=ceil(fret(f)*J);
        end
    elseif(f==0)
        TIMESTAMP(s,ceil(t/2*(60/BPM)/dt))=55;
        TIMESTAMP(s,ceil((t+duration)/2*60/BPM/dt)-2)=-1;
    end
end

function bend(t,duration,s,f,f1)
    if(((f>0 && f<25) && (f1>0 && f1<25)) && f1>f)
        for i=ceil(t/2*(60/BPM)/dt):ceil((t+duration)/2*(60/BPM)/dt)
            if(s==1)
                df=5+f*0.5;
            elseif(s==2)
                df=6+f*0.6;
            elseif(s==3)
                df=10+f*0.8;
            elseif(s==4)
                df=12+f*1.15;
            elseif(s==5)
                df=15+f*1.45;
            elseif(s==6)
                df=20+f*2;
            end
            tfactor=(i-ceil(t/2*(60/BPM)/dt))*2/(ceil((t+duration)/2*(60/BPM)/dt)-ceil(t/2*(60/BPM)/dt));
            if(tfactor>1)
                tfactor=1;
            end
            STRINGT(s,i)=tfactor*df*(f1-f)/4;
        end
    end
end

function vibrato(t,duration,s,f)
    if(f>0 && f<25)
        for i=ceil(t/2*(60/BPM)/dt):ceil((t+duration)/2*(60/BPM)/dt)
            STRINGT(s,i)=abs(sin((i-ceil(t/2*(60/BPM)/dt))/(ceil((t+duration)*(60/BPM)/dt)-ceil(t/2*(60/BPM)/dt))*20)*5);
        end
    end
end

function play_hammeron(t,duration,s,f)
    if(f>0 && f<25)
        TIMESTAMP(s,ceil(t/2*(60/BPM)/dt))=f+60;
        TIMESTAMP(s,ceil((t+duration)/2*60/BPM/dt)-2)=-1;
        for i=ceil(t/2*(60/BPM)/dt):ceil((t+duration)/2*(60/BPM)/dt)
            fp(s,i)=ceil(fret(f)*J);
        end
    end
end

function play_pulloff(t,duration,s,f)
    if(f>0 && f<25)
        TIMESTAMP(s,ceil(t/2*(60/BPM)/dt))=f+90;
        TIMESTAMP(s,ceil((t+duration)/2*60/BPM/dt)-2)=-1;
        for i=ceil(t/2*(60/BPM)/dt):ceil((t+duration)/2*(60/BPM)/dt)
            fp(s,i)=ceil(fret(f)*J);
        end
    elseif(f==0)
        TIMESTAMP(s,ceil(t/2*(60/BPM)/dt))=115;
        TIMESTAMP(s,ceil((t+duration)/2*60/BPM/dt)-2)=-1;
    end
end

function AH(t,duration,s,f)
    if(f>0 && f<25)
        TIMESTAMP(s,ceil(t/2*(60/BPM)/dt))=25;
        TIMESTAMP(s,ceil((t+duration)/2*60/BPM/dt)-2)=-1;
        for i=ceil(t/2*(60/BPM)/dt):ceil((t+duration)/2*(60/BPM)/dt)
            HARMONICP(s,i)=f;
        end
    end
end

function PH(t,duration,s,f)
    if(f>0 && f<25)
        TIMESTAMP(s,ceil(t/2*(60/BPM)/dt))=f;
        TIMESTAMP(s,ceil((t+duration)/2*60/BPM/dt)-2)=-1;
        for i=ceil(t/2*(60/BPM)/dt):ceil((t+duration)/2*(60/BPM)/dt)
            HARMONICP(s,i)=25;
            fp(s,i)=ceil(fret(f)*J);
        end
    end
end

function play(s,f)
    %if(f==0)
    %    fp(s)=0;
    %else
    %    fp(s)=ceil(fret(f)*J);
    %end
    fp0(s)=f;
    %initial conditions for plucked string:
    xp(s)=L*pickpos;
    Hp(s)=.1; %position and amplitude of pluck
end

function hammeron(s,f)
    %fp(s)=ceil(fret(f)*J);
    fp0(s)=f;
    %initial conditions for plucked string:
    xp(s)=L*fret(f);
    Hp(s)=.02; %position and amplitude of pluck
end

function pulloff(s,f,f0)
    %if(f==0)
    %    fp(s)=0;
    %else
    %    fp(s)=ceil(fret(f)*J);
    %end
    fp0(s)=f;
    %initial conditions for plucked string:
    xp(s)=L*fret(f0);
    Hp(s)=.05; %position and amplitude of pluck
end

function release(s)
    Hp(s)=0.0002;
    damp(s)=20;
end

BPM=120;

%Initialize pickup position
pickup = 0.76;
pickup_2 = 0.765;
pickup_3 = 0.95;
pickpos = 0.77;

%frequencies for all 6 strings
f=[82,110,147,196,247,330];
%a copy of the initial frequencies
f_init=f;
%fret distance chart
fret=zeros(1,25);
fret(1)=0.056125;        fret(2)=0.10910185;
fret(3)=0.1591034;       fret(4)=0.20629938;
fret(5)=0.25084722;      fret(6)=0.29289352;
fret(7)=0.33258024691;   fret(8)=0.37004012345;
fret(9)=0.40539660493;   fret(10)=0.43876851851;
fret(11)=0.47026851851;  fret(12)=0.5;
fret(13)=0.5280632716;   fret(14)=0.55455092592;
fret(15)=0.57955246913;  fret(16)=0.60314969135;
fret(17)=0.6254228395;   fret(18)=0.64644598765;
fret(19)=0.66629012345;  fret(20)=0.68502006172;
fret(21)=0.70269753086;  fret(22)=0.71938425925;
fret(23)=0.73513425925;  fret(24)=0.75;
fret(25)=pickpos;

%string parameters to make frequency f1:
L=100;
M=1;
T=zeros(1,6);
for ii=1:6
    T(ii)=M*(2*L*f(ii))^2;
end

tau=2; %decay time (seconds)
%damping constant to make decay time tau:
R=(2*M*L^2)/(tau*pi^2);

%damping factor
damp=zeros(1,6);
for ii=1:6
    damp(ii)=20;
end

J=81;
dx=L/(J-1);
%maximum time step for numerical stability:
%dtmax=zeros(1,6);
dtmax=-(R/T(1))+sqrt((R/T(1))^2+(dx^2/(T(1)/M)));
for ii=1:6
    %dtmax(ii)=-(R/T(ii))+sqrt((R/T(ii))^2+(dx^2/(T(ii)/M)));
    n_dtmax=-(R/T(ii))+sqrt((R/T(ii))^2+(dx^2/(T(ii)/M)));
    if(n_dtmax<dtmax)
        dtmax=n_dtmax;
    end
end
%Now set dt and nskip such that:
%dt<=dtmax, nskip is a positive integer, and dt*nskip = 1/8192.
%Also, make nskip as small as possible, given the above criteria.
nskip=ceil(1/(8192*dtmax));
dt=1/(8192*nskip);
tmax=3; %total time of the simulation in seconds
clockmax=ceil(tmax/dt);

% Action on Timstamp:
% <= 24: Play note
% 25: Open note
% -1: Release
% note+30: Right palm mute
% 59: Left palm mute
% note+60: Hammer on
% note+90: Pull off
TIMESTAMP=zeros(6,clockmax);

% When not 0, bend the string using the corresponding tension
STRINGT=zeros(6,clockmax);

% When not 0, slightly touch the string to make harmonics
HARMONICP=zeros(6,clockmax);

NOISE=zeros(1,clockmax);
for ii=1:clockmax
    NOISE(ii)=(rand-0.5)/100+1;
end

%Currently pressed fret
fp=zeros(6,clockmax);
fp0=zeros(1,6);

lastPluckT=zeros(1,6);
for ii=1:6
    lastPluckT(ii)=tmax+1;
end

H=zeros(6,J);
H0=zeros(6,J);

neck_pickup=1;

xp=zeros(1,6);
Hp=zeros(1,6);
for ii=1:6
    xp(ii)=L*pickpos;
    Hp(ii)=0.1;
end

%-1 denotes x on the tab

%Part of the main riff of Sweet Child O'Mine to demonstrate single note
%play_note(1,16,3,12);
%play_note(2,1,5,15);
%play_note(3,1,4,14);
%play_note(4,1,4,12);
%play_note(5,1,6,15);
%play_note(6,1,4,14);
%play_note(7,1,6,14);
%play_note(8,1,4,14);

play_slide(1,1,5,12,15);

%A simple power chord to demonstrate chord, notice the order of the string
%decides if downpicking or not
%play_chord(1,[2,2,2],[1,2,3],[5,7,7]);
%right_palm_mute(3,1,1,5);
%right_palm_mute(4,1,1,5);
%right_palm_mute(5,1,1,5);
%right_palm_mute(6,1,1,5);
%right_palm_mute(7,1,1,5);
%right_palm_mute(8,1,1,5);

%play_chord(1,[2,2,2],[1,2,3],[-1,-1,-1]);

%double stop
%play_note(1,1,6,12);
%play_note(1,1,5,15);
%bend(1,1,5,15,17);

%play_tremolo(1,1,6,12);
%play_note(1,1,6,12);
%vibrato(1,1,6,12);

%A simple blue phrase to demonstrate hammer on and pull off
%play_note(1,1,3,12);
%play_hammeron(2,1,3,15);
%play_pulloff(3,1,3,12);
%play_hammeron(4,1,3,15);
%play_note(5,1,3,15);
%bend(5,1,3,15,17);

%AH(1,4,6,7);
%PH(1,1,3,4);
%bend(1,1,3,4,6);

count=0;

S=zeros(1,ceil(clockmax/nskip));
tsave=zeros(1,ceil(clockmax/nskip));

for clock=1:clockmax
    t=clock*dt;

    if(mod(clock,nskip)==0)
        count=count+1;
        if(neck_pickup==1)
            S(count)=sum(H(:,ceil(J*pickup))-H0(:,ceil(J*pickup))); %sample the sound
            S(count)=S(count)+sum(H(:,ceil(J*pickup_2))-H0(:,ceil(J*pickup_2))); %sample the sound
        else
            S(count)=sum(H(:,ceil(J*pickup_3))-H0(:,ceil(J*pickup_3)));
        end
        tsave(count)=t; %record sample time
    end

    for str=1:6
        if(TIMESTAMP(str,clock)>0)
            %Play note
            if(TIMESTAMP(str,clock)<25)
                play(str,TIMESTAMP(str,clock));
                lastPluckT(str)=t;
            %Play open note
            elseif(TIMESTAMP(str,clock)==25)
                play(str,0);
                lastPluckT(str)=t;
            %Play right palm mute
            elseif(TIMESTAMP(str,clock)>30 && TIMESTAMP(str,clock)<55)
                play(str,TIMESTAMP(str,clock)-30);
                lastPluckT(str)=t;
                damp(str)=100;
            %Play right palm mute on open note
            elseif(TIMESTAMP(str,clock)==55)
                play(str,0);
                lastPluckT(str)=t;
                damp(str)=100;
            %Play left palm mute
            elseif(TIMESTAMP(str,clock)==59)
                play(str,0);
                lastPluckT(str)=t;
                damp(str)=200;
            elseif(TIMESTAMP(str,clock)>59 && TIMESTAMP(str,clock)<85)
                hammeron(str,TIMESTAMP(str,clock)-60);
                lastPluckT(str)=t;
            elseif(TIMESTAMP(str,clock)>89 && TIMESTAMP(str,clock)<115)
                pulloff(str,TIMESTAMP(str,clock)-90,fp0(str));
                lastPluckT(str)=t;
            elseif(TIMESTAMP(str,clock)==115)
                pulloff(str,0,fp0(str));
                lastPluckT(str)=t;
            end
        elseif(TIMESTAMP(str,clock)<0)
            if(TIMESTAMP(str,clock)==-1)
                release(str);
                lastPluckT(str)=t;
            elseif(TIMESTAMP(str,clock)<-10 && TIMESTAMP(str,clock)>-40)
                play(str,0);
                lastPluckT(str)=t;
            end
        end
        T(str)=M*(2*L*(f_init(str)+STRINGT(str,clock)))^2;
        j=fp(str,clock)+2:(J-1); % list of indices of interior points
        H0(str,j)=H(str,j);
        H(str,j)=0;
        if(t>=lastPluckT(str))
            if(HARMONICP(str,clock)~=0)
                incr=round(1/min([fret(HARMONICP(str,clock)),1-fret(HARMONICP(str,clock))]));
            else
                incr=1;
            end
            for n=1:5
                H(str,j)=H(str,j)+(2*Hp(str)*(L-L*fp(str,clock)/J)^2*sin(xp(str)*n*incr*pi/(L-L*fp(str)/J)))/(xp(str)*(L-L*fp(str,clock)/J-xp(str))*n*incr*n*incr*pi*pi)*sin(n*incr*j/J*pi)*cos(n*incr*pi*(t-lastPluckT(str))*sqrt(T(str)/M)/(L-L*fp(str,clock)/J));
            end
            H(str,j)=H(str,j)/(damp(str)*(t-lastPluckT(str))+1);
            H(str,j)=H(str,j)*NOISE(clock);
        end
    end
    
    %set(Hhandle,'ydata',H) %update movie frame
    %drawnow %show latest frame
end

function wah(audio_in, fs)
    % Input: audio file
    damping = 0.06; %The lower it is, the smaller pass band
    width = 500;
    min_cutoff = 500; % Minimum cut off frequency 
    max_cutoff = 3000;% Maximum cut off frequency
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
    sound(normed, fs);
end

S=gdist(0.9995,S);

%soundsc(S(1:count))

wah(S(1:count),8192)

%plot the soundwave as a function of time:
%figure
plot(tsave(1:count),S(1,1:count))
%plot(STRINGT(6,:))
%plot(NOISE);

end
