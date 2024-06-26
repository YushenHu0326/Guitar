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
            current_df=(i-ceil(t/2*(60/BPM)/dt))/(ceil((t+duration)/2*(60/BPM)/dt)-ceil(t/2*(60/BPM)/dt))*2;
            if current_df>1
                current_df=1;
            end
            current_f=floor(f+current_df*(f0-f));
            fp(s,i)=ceil(fret(current_f)*J);
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

BPM=126;

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
tmax=39; %total time of the simulation in seconds
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

%play_note(1,1,1,2);
%play_note(2,1,1,2);
%play_note(3,1,2,4);
%play_note(4,1,1,2);
%play_note(5,1,2,5);
%play_note(6,1,1,2);
%play_note(7,1,2,4);
%play_note(8,1,1,2);
%play_note(9,1,2,2);
%play_note(10,1,1,5);
%play_note(11,1,1,4);
%play_note(12,1,1,5);
%play_note(13,1,2,2);
%play_note(14,1,1,5);
%play_note(15,1,1,4);
%play_note(16,1,1,0);

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

play_note(3.5,0.5,4,-1);
play_note(4,0.25,4,-1);
%%%
play_note(4.25,0.25,5,7);
play_hammeron(4.5,0.25,5,8);
play_pulloff(4.75,0.25,5,7);
play_note(5,0.5,4,9);
play_note(5.5,0.5,5,7);
play_note(6,1,5,10);
bend(6,1,5,10,12);
play_note(7,0.5,6,8);
play_pulloff(7.5,0.5,6,7);
play_note(8,0.5,6,11);
play_pulloff(8.5,0.5,6,7);
play_note(9,0.5,6,8);
play_note(9.5,0.5,6,11);
play_slide(10,1,6,12,14);
play_note(11,0.5,6,11);
play_note(11.5,0.5,6,12);
%%%
play_slide(12,1,6,14,15);
play_note(13,0.5,6,12);
play_note(13.5,0.5,6,14);
play_note(14,0.5,6,15);
play_pulloff(14.5,0.5,6,14);
play_note(15,0.5,6,12);
play_note(15.5,0.5,5,15);
play_note(16,0.5,6,14);
play_pulloff(16.5,0.5,6,11);
play_note(17,0.25,6,12);
play_hammeron(17.25,0.25,6,14);
play_pulloff(17.5,0.5,6,12);
play_note(18,0.5,5,15);
play_note(18.5,0.5,6,12);
play_note(19,1,5,15);
bend(19,1,5,15,17);
%%%
play_note(20,2,5,15);
bend(20,2,5,15,17);
play_note(22,2,5,15);
bend(22,2,5,15,17);
play_note(24,2,5,15);
bend(24,2,5,15,17);
play_note(26,1,5,15);
bend(26,1,5,15,17);
play_note(27,0.5,5,12);
play_note(27.5,0.5,4,14);
%%%
play_note(28,2,5,15);
bend(28,2,5,15,17);
play_note(30,2,5,15);
bend(30,2,5,15,17);
play_note(32,2,5,15);
bend(32,2,5,15,17);
play_note(34,1,5,15);
bend(34,1,5,15,17);
play_note(35,0.5,5,12);
play_note(35.5,0.5,4,14);
%%%
play_note(36,2,5,17);
bend(36,2,5,17,19);
play_note(38,4,5,17);
bend(38,2,5,17,19);
play_note(42,2,5,17);
bend(42,2,5,17,19);
play_note(42,2,4,17);
bend(44,2,4,17,19);
%%%
play_note(44,3,5,17);
bend(44,3,5,17,19);
play_note(47,1,6,15);
play_note(48,2,5,17);
bend(48,2,5,17,19);
play_note(50,1,5,17);
bend(50,1,5,17,19);
play_note(51,1,5,15);
%%%
play_note(52,3,5,17);
vibrato(52,3,5,17);
play_note(55,2,5,15);
bend(55,2,5,15,17);
play_note(57,1,6,12);
play_note(58,1,5,15);
play_note(59,1,5,12);
%%%
play_note(60,2,4,14);
bend(60,2,4,14,15);
play_note(62,0.5,4,12);
play_note(62.5,1.5,4,14);
vibrato(62.5,1.5,4,14);
play_note(64,0.5,4,12);
play_note(64.5,1.5,4,14);
vibrato(64.5,1.5,4,14);
play_note(66,1,4,12);
play_note(67,2,4,14);
bend(67,2,4,14,16);
%%%
play_note(69,1,4,14);
play_note(70,1,4,12);
play_note(71,1,4,14);
play_note(72,1,4,12);
play_note(73,1,4,14);
bend(73,1,4,14,16);
play_note(73.5,0.5,5,15);
play_note(74,1,4,14);
bend(74,1,4,14,16);
play_note(74.5,0.5,4,12);
play_note(75,1,4,14);
bend(75,1,4,14,16);
play_note(75.5,0.5,5,15);
%%%
play_chord(76,[1.5,1.5],[5,6],[15,15]);
play_note(76.5,1,4,14);
bend(76.5,1,4,14,16);
play_chord(77.5,[1.5,1.5],[5,6],[15,15]);
play_note(78,1,4,14);
bend(78,1,4,14,16);
play_note(79,0.5,6,15);
play_note(79.5,0.5,5,15);
play_note(80,1,4,14);
bend(80,1,4,14,16);
play_note(81,0.5,4,12);
play_note(81.5,0.5,3,14);
play_note(82,0.5,4,14);
play_pulloff(82.5,0.5,4,12);
play_chord(83,[0.5,0.5],[1,2],[-1,-1]);
play_chord(83.5,[0.5,0.5],[1,2],[-1,-1]);
%%%
play_note(84,2,4,12);
play_note(86,0.5,3,14);
play_pulloff(86.5,0.5,3,12);
play_note(87,2,3,14);
play_note(89,0.5,2,-1);
play_note(89.5,0.5,2,-1);
play_note(90,0.5,4,12);
play_hammeron(90.5,0.5,4,14);
play_pulloff(91,0.5,4,12);
play_note(91.5,0.5,3,14);
%%%
play_note(92,2,4,12);
play_note(94,0.5,3,14);
play_pulloff(94.5,0.5,3,12);
play_note(95,2,3,14);
play_note(97,0.5,2,-1);
play_note(97.5,0.5,2,-1);
play_note(98,0.5,4,12);
play_hammeron(98.5,0.5,4,14);
play_pulloff(99,0.5,4,12);
play_note(99.5,0.5,3,14);
%%%
play_note(100,0.5,4,12);
play_note(100.5,0.5,4,12);
play_note(101,0.5,3,14);
play_pulloff(101.5,0.5,3,12);
play_note(102,0.5,4,14);
play_note(102.5,0.5,4,14);
play_pulloff(103,0.5,4,12);
play_note(103.5,0.5,4,14);
play_note(104,0.5,5,12);
play_note(104.5,0.5,5,12);
play_note(105,0.5,4,14);
play_pulloff(105.5,0.5,4,12);
play_note(106,0.5,5,15);
play_note(106.5,0.5,5,15);
play_pulloff(107,0.5,5,12);
play_note(107.5,0.5,4,14);
%%%
play_note(108,0.5,6,12);
play_note(108.5,0.5,6,12);
play_note(109,0.5,5,15);
play_pulloff(109.5,0.5,5,12);
play_note(110,0.5,6,15);
play_note(110.5,0.25,6,15);
play_pulloff(110.75,0.25,6,12);
play_note(111,0.5,6,15);
play_note(111.5,0.5,6,14);
play_note(112,0.5,5,15);
play_note(112.5,0.25,6,12);
play_hammeron(112.75,0.25,6,15);
play_pulloff(113,0.5,6,12);
play_note(113.5,0.5,5,15);
%%%
play_note(114,0.5,5,15);
bend(114,0.5,5,15,17);
play_note(114.5,0.5,6,12);
play_note(115,0.5,5,15);
play_pulloff(115.5,0.5,5,12);
play_note(116,0.5,4,14);
bend(116,0.5,4,14,16);
play_note(116.5,0.5,6,12);
play_note(117,0.5,5,15);
play_pulloff(117.5,0.5,5,12);
play_note(118,0.5,5,15);
bend(118,0.5,5,15,17);
play_note(118.5,0.5,6,12);
play_note(119,0.5,5,15);
play_pulloff(119.5,0.5,5,12);
play_note(120,0.5,4,14);
bend(120,0.5,4,14,16);
play_note(120.5,0.5,6,12);
play_note(121,0.5,5,15);
play_pulloff(121.5,0.5,5,12);
%%%
play_note(122,0.5,5,15);
bend(122,0.5,5,15,17);
play_note(122.5,0.5,6,12);
play_note(123,0.5,5,15);
play_pulloff(123.5,0.5,5,12);
play_note(124,0.5,4,14);
bend(124,0.5,4,14,16);
play_note(124.5,0.5,6,12);
play_note(125,0.5,5,15);
play_pulloff(125.5,0.5,5,12);
play_note(126,0.5,5,15);
bend(126,0.5,5,15,17);
play_note(126.5,0.5,6,12);
play_note(127,0.5,5,15);
play_pulloff(127.5,0.5,5,12);
play_note(128,0.5,4,14);
bend(128,0.5,4,14,16);
play_note(128.5,0.5,6,12);
play_note(129,0.5,5,15);
play_pulloff(129.5,0.5,5,12);
%%%
play_note(130,4,6,15);
bend(130,4,6,15,17);
play_note(134,1.5,6,15);
play_note(135.5,0.5,6,12);
play_hammeron(136,0.5,6,15);
play_pulloff(136.5,0.5,6,12);
play_note(137,0.5,5,15);
play_note(137.5,0.5,6,12);
%%%
play_note(138,0.5,5,15);
bend(138,0.5,5,15,17);
play_note(139.5,0.5,6,12);
play_note(140,0.5,5,15);
play_pulloff(140.5,0.5,5,12);
play_note(141,0.5,4,14);
bend(141,0.5,4,14,16);
play_note(141.5,0.5,4,12);
play_note(142,1,4,14);
bend(142,1,4,14,16);
play_note(143,0.5,5,15);
play_note(143.5,1,4,14);
bend(143.5,1,4,14,16);
play_note(144.5,0.5,4,12);
play_note(145,0.5,3,14);
play_note(145.5,0.5,4,14);
play_pulloff(146,0.5,4,12);
play_note(146.5,1.5,3,14);
play_slide(148,4,4,12,2);

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

S=gdist(0.9995,S);

soundsc(S(1:count))
%plot the soundwave as a function of time:
%figure
plot(tsave(1:count),S(1,1:count))
%plot(STRINGT(6,:))
%plot(NOISE);

end
