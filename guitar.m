% TO DO LIST
% MAKE GUITAR SOUND -check
% PLAY SINGLE NOTE -check
% RELEASE & LET RING -check
% TREMOLO
% BENDING
% VIBRATO
% TEMPO SYSTEM -check
% CHORD -check
% PALM MUTE -check
% HAMMER ON & PULL OFF
% AH
% PH
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
        TIMESTAMP(s,ceil((t+duration)/2*60/BPM*dt))=-1;
    elseif(f==0)
        TIMESTAMP(s,ceil(t/2*(60/BPM)/dt))=25;
        TIMESTAMP(s,ceil((t+duration)/2*60/BPM*dt))=-1;
    end
end

function play_chord(t,duration,s,f)
    interval=0;
    for i=1:size(s,2)
        TIMESTAMP(s(i),ceil((t+interval)/2*(60/BPM)/dt))=f(i);
        TIMESTAMP(s(i),ceil((t+duration+interval)/2*(60/BPM)/dt))=-1;
        interval=interval+0.1;
    end
end

function left_palm_mute(t,duration,s,f)
    if(f>0 && f<25)
        TIMESTAMP(s,ceil(t/2*(60/BPM)/dt))=f+60;
        TIMESTAMP(s,ceil(t+duration)/2*(60/BPM)/dt)=-1;
    end
end

function left_palm_mute_chord(t,duration,s,f)
    interval=0;
    for i=1:size(s,2)
        if(f(i)>0 && f(i)<25)
            TIMESTAMP(s(i),ceil((t+interval)/2*(60/BPM)/dt))=f(i)+60;
            TIMESTAMP(s(i),ceil(t+duration+interval)/2*(60/BPM)/dt)=-1;
            interval=interval+0.1;
        end
    end
end

function right_palm_mute(t,duration,s,f)
    if(f>0 && f<25)
        TIMESTAMP(s,ceil(t/2*(60/BPM)/dt))=f+30;
        TIMESTAMP(s,ceil((t+duration)/2*60/BPM*dt))=-1;
    elseif(f==0)
        TIMESTAMP(s,ceil(t/2*(60/BPM)/dt))=55;
        TIMESTAMP(s,ceil((t+duration)/2*60/BPM*dt))=-1;
    end
end

function play_bend(t,duration,s,f,f1)
    if(((f>0 && f<25) && (f1>0 && f1<25)) && f1>f)
        TIMESTAMP(s,ceil(t/2*(60/BPM)/dt))=f;
        TIMESTAMP(s,ceil((t+duration)/2*60/BPM*dt))=-1;
        for i=ceil(t/2*(60/BPM)/dt):ceil((t+duration)/2*(60/BPM)/dt)
            if(s==1)
                %STRINGT(i)=
            end
        end
    end
end

function play(s,f)
    if(f==0)
        fp(s)=0;
    else
        fp(s)=ceil(fret(f)*J);
    end
    %initial conditions for plucked string:
    xp=L*pickpos;
    Hp=.2; %position and amplitude of pluck

    for jj=fp(s)+1:J
        xpp=fp(s)*dx;
        x=(jj-1)*dx;
        if(x<xp)
            H(s,jj)=Hp*(x-xpp)/(xp-xpp);
        else
            H(s,jj)=Hp*(L-x)/(L-xp);
        end
    end
end

function release(s)
    fp(s)=0;
    H(s,:)=0;
    V(s,:)=V(s,:)/8;
end

clear all

BPM=120;

%frequencies for all 6 strings
f=[82,110,147,196,247,330];
%a copy of the initial frequencies
f_init=f;
%fret distance chart
fret=zeros(1,24);
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
%Currently pressed fret
fp=zeros(1,6);

%string parameters to make frequency f1:
L=100;
M=1;
T=zeros(1,6);
for ii=1:6
    T(ii)=M*(2*L*f(ii))^2;
end

tau=1.2; %decay time (seconds)
%damping constant to make decay time tau:
R=zeros(1,6);
for ii=1:6
    R(ii)=(2*M*L^2)/(tau*pi^2);
end
%a copy of initial damping constant
R_init=R;

J=81;
dx=L/(J-1);
%maximum time step for numerical stability:
%dtmax=zeros(1,6);
dtmax=-(R(1)/T(1))+sqrt((R(1)/T(1))^2+(dx^2/(T(1)/M)));
for ii=1:6
    %dtmax(ii)=-(R/T(ii))+sqrt((R/T(ii))^2+(dx^2/(T(ii)/M)));
    n_dtmax=-(R(ii)/T(ii))+sqrt((R(ii)/T(ii))^2+(dx^2/(T(ii)/M)));
    if(n_dtmax<dtmax)
        dtmax=n_dtmax;
    end
end
%Now set dt and nskip such that:
%dt<=dtmax, nskip is a positive integer, and dt*nskip = 1/8192.
%Also, make nskip as small as possible, given the above criteria.
nskip=ceil(1/(8192*dtmax));
dt=1/(8192*nskip);
tmax=10; %total time of the simulation in seconds
clockmax=ceil(tmax/dt);

% Action on Timstamp:
% <= 24: Play note
% 25: Open note
% -1: Release
% note+30: Right palm mute
% note+60: Left palm mute
TIMESTAMP=zeros(6,clockmax);

% When not 0, bend the string using the corresponding tension
STRINGT=zeros(6,clockmax);

H=zeros(6,J);
V=zeros(6,J);

%Initialize pickup position
pickup = 0.81;
pickup_2 = 0.82;
pickup_3 = 0.87;
pickpos = 0.9;


%Part of the main riff of Sweet Child O'Mine to demonstrate single note
%play_note(1,1,3,12);
%play_note(2,1,5,15);
%play_note(3,1,4,14);
%play_note(4,1,4,12);
%play_note(5,1,6,15);
%play_note(6,1,4,14);
%play_note(7,1,6,14);
%play_note(8,1,4,14);

%A simple power chord to demonstrate chord, notice the order of the string
%decides if downpicking or not
%play_chord(1,1,[1,2,3],[5,7,7]);
left_palm_mute_chord(1,0.5,[1,2,3],[5,7,7]);
left_palm_mute_chord(1.5,0.5,[1,2,3],[5,7,7]);
play_chord(2,0.5,[1,2,3],[5,7,7]);
play_chord(2.5,0.5,[1,2,3],[5,7,7]);

count=0;

S=zeros(1,ceil(clockmax/nskip));
tsave=zeros(1,ceil(clockmax/nskip));

for clock=1:clockmax
    t=clock*dt;
    for str=1:6
        if(TIMESTAMP(str,clock)>0)
            %Play note
            if(TIMESTAMP(str,clock)<25)
                play(str,TIMESTAMP(str,clock));
            %Play open note
            elseif(TIMESTAMP(str,clock)==25)
                play(str,0);
            %Play right palm mute
            elseif(TIMESTAMP(str,clock)>30 && TIMESTAMP(str,clock)<55)
                play(str,TIMESTAMP(str,clock)-30);
                R(str)=(2*M*L^2)/(0.2*pi^2);
            %Play right palm mute on open note
            elseif(TIMESTAMP(str,clock)==55)
                play(str,0);
                R(str)=(2*M*L^2)/(0.2*pi^2);
            %Play left palm mute
            elseif(TIMESTAMP(str,clock)>60 && TIMESTAMP(str,clock)<85)
                play(str,TIMESTAMP(str,clock)-60);
                R(str)=(2*M*L^2)/(0.06*pi^2);
            end
        elseif(TIMESTAMP(str,clock)<0)
            if(TIMESTAMP(str,clock)==-1)
                release(str);
                R(str)=R_init(str);
            end
        end
        j=fp(str)+2:(J-1); % list of indices of interior points
        V(str,j)=V(str,j)+(dt/dx^2)*(T(str)/M)*(H(str,j+1)-2*H(str,j)+H(str,j-1))+(dt/dx^2)*(R(str)/M)*(V(str,j+1)-2*V(str,j)+V(str,j-1));
        H(str,j)=H(str,j)+dt*V(str,j);
    end
    if(mod(clock,nskip)==0)
        count=count+1;
        S(count)=sum(H(:,ceil(J*pickup))); %sample the sound
        S(count)=S(count)+sum(H(:,ceil(J*pickup_2))); %sample the sound
        tsave(count)=t; %record sample time
    end
    %set(Hhandle,'ydata',H) %update movie frame
    %drawnow %show latest frame
end

S=gdist(0.99,S);

soundsc(S(1:count))
%plot the soundwave as a function of time:
%figure
%plot(tsave(1:count),S(1,1:count))

end
