function guitar()

function [x]=gdist(a,x)
    k = 2*a/(1-a);
    x = (1+k)*(x)./(1+k*abs(x));
end

function play(s,f)
    fp(s)=ceil(fret(f)*J);
    %initial conditions for plucked string:
    xp=L*pickpos;
    Hp=1; %position and amplitude of pluck

    for jj=fp(s)+1:J
        x=(jj-1)*dx;
        if(x<xp)
            H(s,jj)=Hp*x/xp;
        else
            H(s,jj)=Hp*(L-x)/(L-xp);
        end
    end
end

clear all

%frequencies for all 6 strings
f=[82,110,147,196,247,330];
f1=f(1); f2=f(2); f3=f(3); f4=f(4); f5=f(5); f6=f(6);
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
R=(2*M*L^2)/(tau*pi^2);
J=81;
dx=L/(J-1);
%maximum time step for numerical stability:
dtmax=zeros(1,6);
for ii=1:6
    dtmax(ii)=-(R/T(ii))+sqrt((R/T(ii))^2+(dx^2/(T(ii)/M)));
end
%Now set dt and nskip such that:
%dt<=dtmax, nskip is a positive integer, and dt*nskip = 1/8192.
%Also, make nskip as small as possible, given the above criteria.
nskip=zeros(1,6);
dt=zeros(1,6);
for ii=1:6
    nskip(ii)=ceil(1/(8192*dtmax(ii)));
    dt(ii)=1/(8192*nskip(ii));
end
tmax=10; %total time of the simulation in seconds
clockmax=zeros(1,6);
for ii=1:6
    clockmax(ii)=ceil(tmax/dt(ii));
end

H=zeros(6,J);
V=zeros(6,J);

%Initialize pickup position
pickup = 0.9;
pickpos = 0.9;

play(1,2);
play(2,7);
%play(3,7);

count=zeros(1,6);

S1=zeros(1,ceil(clockmax(1)/nskip(1)));
tsave1=zeros(1,ceil(clockmax(1)/nskip(1)));

S2=zeros(1,ceil(clockmax(2)/nskip(2)));
tsave2=zeros(1,ceil(clockmax(2)/nskip(2)));

S3=zeros(1,ceil(clockmax(3)/nskip(3)));
tsave3=zeros(1,ceil(clockmax(3)/nskip(3)));

S4=zeros(1,ceil(clockmax(4)/nskip(4)));
tsave4=zeros(1,ceil(clockmax(4)/nskip(4)));

S5=zeros(1,ceil(clockmax(5)/nskip(5)));
tsave5=zeros(1,ceil(clockmax(5)/nskip(5)));

S6=zeros(1,ceil(clockmax(6)/nskip(6)));
tsave6=zeros(1,ceil(clockmax(6)/nskip(6)));

for str=1:6
    j=fp(str)+2:(J-1); % list of indices of interior points
    for clock=1:clockmax
        t=clock*dt(str);
        V(str,j)=V(str,j)+(dt(str)/dx^2)*(T(str)/M)*(H(str,j+1)-2*H(str,j)+H(str,j-1))+(dt(str)/dx^2)*(R/M)*(V(str,j+1)-2*V(str,j)+V(str,j-1));
        H(str,j)=H(str,j)+dt(str)*V(str,j);
        if(mod(clock,nskip(str))==0)
            count(str)=count(str)+1;
            if(str==1)
                S1(count(str))=H(str,ceil(J*pickup)); %sample the sound
                tsave1(count(str))=t; %record sample time
            elseif(str==2)
                S2(count(str))=H(str,ceil(J*pickup));
                tsave2(count(str))=t;
            elseif(str==3)
                S3(count(str))=H(str,ceil(J*pickup));
                tsave3(count(str))=t;
            elseif(str==4)
                S4(count(str))=H(str,ceil(J*pickup));
                tsave4(count(str))=t;
            elseif(str==5)
                S5(count(str))=H(str,ceil(J*pickup));
                tsave5(count(str))=t;
            elseif(str==6)
                S6(count(str))=H(str,ceil(J*pickup));
                tsave6(count(str))=t;
            end
        end
        %set(Hhandle,'ydata',H) %update movie frame
        %drawnow %show latest frame
    end
end

S1=gdist(0.99,S1);

soundsc(S2(1:count))
%plot the soundwave as a function of time:
%figure
%plot(tsave(1:count),S(1,1:count))

end
