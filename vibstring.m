%vibstring.m
clear all
f1=440 % frequency in Hz (cycles/second)
%string parameters to make frequency f1:
L=1;M=1;T=M*(2*L*f1)^2;
tau = 1.2; %decay time (seconds)
%damping constant to make decay time tau:
R= (2*M*L^2)/(tau*pi^2);
J=81;dx=L/(J-1);
%maximum time step for numerical stability:
dtmax=-(R/T)+sqrt((R/T)^2+(dx^2/(T/M)));
%Now set dt and nskip such that:
%dt<=dtmax, nskip is a positive integer, and dt*nskip = 1/8192.
%Also, make nskip as small as possible, given the above criteria.
nskip = ceil(1/(8192*dtmax));
dt=1/(8192*nskip);
tmax= 4; %total time of the simulation in seconds
clockmax=ceil(tmax/dt);
%initial conditions for plucked string:
V=zeros(1,J);
xp=L/3;Hp=1; %position and amplitude of pluck
for jj=1:J
x=(jj-1)*dx;
if(x<xp)
H(jj)=Hp*x/xp;
else
H(jj)=Hp*(L-x)/(L-xp);
end
end
%Initialize graphics to see a movie of the string:
set(gcf,'double','on') %turn on double-buffering
%make initial plot:
Hhandle=plot(0:dx:L,H,'linewidth',2);
axis([-0.2*L,L+0.2*L,-Hp,Hp]) %set axis limits
axis manual %freeze axes
drawnow
count=0;
S=zeros(1,ceil(clockmax/nskip));
tsave = zeros(1,ceil(clockmax/nskip));
j=2:(J-1); % list of indices of interior points
for clock=1:clockmax
t=clock*dt;
V(j)=V(j)+(dt/dx^2)*(T/M)*(H(j+1)-2*H(j)+H(j-1)) ...
+(dt/dx^2)*(R/M)*(V(j+1)-2*V(j)+V(j-1));
H(j)=H(j)+dt*V(j);
if(mod(clock,nskip)==0)
count=count+1;
S(count)=H(2); %sample the sound
tsave(count)=t; %record sample time
end
set(Hhandle,'ydata',H) %update movie frame
drawnow %show latest frame
end
soundsc(S(1:count)) %play the recorded soundwave
%plot the soundwave as a function of time:
figure(2)
plot(tsave(1:count),S(1:count))