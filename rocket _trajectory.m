clc; clear; close all;
% all quantities are measured in the international system

% plot selection (1-on; 0-off):
Gv = 1; % the components of velocity as functions of time
Gc = 1; % laws of motion
Gt = 1; % trajectory (the ballistic curve)
Gd = 1; % dynamical representation
if Gd==1, Gv=0; Gc=0; Gt=0; end

% Physical parameters:
g = 9.80665; % gravitational acceleration
ro = 7850; % density of steel
r = 0.13; % radius of the projectile
m = 4/3 * pi * r^3 * ro; % mass of the projectile
G = m*g; % weight of the projectile

% Initial conditions:
v0 = 1100; % initial velocity
alpha0 = 43; % angle of launch
eta = 1.81 * 1e-5; % viscosity coefficient
b1 = 6 * pi * eta * r; % linear term coefficient
c = 0.469; % shape coefficient
ro0 = 1.22; % density of air
b2 = c * 4 * pi * r^2 * ro0/2; % quadratic term coefficient

% defining the time
t0 = 0; tf = 2 * v0 /g*sind(alpha0);
N = 1500; % number of moments of time
t = linspace(t0,tf,N); dt = t(2) - t(1);

% Starting values:
vx = zeros(1,N); vy = vx;
x = zeros(1,N); y = x;
vx(1) = v0 * cosd(alpha0);
vy(1) = v0 * sind(alpha0);
for i = 1:N-1
    aux = 1 - dt*(b1 + b2*sqrt(vx(i)^2 + vy(i)^2))/m;
    vx(i+1)=vx(i)*aux;
    vy(i+1)=vy(i)*aux - g*dt;
    x(i+1)=x(i) + vx(i)*dt;
    y(i+1)=y(i) + vy(i)*dt;
    if y(i+1)<0, break; end
end
t = t(1:i); vx=vx(1:i); vy=vy(1:i); x=x(1:i); y=y(1:i); % elimination of surplus values
if Gv==1 % laws of velocity
    figure(1);
    plot(t,vx,'-r',t,vy,'-b');
    xlabel('t(s)'); ylabel('v(m/s)'); grid;
    title('The components of velocity as functions of time');
    legend('vx','vy');
end
if Gc==1 % laws of motion
    figure(2);
    plot(t,x/1e3,'-r',t,y/1e3,'-b');
    xlabel('t(s)'); ylabel('coord(km)'); grid;
    title('Coordinates as functions of time');
    legend('x','y','Location','northwest');
end
if Gt==1 % traiectoria
    figure(3);
    plot(x/1e3,y/1e3,'-k','LineWidth',2);
    xlabel('x(km)'); ylabel('y(km)'); grid;
    title('The ballistic curve');
    axis equal; axis tight;
end

% Display of quantities of interest:
tf = t(i); % time of flight
b = x(i); % projectile range
h=max(y); % max altitude
tu=t(y==h); % time of ascent
tc=tf-tu; % time of descent
Q=1/2*m*(v0^2-vx(i)^2-vy(i)^2); % heat produced by friction
afis=['time of flight: ', num2str(tf),' s']; disp(afis);
afis=['projectile range: ', num2str(b/1e3),' km']; disp(afis);
afis=['max altitude: ', num2str(h/1e3),' km']; disp(afis);
afis=['time of ascent: ', num2str(tu),' s']; disp(afis);
afis=['time of descent: ', num2str(tc),' s']; disp(afis);
afis=['heat produced by friction: ', num2str(Q/1e6),' MJ']; disp(afis);

if Gd==1 % dynamical simulation
    figure(4);
    set(4,'Position',[50 50 850 600]);
    tic; simt=0; % starts the timer
    while simt<tfx
        plot(x/1e3,y/1e3,'-c'); hold on;
        xlabel('x(km)'); ylabel('y(km)'); grid;
        title('Motion simulation');
        axis equal; axis tight;
        index=abs(t-simt)==min(abs(t-simt)); % cauta cel mai apropiat t din discretizare
        plot(x(index)/1e3,y(index)/1e3,'.b','MarkerSize',10); hold off
        text(b/2/1e3,h/3/1e3,['vx=',num2str(round(round(vx(index)))),' m/s']);
        text((b/2-b/5)/1e3,h/3/1e3,['t=',num2str(round(t(index))),' s']);
        text((b/2+b/5)/1e3,h/3/1e3,['vy=',num2str(round(round(vy(index)))),' m/s']);
        pause(1e-3);
        simt=toc;
    end
end
