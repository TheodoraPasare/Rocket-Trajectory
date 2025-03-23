clc; clear; close all;
% toate marimile sunt exprimate in unitati SI

% Selectori pentru grafice (1-on; 0-off):
Gv = 1; % componentele vitezei ca functii de timp
Gc = 1; % legile de miscare
Gt = 1; % traiectoria (curba balistica)
Gd = 1; % reprezentarea dinamica
if Gd==1, Gv=0; Gc=0; Gt=0; end

% Parametrii fizici:
g = 9.80665; % acceleratia gravitationala
ro = 7850; % densitatea otelului
r = 0.13; % raza proiectilului
m = 4/3 * pi * r^3 * ro; % masa proiectilului
G = m*g; % greutatea proiectilului

% Conditii initiale:
v0 = 1100; % viteza initiala
alpha0 = 43; % unghiul de lansare
eta = 1.81 * 1e-5; % coeficientul de vascozitate
b1 = 6 * pi * eta * r; % coeficient termen liniar
c = 0.469; % coeficient de forma
ro0 = 1.22; % densitatea aerului
b2 = c * 4 * pi * r^2 * ro0/2; % coeficient termen patratic

% Definirea intervalului de timp de inters
t0 = 0; tf = 2 * v0 /g*sind(alpha0);
N = 1500; % numarul momentelor de timp
t = linspace(t0,tf,N); dt = t(2) - t(1);

% Prealocare si valori de inceput:
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
t = t(1:i); vx=vx(1:i); vy=vy(1:i); x=x(1:i); y=y(1:i); % eliminarea valorilor in surplus
if Gv==1 % legile de viteza
    figure(1);
    plot(t,vx,'-r',t,vy,'-b');
    xlabel('t(s)'); ylabel('v(m/s)'); grid;
    title('Componentele vitezei ca functii de timp');
    legend('vx','vy');
end
if Gc==1 % legile de miscare
    figure(2);
    plot(t,x/1e3,'-r',t,y/1e3,'-b');
    xlabel('t(s)'); ylabel('coord(km)'); grid;
    title('Coordonatele ca functii de timp');
    legend('x','y','Location','northwest');
end
if Gt==1 % traiectoria
    figure(3);
    plot(x/1e3,y/1e3,'-k','LineWidth',2);
    xlabel('x(km)'); ylabel('y(km)'); grid;
    title('Curba balistica');
    axis equal; axis tight;
end

% Afisarea unor marimi de interes:
tf = t(i); % timpul de zbor
b = x(i); % bataia
h=max(y); % altitudinea maxima
tu=t(y==h); % timpul de urcare
tc=tf-tu; % timpul de coborare
Q=1/2*m*(v0^2-vx(i)^2-vy(i)^2); % caldura produsa prin frecare
afis=['Timpul de zbor: ', num2str(tf),' s']; disp(afis);
afis=['Bataia proiectilului: ', num2str(b/1e3),' km']; disp(afis);
afis=['Altitudinea maxima: ', num2str(h/1e3),' km']; disp(afis);
afis=['Timpul de urcare: ', num2str(tu),' s']; disp(afis);
afis=['Timpul de coborare: ', num2str(tc),' s']; disp(afis);
afis=['Caldura produsa: ', num2str(Q/1e6),' MJ']; disp(afis);

if Gd==1 % simularea dinamica
    figure(4);
    set(4,'Position',[50 50 850 600]);
    tic; simt=0; % porneste cronometrul si retine timpul initial
    while simt<tfx
        plot(x/1e3,y/1e3,'-c'); hold on;
        xlabel('x(km)'); ylabel('y(km)'); grid;
        title('Simularea miscarii');
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
