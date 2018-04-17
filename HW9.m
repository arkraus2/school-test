%% Homework 9
function HW9()
close all;
clc;
sv=0;

% M=[256,512,1024,2048,4096,8192]; sv=1;
M=[256,1024];
% M=256;

xl=3;
xr=7;

outputTime=[0.1,0.5,1,3];
% outputTime=1;
timeOfInterest=1;

leg={'0 s','0.1 s','0.5 s','1.0 s','3.0 s'};

CFL=0.1;

L=xr-xl;

for mind=1:length(M)
    t=0;
    dx=L/M(mind);
    x=xl-dx*3/2:dx:xr+dx*3/2;
    u=zeros(1,length(x));
    
    for ii=3:length(x)-2
        if x(ii)<4.5
            u(ii)=1/4+1/2*sin(pi/4*(x(ii)-3));
        elseif x(ii)<=5.5
            u(ii)=1/4+1/2*sin(pi/4*(x(ii)-3))+(1+cos(2*pi*x(ii)))*cos(8*pi*x(ii));
        else
            u(ii)=1/4+1/2*sin(pi/4*(x(ii)-3));
        end
    end
    u=HW9_bc(u);
    
    pl=0;
    if M(mind)==256 || M(mind)==1024
        pl=1;
        figure();
        plot(x,u)
        hold on;
    end
    
    
    while t<outputTime(end)
%     for m=1:4
        out=0;
        E=1/2*u.^2;
        alpha=HW9_alpha(u,E);
        dt=CFL*dx/max(abs(alpha));
        if t<max(outputTime)
            for ii=1:length(outputTime)
                if t<outputTime(ii) && t+dt>=outputTime(ii)
                    dt=outputTime(ii)-t;
                    out=1;
                end
            end        
        end
        
        sigma=HW9_sigma(alpha,dt,dx);
        G=HW9_G(sigma,u);
        beta=HW9_beta(u,G);
        
        p=zeros(1,length(beta));
        for ii=1:length(p)
            p(ii)=HW9_psi(alpha(ii)+beta(ii));
        end
        
        phi=zeros(1,length(p));
        for ii=1:length(phi)
            phi(ii)=(G(ii+1)+G(ii))-p(ii)*(u(ii+1)-u(ii));
        end
        
        h=zeros(1,length(phi));
        for ii=1:length(h)
            h(ii)=1/2*((E(ii+1)+E(ii))+phi(ii));
        end
        
        for ii=3:length(x)-2
            u(ii)=u(ii)-dt/dx*(h(ii)-h(ii-1));
        end
        u=HW9_bc(u);
        
        t=t+dt;
        if out==1 && pl==1
            plot(x,u)
        end
        
        if t==timeOfInterest
            switch M(mind)
                case 256
                    u256=u;
                case 512
                    u512=u;
                case 1024
                    u1024=u;
                case 2048
                    u2048=u;
                case 4096
                    u4096=u;
                case 8192
                    u8192=u;
            end
        end
        
    end
    if pl==1
        legend(leg)
        xlim([3,7])
        xlabel('x')
        ylabel('u')
        title(sprintf('M = %4.0f',M(mind)))
    end
end
if sv==1
    save('HW9_vars.mat');
end
end
%% Boundary Conditions
function u=HW9_bc(u)

u(1)=u(length(u)-3);
u(2)=u(length(u)-2);
u(length(u)-1)=u(3);
u(length(u))=u(4);

end
%% Alpha
function alpha=HW9_alpha(u,E)
len=length(u)-1;
alpha=zeros(1,len);

eps=1e-12;

for ii=1:len
    if abs(u(ii+1)-u(ii))>=eps
        alpha(ii)=(E(ii+1)-E(ii))/(u(ii+1)-u(ii));
    else
        alpha(ii)=(u(ii+1)+u(ii))/2;
    end
end

end
%% Psi
function out=HW9_psi(y)

eps=0.1;

abs_y=abs(y);
if abs_y>=eps
    out=abs_y;
else
    out=(y^2+eps^2)/(2*eps);
end

end
%% Sigma
function sigma=HW9_sigma(a,dt,dx)
len=length(a);
sigma=zeros(1,len);

for ii=1:len
    sigma(ii)=1/2*(HW9_psi(a(ii))-dt/dx*a(ii)^2);
end

end
%% G
function G=HW9_G(sigma,u)
len=length(u);
G=zeros(1,len);

for ii=3:len-2
    up=u(ii+1)-u(ii);
    um=u(ii)-u(ii-1);
    
    S=sign(up);
    
    G(ii)=S*max([0,min([sigma(ii)*abs(up),S*sigma(ii-1)*um])]);
    
end

G=HW9_bc(G);

end
%% Beta
function beta=HW9_beta(u,G)
len=length(u)-1;
beta=zeros(1,len);

eps=1e-12;

for ii=1:len
    if abs(u(ii+1)-u(ii))>=eps
        beta(ii)=(G(ii+1)-G(ii))/(u(ii+1)-u(ii));
    else
        beta(ii)=0;
    end
end

end
%% GCI Calculations
function GCI_calc()
load('HW9_vars.mat');

M=[256,512,1024,2048,4096,8192];


xl=3;
xr=7;

L=xr-xl;

for ii=1:length(M)
    dx=L/M(ii);
    x=xl-dx*3/2:dx:xr+dx*3/2;

    xpos_ind=find(x>6,1);
    xneg_ind=xpos_ind-1;
    
    switch ii
        case 1
            phi=u256;
        case 2
            phi=u512;
        case 3
            phi=u1024;
        case 4
            phi=u2048;
        case 5
            phi=u4096;
        case 6
            phi=u8192;
    end
    
    u_x6=(phi(xpos_ind)+phi(xneg_ind))/2;

    fprintf('M = %4.0f\n',M(ii))
    fprintf('u @x=6: %1.16e\n\n',u_x6)
end
end
%% GCI Analysis
function GCI_analysis()
M_1=256;
phi_1=7.9795954505809208e-01;

M_2=M_1*2;
phi_2=8.1029018817039589e-01;

M_3=M_2*2;
phi_3=8.1047378419580784e-01;

M_4=M_3*2;
phi_4=8.1055049349128838e-01;

M_5=M_4*2;
phi_5=8.1057359825292274e-01;

M_6=M_5*2;
phi_6=8.1057952153787038e-01;


% f1=phi_3; f2=phi_2; f3=phi_1; outM=M_3;
f1=phi_4; f2=phi_3; f3=phi_2; outM=M_4;
% f1=phi_5; f2=phi_4; f3=phi_3; outM=M_5;
% f1=phi_6; f2=phi_5; f3=phi_4; outM=M_6;

p=log(abs(f3-f2)/abs(f2-f1))/log(2);
% p=log((f3-f2)/(f2-f1))/log(2);

f=f1+(f1-f2)/(2^p-1);
GCI_12=1.25*abs((f1-f2)/f1)/(2^p-1)*100;
GCI_23=1.25*abs((f2-f3)/f2)/(2^p-1)*100;
approx1=GCI_12/GCI_23*2^p;

fprintf('M = %4.0f\n',outM)
fprintf('p = %1.8f\n',p)
fprintf('phi at x=0 and t=1.25: %1.16e\n',f)
fprintf('The error band is %0.5f%%\n',GCI_12)
end