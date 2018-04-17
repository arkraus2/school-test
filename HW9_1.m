%% Homework 9

close all;
clc;
sv=0;

% M=[256,512,1024,2048,4096,8192]; sv=1;
M=[256,512,1024,2048]; sv=1;
% M=[256,1024];
% M=256;

xl=3;
xr=7;

outputTime=[0.1,0.5,1,3];
% outputTime=1;
timeOfInterest=1;

leg={'0 s','0.1 s','0.5 s','1.0 s','3.0 s'};

% CFL=0.1;
CFL=0.8;

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