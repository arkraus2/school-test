%% Homework 8

clc;

M=[256,512,1024,2048,4096,8192];
% M=256;

a=@(x,t) 2.3+1.7*sin(2*pi*x)+1.5*sin(5*pi*t);
% a=@(x,t) 2.3-1.7*sin(2*pi*x)-1.5*sin(5*pi*t);

xl=-1;
xr=3;

% outputTime=[0.25,0.5,1,1.25,1.5,2.1];
outputTime=1.25;

L=xr-xl;

alpha10=1;
alpha20=-3/4;
alpha21=1/4;
alpha30=-1/12;
alpha31=-1/12;
alpha32=2/3;

for mind=1:length(M)
    h=L/M(mind);
    x=xl-h*5/2:h:xr+h*5/2;

    t=0;

    phi=zeros(1,M(mind)+6);
    phi=HW8_bc(phi,x,xl,t);
    while t<outputTime(end)
    % for m=1:4

        out=0;
        dt=h/max(a(x,t))*0.8;
        if t<max(outputTime)
            for ii=1:length(outputTime)
                if t<outputTime(ii) && t+dt>=outputTime(ii)
                    dt=outputTime(ii)-t;
                    out=1;
                end
            end        
        end

        phi1=phi;
        phi2=phi;
        phi3=phi;

        % Step 1
        for ii=4:M(mind)+3
            if a(x(ii),t)>=0
                pos=1;
            else
                pos=-1;
            end

            phi1(ii)=phi(ii)-alpha10*(a(x(ii),t)*dt*WENO5(phi,ii,pos,h));
        end
        phi1=HW8_bc(phi1,x,xl,t);

        % Step 2
        for ii=4:M(mind)+3
            if a(x(ii),t)>=0
                pos=1;
            else
                pos=-1;
            end

            phi2(ii)=phi1(ii)-alpha20*(a(x(ii),t)*dt*WENO5(phi,ii,pos,h))-alpha21*(a(x(ii),t)*dt*WENO5(phi1,ii,pos,h));
        end
        phi2=HW8_bc(phi2,x,xl,t);

        % Step 3
        for ii=4:M(mind)+3
            if a(x(ii),t)>=0
                pos=1;
            else
                pos=-1;
            end

            phi3(ii)=phi2(ii)-alpha30*(a(x(ii),t)*dt*WENO5(phi,ii,pos,h))-alpha31*(a(x(ii),t)*dt*WENO5(phi1,ii,pos,h))-...
                alpha32*(a(x(ii),t)*dt*WENO5(phi2,ii,pos,h));
        end
        
        t=t+dt;
        phi3=HW8_bc(phi3,x,xl,t);

        phi=phi3;
        
%         phi=HW8_bc(phi,x,xl,t);

        if out==1
            figure();
            plot(x,phi);
            ylim([-0.1,1.1]);
            xlim([-1,3]);
            xlabel('x')
            ylabel('\phi')
            title(sprintf('t = %1.2f, M = %4.0f',t,M(mind)))
        end
        if t==1.25
            switch mind
                case 1
                    phi256=phi;
                case 2
                    phi512=phi;
                case 3
                    phi1024=phi;
                case 4
                    phi2048=phi;
                case 5
                    phi4096=phi;
                case 6
                    phi8192=phi;
            end
        end
    end
end