clc;

M=256*2*2*2*2*2;

a=@(x,t) 2.3+1.7*sin(2*pi*x)+1.5*sin(5*pi*t);

xl=-1;
xr=3;

% outputTime=[0.25,0.5,1,1.25,1.5,2.1];
outputTime=1.25;

L=xr-xl;
h=L/M;
x=xl-h*5/2:h:xr+h*5/2;

t=0;

phi=zeros(1,M+6);
phi=HW8_bc(phi,x,xl,t);


alpha10=1;
alpha20=-3/4;
alpha21=1/4;
alpha30=-1/12;
alpha31=-1/12;
alpha32=2/3;

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
    for ii=4:M+3
        if a(x(ii),t)>=0
            pos=1;
        else
            pos=-1;
        end
        
        phi1(ii)=phi(ii)-alpha10*(a(x(ii),t)*dt*WENO5(phi,ii,pos,h));
    end
    phi1=HW8_bc(phi1,x,xl,t);
    
    % Step 2
    for ii=4:M+3
        if a(x(ii),t)>=0
            pos=1;
        else
            pos=-1;
        end
        
        phi2(ii)=phi1(ii)-alpha20*(a(x(ii),t)*dt*WENO5(phi,ii,pos,h))-alpha21*(a(x(ii),t)*dt*WENO5(phi1,ii,pos,h));
    end
    phi2=HW8_bc(phi2,x,xl,t);
    
    % Step 3
    for ii=4:M+3
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
    
%     phi=HW8_bc(phi,x,xl,t);
    
%     if out==1
%         figure();
%         plot(x,phi);
%         xlim([-0.1,1.1]);
% %         ylim([0,1]);
%         xlabel('x')
%         ylabel('\phi')
%         title(sprintf('t = %1.2f',t))
%     end
end

xpos_ind=find(x>0,1);
xneg_ind=xpos_ind-1;

% phi_x0=phi(xpos_ind)-(x(xpos_ind)*(phi(xpos_ind)-phi(xneg_ind)))/(x(xpos_ind)-x(xneg_ind));
phi_x0=(phi(xpos_ind)+phi(xneg_ind))/2;

fprintf('M = %4.0f\n',M)
fprintf('phi @x=0: %1.16e\n',phi_x0)

x0=(x(xpos_ind)+x(xneg_ind))/2;
fprintf('x=0: %1.16e\n',x0)

fprintf('t: %1.16f\n',t)