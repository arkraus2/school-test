%% Homework 8
function HW8()
clc;

M=[256,8192];

a=@(x,t) 2.3+1.7*sin(2*pi*x)+1.5*sin(5*pi*t);
% a=@(x,t) 2.3-1.7*sin(2*pi*x)-1.5*sin(5*pi*t);

xl=-1;
xr=3;

outputTime=[0.25,0.5,1,1.25,1.5,2.1];

L=xr-xl;

alpha10=1;
alpha20=-3/4;
alpha21=1/4;
alpha30=-1/12;
alpha31=-1/12;
alpha32=2/3;

for mind=1:2
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
        

        if out==1
            figure();
            plot(x,phi);
            ylim([-0.1,1.1]);
            xlim([-1,3]);
            xlabel('x')
            ylabel('\phi')
            title(sprintf('t = %1.2f, M = %4.0f',t,M(mind)))
        end
    end
end

GCI_calc();
end
%% Boundary Conditions
function phi=HW8_bc(phi,x,xl,t)

phi(1)=phi(6)-((x(6)-x(1))*(phi(6)-Left(t)))/(x(6)-xl);
phi(2)=phi(5)-((x(5)-x(2))*(phi(5)-Left(t)))/(x(5)-xl);
phi(3)=phi(4)-((x(4)-x(3))*(phi(4)-Left(t)))/(x(4)-xl);

phi(end)=phi(end-5);
phi(end-1)=phi(end-4);
phi(end-2)=phi(end-3);
end
%% Left Boundary at t
function phi=Left(t)

if t<=1/4
    phi=1;
elseif t<=1/2
    phi=0;
elseif t<=1
    phi=(1-cos(4*pi*t))/2;
else
    phi=0;
end

end
%% WENO5
function out=WENO5(phi,ind,pos,h)

out=1/(12*h)*(-dP(phi,ind-2)+7*dP(phi,ind-1)+7*dP(phi,ind)-dP(phi,ind+1));

if pos>0
    a=dMdP(phi,ind-2)/h;
    b=dMdP(phi,ind-1)/h;
    c=dMdP(phi,ind)/h;
    d=dMdP(phi,ind+1)/h;
    out=out-psiWENO(a,b,c,d);
else
    a=dMdP(phi,ind+2)/h;
    b=dMdP(phi,ind+1)/h;
    c=dMdP(phi,ind)/h;
    d=dMdP(phi,ind-1)/h;
    out=out+psiWENO(a,b,c,d);
end

end
%% Delta+
function dp=dP(phi,ind)
dp=phi(ind+1)-phi(ind);
end
%% Delta- Delta+
function dmdp=dMdP(phi,ind)
dmdp=phi(ind+1)-2*phi(ind)+phi(ind-1);
end
%% psiWENO
function psi=psiWENO(a,b,c,d)

IS_0=13*(a-b)^2+3*(a-3*b)^2;
IS_1=13*(b-c)^2+3*(b+c)^2;
IS_2=13*(c-d)^2+3*(3*c-d)^2;

epsilon=1e-6;

alpha0=1/(epsilon+IS_0)^2;
alpha1=6/(epsilon+IS_1)^2;
alpha2=3/(epsilon+IS_2)^2;

w0=alpha0/(alpha0+alpha1+alpha2);
w2=alpha2/(alpha0+alpha1+alpha2);

psi=1/3*w0*(a-2*b+c)+1/6*(w2-1/2)*(b-2*c+d);

end
%% GCI Analysis
function GCI_A()
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
    phi3=HW8_bc(phi3,x,xl,t);
    
    phi=phi3;
    t=t+dt;
    
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

phi_x0=phi(xpos_ind)-(x(xpos_ind)*(phi(xpos_ind)-phi(xneg_ind)))/(x(xpos_ind)-x(xneg_ind));
fprintf('phi @x=0: %1.16e\n',phi_x0)
end
%% GCI Calculation
function GCI_calc()
M_1=256;
phi_1=4.6603455097273294e-01;

M_2=M_1*2;
phi_2=4.6967762130439206e-01;

M_3=M_2*2;
phi_3=4.7138255039699994e-01;

M_4=M_3*2;
phi_4=4.7230196350605191e-01;

M_5=M_4*2;
phi_5=4.7274427087938542e-01;

M_6=M_5*2;
phi_6=4.7296689307598505e-01;


f1=phi_6; f2=phi_5; f3=phi_4;

p=log(abs(f3-f2)/abs(f2-f1))/log(2);
f=f1+(f1-f2)/(2^p-1);
GCI_12=1.25*abs((f1-f2)/f1)/(2^p-1)*100;
GCI_23=1.25*abs((f2-f3)/f2)/(2^p-1)*100;
approx1=GCI_12/GCI_23*2^p;

fprintf('phi at x=0 and t=1.25: %1.16e\n',f)
fprintf('The error band is %0.5f%%\n',GCI_12)
end