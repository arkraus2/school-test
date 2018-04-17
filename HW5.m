function d=HW5()
% clc;

M=8;
sz=M+1;

alpha=0.1;

xl=-1;
xr=1;

max_t=6.25;

L=xr-xl;
h=L/M;

x=(xl:h:xr)';

dt=2*h^2/alpha;

mult=alpha*dt/(2*h^2);

T0=2*ones(sz,1);

a=ones(sz,1)*-alpha*dt/(2*h^2);
b=ones(sz,1)*(1+alpha*dt/h^2);
b(sz)=1+alpha*dt/(2*h^2);
c=a;
d=ones(sz,1);

t=0;
for ii=1:sz
    if ii==1
        d(ii)=(1-2*mult)*T0(ii)+mult*T0(ii+1)+mult*(2-sin(3*pi/2*(t+dt)));
    elseif ii==sz
        d(ii)=mult*T0(ii-1)+(1-2*mult)*T0(ii);
    else
        d(ii)=mult*T0(ii-1)+(1-2*mult)*T0(ii)+mult*T0(ii+1);
    end
end

q0=q(x,0,alpha);
q1=q(x,dt,alpha);
q_avg=1/2*(q0+q1);

d=d+q_avg;

end
%% Tri-Diagonal Solver
function d=tri_diagonal(a,b,c,d)
P=length(d);

for ii = 2:P
    b(ii)=b(ii) - c(ii-1)*a(ii)/b(ii-1);
    d(ii)=d(ii) - d(ii-1)*a(ii)/b(ii-1);
end

d(P) = d(P)/b(P);
for jj = P-1:-1:1
    d(jj) = (d(jj) - c(jj)*d(jj+1))/b(jj);
end

end
%% Source Term
function out=q(x,t,a)
if t==0
    out=x*0;
else
    out=3*pi/(2*sqrt(t))*sin(pi/2*sqrt(t))*cos(pi/2*sqrt(t))*(x.^3-x.^2-x+1)+...
        3*pi/2*cos(3*pi/2*t)*sin(pi/2*x)+a*((sin(pi/2*sqrt(t))).^2.*(6-18*x)+pi^2/4*sin(3*pi/2*t)*sin(pi/2*x));
end
end