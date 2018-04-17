clear; clc;
% Given Function
f = @(x) (1+ 2 * x.^2 .* cos(x))./(x.^(2.4));

% Analytical Derivative of Given Function
f_prime = @(x) -1*(2.4*x.^(-3.4)+0.8*x.^(-1.4).*cos(x)+2*x.^(-0.4).*sin(x));

% Bounds of the interval
xl=1;
xr=3;

M=2^4;
N=M+1;

L=xr-xl;    % Length of the interval
h=L/M;     % Mesh size

% Set the uneven grid spacing
a=.3;
b=1-4*a;
c=3*a;

xi=(xl:h:xr)';
x=a*xi.^2+b*xi+c;

cent=ones(N,1);
cent(1,1)=f_prime(x(1))-1/h*(f(x(2))-f(x(1)))*(1/(2*a*sqrt(x(1)/a+(b-2*c)/(2*a))));
cent(N,1)=f_prime(x(N))-1/h*(f(x(N))-f(x(N-1)))*(1/(2*a*sqrt(x(1)/a+(b-2*c)/(2*a))));
for ii=2:M
    cent(ii,1)=f_prime(x(ii))-1/(h*2)*(f(x(ii+1))-f(x(ii-1)))*(1/(2*a*sqrt(x(1)/a+(b-2*c)/(2*a))));
end

cent_norm(1,1)=max(abs(cent(:,1)));
cent_norm(2,1)=sum(abs(cent(:,1)))/(N+1);
cent_norm(3,1)=sqrt(sum(cent(:,1).^2)/(N+1));