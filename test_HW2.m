clear; clc;
% Given Function
f = @(x) (1+ 2 * x.^2 .* cos(x))./(x.^(2.4));

% Analytical Derivative of Given Function
f_prime = @(x) -1*(2.4*x.^(-3.4)+0.8*x.^(-1.4).*cos(x)+2*x.^(-0.4).*sin(x));

% Bounds of the interval
x0=1;
x3=3;

M=2^4;  % Number of elements
N=M+1;  % Number of nodes

L=x3-x0;    % Length of the interval
h=L/M;     % Mesh size

% Set the uneven grid spacing
% Using 3 lines of the form x=m*xi+b
% xi=(x-b)/m -> xi'=1/m
m1=.25;
z1=.5;

m3=.25;
z2=.5;

b1=x0-m1*x0;
b3=x3-m3*x3;

m2=(x3-x0-m1*z1-m3*z2)/(x3-x0-z2-z1);
b2=m1*z1+x0-m2*(x0+z1);

xi=(x0:h:x3)';
x=ones(length(xi),1);
for ii=1:length(xi)
    if xi(ii)<=x0+z1
        x(ii)=m1*xi(ii)+b1;
    elseif xi(ii)<=x3-z2
        x(ii)=m2*xi(ii)+b2;
    else
        x(ii)=m3*xi(ii)+b3;
    end
end

% Using central difference
cent=ones(N,1);
cent(1,1)=f_prime(x(1))-1/h*(f(x(2))-f(x(1)))*(1/m1);   % Forward difference at left bound
cent(N,1)=f_prime(x(N))-1/h*(f(x(N))-f(x(N-1)))*(1/m3); % Backward difference at right bound
for ii=2:M
    if xi(ii)<=x0+z1
        if xi(ii+1)>=x0+z1
            % Backward difference at discontinuous derivative
            cent(ii,1)=f_prime(x(ii))-1/(h)*(f(x(ii))-f(x(ii-1)))*(1/m1);
        else
            cent(ii,1)=f_prime(x(ii))-1/(h*2)*(f(x(ii+1))-f(x(ii-1)))*(1/m1);
        end
    elseif xi(ii)<=x3-z2
        if xi(ii+1)>=x3-z2
            % Backward difference at discontinuous derivative
            cent(ii,1)=f_prime(x(ii))-1/(h)*(f(x(ii))-f(x(ii-1)))*(1/m2);
        else
            cent(ii,1)=f_prime(x(ii))-1/(h*2)*(f(x(ii+1))-f(x(ii-1)))*(1/m2);
        end
    else
        cent(ii,1)=f_prime(x(ii))-1/(h*2)*(f(x(ii+1))-f(x(ii-1)))*(1/m3);
    end
end

cent_norm(1,1)=max(abs(cent(:,1)));
cent_norm(2,1)=sum(abs(cent(:,1)))/(N+1);
cent_norm(3,1)=sqrt(sum(cent(:,1).^2)/(N+1));

hi=zeros(N,1);
for ii=1:M
    hi(ii+1,1)=x(ii+1)-x(ii);
end
plot(x,hi,'*')