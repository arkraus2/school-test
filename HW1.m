%% Homework 1 Problem 2
function HW1()
clear; clc;
% Given Function
f = @(x) (1+ 2 * x.^2 .* cos(x))./(x.^(2.4));

% Analytical Derivative of Given Function
f_prime = @(x) -1*(2.4*x.^(-3.4)+0.8*x.^(-1.4).*cos(x)+2*x.^(-0.4).*sin(x));

% k and h
k=-3:1:25;
h=2.^(-k);

% Value to find the derivative at: x = 33.3
x=33.3;

% first-order forward difference
forward_1 = (f(x+h) - f(x))./h;

% second-order central difference
central_2 = (f(x+h) - f(x-h)) ./ (2*h);

% sixth-order central difference
central_6 = (-f(x-3*h)+9*f(x-2*h)-45*f(x-h)+45*f(x+h)-9*f(x+2*h)+f(x+3*h))./(60*h);

%------------------------------------------------------------------------
% Error calculations
forward_1_e=abs(f_prime(x)-forward_1);
central_2_e=abs(f_prime(x)-central_2);
central_6_e=abs(f_prime(x)-central_6);

% Log 10 of h and errors
h_log=log10(h);
forward_1_e_log = log10(forward_1_e);
central_2_e_log = log10(central_2_e);
central_6_e_log = log10(central_6_e);

[m_f1_log,b_f1_log] = lin_reg(h_log(3:28),forward_1_e_log(3:28));
[m_c2_log,b_c2_log] = lin_reg(h_log(3:21),central_2_e_log(3:21));
[m_c6_log,b_c6_log] = lin_reg(h_log(3:10),central_6_e_log(3:10));

% fprintf('f1_slope = %g\n',m_f1_log)
% fprintf('c2_slope = %g\n',m_c2_log)
% fprintf('c6_slope = %g\n',m_c6_log)

f1_log=m_f1_log.*h_log(3:28)+b_f1_log;
c2_log=m_c2_log.*h_log(3:21)+b_c2_log;
c6_log=m_c6_log.*h_log(3:10)+b_c6_log;

f1=zeros(1,length(f1_log));
c2=zeros(1,length(c2_log));
c6=zeros(1,length(c6_log));

for ii=1:length(f1_log)
    f1(ii)=10^f1_log(ii);
end

for ii=1:length(c2_log)
    c2(ii)=10^c2_log(ii);
end

for ii=1:length(c6_log)
    c6(ii)=10^c6_log(ii);
end

loglog(h,forward_1_e,'-xk',h,central_2_e,'-xb',h,central_6_e,'-xr');
hold on;
loglog(h(3:28),f1,'-*r',h(3:21),c2,'-*k',h(3:10),c6,'-*b');
xlabel('h');
ylabel('|e|');
title('|e| vs. h');
legend('first-order forward difference','second-order central difference',...
    'sixth-order central difference',...
    sprintf('slope = %g',m_f1_log),sprintf('slope = %g',m_c2_log),...
    sprintf('slope = %g',m_c6_log),'Location','northwest');
end
%% Linear regression
function [m,b] = lin_reg(x,y)
size=length(x);

x_avg=mean(x);
y_avg=mean(y);

m_num=0;
m_den=0;

for jj=1:size
    x_comp=(x(jj)-x_avg);
    m_num=m_num+x_comp*(y(jj)-y_avg);
    m_den=m_den+x_comp*x_comp;
end

m=m_num/m_den;
b=y_avg-m*x_avg;
end