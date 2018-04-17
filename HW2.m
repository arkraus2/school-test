%% Homework 2
function HW2()
clear; clc;

% Bounds of the interval
xl=1;
xr=3;

k=2:10;
M=2.^k;	% Number of elements
N=M+1;  % Number of nodes

L=xr-xl;    % Length of the interval
h=L./M;     % Mesh size
h_len=length(h);

% First order backward difference
back=ones(max(N),h_len);
for jj=1:h_len
    for ii=1:N(jj)
        if ii==1
            back(ii,jj)=f_prime(xl)-1/h(jj)*(f(xl+h(jj))-f(xl));
        else
            back(ii,jj)=f_prime(xl+(ii-1)*h(jj))-1/h(jj)*(f(xl+(ii-1)*h(jj))-f(xl+(ii-2)*h(jj)));
        end
    end
end

% Second order central difference
cent=ones(max(N),h_len);
for jj=1:h_len
    for ii=1:N(jj)
        if ii==1
            cent(ii,jj)=f_prime(xl)-1/h(jj)*(f(xl+h(jj))-f(xl));
        elseif ii<N(jj)
            cent(ii,jj)=f_prime(xl+(ii-1)*h(jj))-1/(h(jj)*2)*(f(xl+(ii)*h(jj))-f(xl+(ii-2)*h(jj)));
        else
            cent(ii,jj)=f_prime(xl+(ii-1)*h(jj))-1/h(jj)*(f(xl+(ii-1)*h(jj))-f(xl+(ii-2)*h(jj)));
        end
    end
end

% Fourth order PADE
pade=ones(max(N),h_len);
for jj=1:h_len
    a=ones(N(jj),1);
    a(N(jj),1)=2;
    b=ones(N(jj),1)*4;
    b(1,1)=1;
    b(N(jj),1)=1;
    c=ones(N(jj),1);
    c(1,1)=2;
    d=ones(N(jj),1);
    for ii=1:N(jj)
        if ii==1
            d(ii,1)=1/h(jj)*(-5/2*f(xl)+2*f(xl+h(jj))+1/2*f(xl+2*h(jj)));
        elseif ii~=N(jj)
            d(ii,1)=3/h(jj)*(f(xl+ii*h(jj))-f(xl+(ii-2)*h(jj)));
        else
            d(ii,1)=1/h(jj)*(5/2*f(xr)-2*f(xr-h(jj))-1/2*f(xr-2*h(jj)));
        end
    end
    pade(1:N(jj),jj)=tri_diagonal(a,b,c,d);
    for ii=1:N(jj)
        pade(ii,jj)=f_prime(xl+(ii-1)*h(jj))-pade(ii,jj);
    end
end

% Norms
back_norm=ones(3,h_len);
cent_norm=ones(3,h_len);
pade_norm=ones(3,h_len);
for jj=1:h_len
    back_norm(1,jj)=max(abs(back(1:N(jj),jj)));
    back_norm(2,jj)=sum(abs(back(1:N(jj),jj)))/(N(jj)+1);
    back_norm(3,jj)=sqrt(sum(back(1:N(jj),jj).^2)/(N(jj)+1));
    
    cent_norm(1,jj)=max(abs(cent(1:N(jj),jj)));
    cent_norm(2,jj)=sum(abs(cent(1:N(jj),jj)))/(N(jj)+1);
    cent_norm(3,jj)=sqrt(sum(cent(1:N(jj),jj).^2)/(N(jj)+1));
    
    pade_norm(1,jj)=max(abs(pade(1:N(jj),jj)));
    pade_norm(2,jj)=sum(abs(pade(1:N(jj),jj)))/(N(jj)+1);
    pade_norm(3,jj)=sqrt(sum(pade(1:N(jj),jj).^2)/(N(jj)+1));
end

% Tables
back_array=ones(4,h_len);
cent_array=ones(4,h_len);
pade_array=ones(4,h_len);

back_array(2:4,:)=back_norm;
cent_array(2:4,:)=cent_norm;
pade_array(2:4,:)=pade_norm;

for jj=1:h_len
    back_array(1,jj)=M(jj);
    cent_array(1,jj)=M(jj);
    pade_array(1,jj)=M(jj);
end

back_array=back_array';
cent_array=cent_array';
pade_array=pade_array';

for jj=5:7
    back_array(1,jj)=NaN;
    cent_array(1,jj)=NaN;
    pade_array(1,jj)=NaN;
end

for ii=2:h_len
    back_array(ii,5)=log(back_array(ii-1,2)/back_array(ii,2))/log(h(ii-1)/h(ii));
    back_array(ii,6)=log(back_array(ii-1,3)/back_array(ii,3))/log(h(ii-1)/h(ii));
    back_array(ii,7)=log(back_array(ii-1,4)/back_array(ii,4))/log(h(ii-1)/h(ii));
    
    cent_array(ii,5)=log(cent_array(ii-1,2)/cent_array(ii,2))/log(h(ii-1)/h(ii));
    cent_array(ii,6)=log(cent_array(ii-1,3)/cent_array(ii,3))/log(h(ii-1)/h(ii));
    cent_array(ii,7)=log(cent_array(ii-1,4)/cent_array(ii,4))/log(h(ii-1)/h(ii));
    
    pade_array(ii,5)=log(pade_array(ii-1,2)/pade_array(ii,2))/log(h(ii-1)/h(ii));
    pade_array(ii,6)=log(pade_array(ii-1,3)/pade_array(ii,3))/log(h(ii-1)/h(ii));
    pade_array(ii,7)=log(pade_array(ii-1,4)/pade_array(ii,4))/log(h(ii-1)/h(ii));
end

back_table=array2table(back_array,'VariableNames',{'M','Loo_Norm','L1_Norm','L2_Norm',...
    'Order_Loo','Order_L1','Order_L2'});
cent_table=array2table(cent_array,'VariableNames',{'M','Loo_Norm','L1_Norm','L2_Norm',...
    'Order_Loo','Order_L1','Order_L2'});
pade_table=array2table(pade_array,'VariableNames',{'M','Loo_Norm','L1_Norm','L2_Norm',...
    'Order_Loo','Order_L1','Order_L2'});

disp('First Order Backward Difference')
disp(back_table)
disp('Second Order Central Difference')
disp(cent_table)
disp('Fourth Order Pade')
disp(pade_table)

% Plot 
% k(3)=4
ind=3;
x=xl:h(ind):xr;
b=back(1:N(ind),ind);
c=cent(1:N(ind),ind);
p=pade(1:N(ind),ind);
figure()
plot(x,b,'-*')
title('First Order Backward Difference @ k=4');
xlabel('x');
ylabel('error');
figure()
plot(x,c,'-*')
title('Second Order Central Difference @ k=4');
xlabel('x');
ylabel('error');
figure()
plot(x,p,'-*')
title('Fourth Order PADE @ k=4');
xlabel('x');
ylabel('error');

Bonus();
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
%% Given function f
function out=f(x)
out=(1+ 2 * x.^2 .* cos(x))./(x.^(2.4));
end
%% Analytical derivative of the given function f
function out=f_prime(x)
out=-1*(2.4*x.^(-3.4)+0.8*x.^(-1.4).*cos(x)+2*x.^(-0.4).*sin(x));
end
%% Bonus Problem
function Bonus()
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

cent_array=M;
cent_array(1,2:4)=cent_norm';

cent_table=array2table(cent_array,'VariableNames',{'M','Loo_Norm','L1_Norm','L2_Norm'});
disp('Uneven Spacing Norm')
disp(cent_table)

hi=zeros(N,1);
for ii=1:M
    hi(ii+1,1)=x(ii+1)-x(ii);
end

figure()
plot(x,hi,'*')
title('Unequal Mesh Spacing');
xlabel('x');
ylabel('Mesh Spacing');

figure()
plot(x,cent,'-*')
title('Error With Unequal Mesh Spacing');
xlabel('x');
ylabel('error');
end