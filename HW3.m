%% Homework 3
function HW3()
clc;
% tic;
% p1();
p2();
% p3();
% p4();
% toc
end
%% Problem 1
function p1()
M=256;	%Number of elements
xl=-1;  %Left bound
xr=1;   %Right bound
lbc=0;  %Left Dirichlet Boundary Condition
rbc=0;  %Right Dirichlet Boundary Condition

N=M+1;	%Number of nodes
L=xr-xl;    %Length of interval
h=L/M;  %Mesh spacing

x=xl:h:xr;  %x-axis grid

% Initial Guess
phi0=1/2*sin(pi*x)-1/10*sin(32*pi*x);
phi0(1)=lbc;
phi0(N)=rbc;

% Right Hand Side
f1=cos(pi*x).*(-(pi^2/4)*x.^3-2*pi^2*x.^2+3/2*x+4)-...
    sin(pi*x).*(3/2*pi*x.^2+8*pi*x);

rho=cos(pi/M);
w=2/(1+sqrt(1-rho^2));

fprintf('Optimum over-relaxation factor for problem 1: %c = %.16f\n',969,w);

num_it=1000;
ind=num_it+1;

phi_pj=zeros(ind,length(phi0));
phi_gs=zeros(ind,length(phi0));
phi_sor=zeros(ind,length(phi0));

phi_pj(1,:)=phi0;
phi_gs(1,:)=phi0;
phi_sor(1,:)=phi0;

for ii=2:ind
    phi_pj(ii,:)=point_Jacobi(phi_pj(ii-1,:),f1,h,1);
    phi_gs(ii,:)=Gauss_Seidel_1D(phi_gs(ii-1,:),f1,h,1);
    phi_sor(ii,:)=SOR_1D(phi_sor(ii-1,:),f1,h,w,1);
end

figure();
plot(x,phi_pj(1,:),'k',x,phi_pj(201,:),x,phi_pj(401,:),x,phi_pj(601,:),x,phi_pj(801,:),x,phi_pj(1001,:))
title('Point Jacobi');
ylim([-1,2.5]);
xlabel('x');
ylabel('\phi');
legend('Initial Guess','200 iterations','400 iterations','600 iterations',...
    '800 iterations','1000 iterations','Location', 'northeastoutside');

figure();
plot(x,phi_gs(1,:),'k',x,phi_gs(201,:),x,phi_gs(401,:),x,phi_gs(601,:),x,phi_gs(801,:),x,phi_gs(1001,:))
title('Gauss-Seidel');
ylim([-1,2.5]);
xlabel('x');
ylabel('\phi');
legend('Initial Guess','200 iterations','400 iterations','600 iterations',...
    '800 iterations','1000 iterations','Location', 'northeastoutside');

figure();
plot(x,phi_sor(1,:),'k',x,phi_sor(201,:),x,phi_sor(401,:),x,phi_sor(601,:),x,phi_sor(801,:),x,phi_sor(1001,:))
title(sprintf('SOR with \\omega = %g',w));
ylim([-1,2.5]);
xlabel('x');
ylabel('\phi');
legend('Initial Guess','200 iterations','400 iterations','600 iterations',...
    '800 iterations','1000 iterations','Location', 'northeastoutside');

its=401;

r_pj=zeros(its,length(phi0));
r_gs=zeros(its,length(phi0));
r_sor=zeros(its,length(phi0));

norm_pj=zeros(its,1);
norm_gs=zeros(its,1);
norm_sor=zeros(its,1);

%Residuals
for jj=1:its
    r_pj(jj,:)=residual_1D(phi_pj(jj,:),f1,h);
    r_gs(jj,:)=residual_1D(phi_gs(jj,:),f1,h);
    r_sor(jj,:)=residual_1D(phi_sor(jj,:),f1,h);
end

%Norms of the Residuals
iterations=0:(its-1);
for kk=1:its
    norm_pj(kk,1)=max(abs(r_pj(kk,:)));
    norm_gs(kk,1)=max(abs(r_gs(kk,:)));
    norm_sor(kk,1)=max(abs(r_sor(kk,:)));
end

figure()
semilogy(iterations,norm_pj,iterations,norm_gs,iterations,norm_sor)
title('Problem 1');
xlabel('Iterations');
ylabel('Residual Norm');
legend('Point Jacobi','Gauss-Seidel','SOR');
end
%% Problem 2
function p2()
m=80;    %Number of elements in x
n=m;    %Number of elements in y
M=m+1;	%Number of points in x
N=n+1;    %Number of points in y

xl=-1;  %Left bound
xr=1;   %Right bound
lbc=0;  %Left Dirichlet Boundary Condition
rbc=0;  %Right Dirichlet Boundary Condition

yb=-1;  %Bottom bound
yt=1;   %Top bound
bbc=0;  %Bottom Dirichlet Boundary Condition
tbc=0;  %Top Dirichlet Boundary Condition

Lx=xr-xl;   %Length of x interval
hx=Lx/(M-1);  %Mesh spacing in x

x=xl:hx:xr;  %x-axis grid

Ly=yt-yb;   %Length of y interval
hy=Ly/(N-1);  %Mesh spacing in y

y=yb:hy:yt;  %y-axis grid

% Initial Guess
phi0=zeros(M,N);
for ii=1:M
    for jj=1:N
        if ii==1
            phi0(ii,jj)=lbc;
        elseif ii==M
            phi0(ii,jj)=rbc;
        elseif jj==1
            phi0(ii,jj)=bbc;
        elseif jj==N
            phi0(ii,jj)=tbc;
        else
            phi0(ii,jj)=2/5*sin(8*pi*x(ii))*sin(2*pi*y(jj));
        end
    end
end

% Right Hand Side
f2=zeros(M,N);
for ii=2:M-1
    for jj=2:N-1
        f2(ii,jj)=sqrt(y(jj)+1)*sin(pi/2*(y(jj)+1))*...
            (-pi^2*x(ii)^3*cos(pi/2*x(ii))-6*pi*x(ii)^2*sin(pi/2*x(ii))+...
            12*x(ii)*cos(pi/2*x(ii)))-...
            x(ii)^3*cos(pi/2*x(ii))*(sin(pi/2*(y(jj)+1))/(2*(y(jj)+1)^(3/2))-...
            pi*cos(pi/2*(y(jj)+1))/sqrt(y(jj)+1));
    end
end

rho=1/2*(cos(pi/(M-1))+cos(pi/(N-1)));
w=2/(1+sqrt(1-rho^2));

fprintf('Optimum over-relaxation factor for problem 2: %c = %.16f\n',969,w);


num_it=400;
its=num_it+1;

phi_gs=zeros(M,N,its);
phi_sor=zeros(M,N,its);

phi_gs(:,:,1)=phi0;
phi_sor(:,:,1)=phi0;

for ii=2:its
    phi_gs(:,:,ii)=Gauss_Seidel_2D(phi_gs(:,:,ii-1),f2,hx,1);
    phi_sor(:,:,ii)=SOR_2D(phi_sor(:,:,ii-1),f2,hx,w,1);
end

r_gs=zeros(M,N,its);
r_sor=zeros(M,N,its);



for jj=1:its
    r_gs(:,:,jj)=residual_2D(phi_gs(:,:,jj),f2,hx);
    r_sor(:,:,jj)=residual_2D(phi_sor(:,:,jj),f2,hx);
end

%Norms of the Residuals
r_it=400;
iterations=0:r_it;
norm_gs=zeros(r_it+1,1);
norm_sor=zeros(r_it+1,1);

for kk=1:r_it+1
    norm_gs(kk,1)=max(max(abs(r_gs(:,:,kk))));
    norm_sor(kk,1)=max(max(abs(r_sor(:,:,kk))));
end

pl_x=ones(length(y),length(x)).*x;
pl_y=trans(ones(length(y),length(x)).*y);

figure();
surf(pl_x,pl_y,phi_gs(:,:,1))
title('Gauss-Seidel Initial Iteration');
zlim([-0.4,0.4])
xlabel('x');
ylabel('y');
zlabel('\phi');
figure();
surf(pl_x,pl_y,phi_gs(:,:,11))
title('Gauss-Seidel 10th Iteration');
zlim([-0.4,0.4])
xlabel('x');
ylabel('y');
zlabel('\phi');
figure();
surf(pl_x,pl_y,phi_gs(:,:,21))
title('Gauss-Seidel 20th Iteration');
zlim([-0.4,0.4])
xlabel('x');
ylabel('y');
zlabel('\phi');
figure();
surf(pl_x,pl_y,phi_gs(:,:,101))
title('Gauss-Seidel 100th Iteration');
zlim([-0.4,0.4])
xlabel('x');
ylabel('y');
zlabel('\phi');
figure();
surf(pl_x,pl_y,phi_gs(:,:,401))
title('Gauss-Seidel 400th Iteration');
zlim([-0.4,0.4])
xlabel('x');
ylabel('y');
zlabel('\phi');

figure();
surf(pl_x,pl_y,phi_sor(:,:,1))
title('Succesive Over-Relaxation Initial Iteration');
zlim([-0.4,0.4])
xlabel('x');
ylabel('y');
zlabel('\phi');
figure();
surf(pl_x,pl_y,phi_sor(:,:,11))
title('Succesive Over-Relaxation 10th Iteration');
zlim([-0.4,0.4])
xlabel('x');
ylabel('y');
zlabel('\phi');
figure();
surf(pl_x,pl_y,phi_sor(:,:,21))
title('Succesive Over-Relaxation 20th Iteration');
zlim([-0.4,0.4])
xlabel('x');
ylabel('y');
zlabel('\phi');
figure();
surf(pl_x,pl_y,phi_sor(:,:,101))
title('Succesive Over-Relaxation 100th Iteration');
zlim([-0.4,0.4])
xlabel('x');
ylabel('y');
zlabel('\phi');
figure();
surf(pl_x,pl_y,phi_sor(:,:,401))
title('Succesive Over-Relaxation 400th Iteration');
zlim([-0.4,0.4])
xlabel('x');
ylabel('y');
zlabel('\phi');

figure()
semilogy(iterations,norm_gs,iterations,norm_sor)
title('Problem 2');
xlabel('Iterations');
ylabel('Residual Norm');
legend('Gauss-Seidel','SOR');
end
%% Problem 3
function p3()
num_it=1000;

M=2.^[7,8,9];   %Number of elements
N=M+1;          %Number of nodes

xl=-1;  %Left bound
xr=1;   %Right bound
lbc=0;  %Left Dirichlet Boundary Condition
rbc=0;  %Right Dirichlet Boundary Condition

L=xr-xl;    %Length of interval
h=L./M;  %Mesh spacing


x_1=xl:h(1):xr;  %x-axis grid
x_2=xl:h(2):xr;
x_3=xl:h(3):xr;

% Initial Guess
phi0_1=1/2*sin(pi*x_1)-1/10*sin(32*pi*x_1);
phi0_2=1/2*sin(pi*x_2)-1/10*sin(32*pi*x_2);
phi0_3=1/2*sin(pi*x_3)-1/10*sin(32*pi*x_3);
phi0_1(1)=lbc;
phi0_1(N(1))=rbc;
phi0_2(1)=lbc;
phi0_2(N(2))=rbc;
phi0_3(1)=lbc;
phi0_3(N(3))=rbc;

% Right Hand Side
f3_1=cos(pi*x_1).*(-(pi^2/4)*x_1.^3-2*pi^2*x_1.^2+3/2*x_1+4)-...
    sin(pi*x_1).*(3/2*pi*x_1.^2+8*pi*x_1);
f3_2=cos(pi*x_2).*(-(pi^2/4)*x_2.^3-2*pi^2*x_2.^2+3/2*x_2+4)-...
    sin(pi*x_2).*(3/2*pi*x_2.^2+8*pi*x_2);
f3_3=cos(pi*x_3).*(-(pi^2/4)*x_3.^3-2*pi^2*x_3.^2+3/2*x_3+4)-...
    sin(pi*x_3).*(3/2*pi*x_3.^2+8*pi*x_3);

% Exact Answer
phi_e_128=2+x_1/4+cos(pi*x_1).*(x_1.^3/4+2*x_1.^2);
phi_e_256=2+x_2/4+cos(pi*x_2).*(x_2.^3/4+2*x_2.^2);
phi_e_512=2+x_3/4+cos(pi*x_3).*(x_3.^3/4+2*x_3.^2);

% Phi
ind=num_it+1;

phi_128=zeros(ind,N(1));
phi_256=zeros(ind,N(2));
phi_512=zeros(ind,N(3));

phi_128(1,:)=phi0_1;
phi_256(1,:)=phi0_2;
phi_512(1,:)=phi0_3;

for ii=2:ind
    phi_128(ii,:)=Gauss_Seidel_1D(phi_128(ii-1,:),f3_1,h(1),1);
    phi_256(ii,:)=Gauss_Seidel_1D(phi_256(ii-1,:),f3_2,h(2),1);
    phi_512(ii,:)=Gauss_Seidel_1D(phi_512(ii-1,:),f3_3,h(3),1);
end

phi_128_end=phi_128(ind,:);
phi_256_end=phi_256(ind,:);
phi_512_end=phi_512(ind,:);

% Error
err_128=phi_e_128-phi_128_end;
err_256=phi_e_256-phi_256_end;
err_512=phi_e_512-phi_512_end;

% Norms
norm_128=zeros(1,3);
norm_256=zeros(1,3);
norm_512=zeros(1,3);

norm_128(1,1)=max(abs(err_128));
norm_256(1,1)=max(abs(err_256));
norm_512(1,1)=max(abs(err_512));

norm_128(1,2)=sum(abs(err_128))/(N(1)+1);
norm_256(1,2)=sum(abs(err_256))/(N(2)+1);
norm_512(1,2)=sum(abs(err_512))/(N(3)+1);

norm_128(1,3)=sqrt(sum(err_128.^2)/(N(1)+1));
norm_256(1,3)=sqrt(sum(err_256.^2)/(N(2)+1));
norm_512(1,3)=sqrt(sum(err_512.^2)/(N(3)+1));

% Order of the error norms
o_128=zeros(1,3);
o_256=zeros(1,3);
o_512=zeros(1,3);

o_128(1,1)=NaN;
o_256(1,1)=NaN;
o_512(1,1)=NaN;

for ii=1:3
    o_128(1,ii)=NaN;
    o_256(1,ii)=log(norm_128(1,ii)/norm_256(1,ii))/log(h(1)/h(2));
    o_512(1,ii)=log(norm_256(1,ii)/norm_512(1,ii))/log(h(2)/h(3));
end

norms=[norm_128;norm_256;norm_512];
order=[o_128;o_256;o_512];

array=[M',norms,order];

table=array2table(array,'VariableNames',{'M','Loo_Norm','L1_Norm','L2_Norm',...
    'Order_Loo','Order_L1','Order_L2'});

disp('Problem 3 Table')
disp(table)

r_128=zeros(ind,N(1));
r_256=zeros(ind,N(2));
r_512=zeros(ind,N(3));

norm_128=zeros(ind,1);
norm_256=zeros(ind,1);
norm_512=zeros(ind,1);

% Residuals
for jj=1:ind
    r_128(jj,:)=residual_1D(phi_128(jj,:),f3_1,h(1));
    r_256(jj,:)=residual_1D(phi_256(jj,:),f3_2,h(2));
    r_512(jj,:)=residual_1D(phi_512(jj,:),f3_3,h(3));
end

iterations=0:num_it;
for kk=1:ind
    norm_128(kk,1)=max(abs(r_128(kk,:)));
    norm_256(kk,1)=max(abs(r_256(kk,:)));
    norm_512(kk,1)=max(abs(r_512(kk,:)));
end

figure()
semilogy(iterations,norm_128,iterations,norm_256,iterations,norm_512)
title('Problem 3');
xlabel('Iterations');
ylabel('Residual Norm');
legend('128 Elements','256 Elements','512 Elements');
end
%% Problem 4
function p4()
num_it=1e6;

M=2.^[7,8,9];   %Number of elements
N=M+1;          %Number of nodes

xl=-1;  %Left bound
xr=1;   %Right bound
lbc=0;  %Left Dirichlet Boundary Condition
rbc=0;  %Right Dirichlet Boundary Condition

L=xr-xl;    %Length of interval
h=L./M;  %Mesh spacing


x_1=xl:h(1):xr;  %x-axis grid
x_2=xl:h(2):xr;
x_3=xl:h(3):xr;

% Initial Guess
phi0_1=1/2*sin(pi*x_1)-1/10*sin(32*pi*x_1);
phi0_2=1/2*sin(pi*x_2)-1/10*sin(32*pi*x_2);
phi0_3=1/2*sin(pi*x_3)-1/10*sin(32*pi*x_3);
phi0_1(1)=lbc;
phi0_1(N(1))=rbc;
phi0_2(1)=lbc;
phi0_2(N(2))=rbc;
phi0_3(1)=lbc;
phi0_3(N(3))=rbc;

% Right Hand Side
f3_1=cos(pi*x_1).*(-(pi^2/4)*x_1.^3-2*pi^2*x_1.^2+3/2*x_1+4)-...
    sin(pi*x_1).*(3/2*pi*x_1.^2+8*pi*x_1);
f3_2=cos(pi*x_2).*(-(pi^2/4)*x_2.^3-2*pi^2*x_2.^2+3/2*x_2+4)-...
    sin(pi*x_2).*(3/2*pi*x_2.^2+8*pi*x_2);
f3_3=cos(pi*x_3).*(-(pi^2/4)*x_3.^3-2*pi^2*x_3.^2+3/2*x_3+4)-...
    sin(pi*x_3).*(3/2*pi*x_3.^2+8*pi*x_3);

% Exact Answer
phi_e_128=2+x_1/4+cos(pi*x_1).*(x_1.^3/4+2*x_1.^2);
phi_e_256=2+x_2/4+cos(pi*x_2).*(x_2.^3/4+2*x_2.^2);
phi_e_512=2+x_3/4+cos(pi*x_3).*(x_3.^3/4+2*x_3.^2);

% Phi
ind=num_it+1;

phi_128=zeros(ind,N(1));
phi_256=zeros(ind,N(2));
phi_512=zeros(ind,N(3));

phi_128(1,:)=phi0_1;
phi_256(1,:)=phi0_2;
phi_512(1,:)=phi0_3;

for ii=2:ind
    phi_128(ii,:)=Gauss_Seidel_1D(phi_128(ii-1,:),f3_1,h(1),1);
    phi_256(ii,:)=Gauss_Seidel_1D(phi_256(ii-1,:),f3_2,h(2),1);
    phi_512(ii,:)=Gauss_Seidel_1D(phi_512(ii-1,:),f3_3,h(3),1);
end

phi_128_end=phi_128(ind,:);
phi_256_end=phi_256(ind,:);
phi_512_end=phi_512(ind,:);

% Error
err_128=phi_e_128-phi_128_end;
err_256=phi_e_256-phi_256_end;
err_512=phi_e_512-phi_512_end;

% Norms
norm_128=zeros(1,3);
norm_256=zeros(1,3);
norm_512=zeros(1,3);

norm_128(1,1)=max(abs(err_128));
norm_256(1,1)=max(abs(err_256));
norm_512(1,1)=max(abs(err_512));

norm_128(1,2)=sum(abs(err_128))/(N(1)+1);
norm_256(1,2)=sum(abs(err_256))/(N(2)+1);
norm_512(1,2)=sum(abs(err_512))/(N(3)+1);

norm_128(1,3)=sqrt(sum(err_128.^2)/(N(1)+1));
norm_256(1,3)=sqrt(sum(err_256.^2)/(N(2)+1));
norm_512(1,3)=sqrt(sum(err_512.^2)/(N(3)+1));

% Order of the error norms
o_128=zeros(1,3);
o_256=zeros(1,3);
o_512=zeros(1,3);

o_128(1,1)=NaN;
o_256(1,1)=NaN;
o_512(1,1)=NaN;

for ii=1:3
    o_128(1,ii)=NaN;
    o_256(1,ii)=log(norm_128(1,ii)/norm_256(1,ii))/log(h(1)/h(2));
    o_512(1,ii)=log(norm_256(1,ii)/norm_512(1,ii))/log(h(2)/h(3));
end

norms=[norm_128;norm_256;norm_512];
order=[o_128;o_256;o_512];

array=[M',norms,order];

table=array2table(array,'VariableNames',{'M','Loo_Norm','L1_Norm','L2_Norm',...
    'Order_Loo','Order_L1','Order_L2'});

disp('Problem 4 Table')
disp(table)

r_128=zeros(ind,N(1));
r_256=zeros(ind,N(2));
r_512=zeros(ind,N(3));

norm_128=zeros(ind,1);
norm_256=zeros(ind,1);
norm_512=zeros(ind,1);

% Residuals
for jj=1:ind
    r_128(jj,:)=residual_1D(phi_128(jj,:),f3_1,h(1));
    r_256(jj,:)=residual_1D(phi_256(jj,:),f3_2,h(2));
    r_512(jj,:)=residual_1D(phi_512(jj,:),f3_3,h(3));
end

iterations=0:num_it;
for kk=1:ind
    norm_128(kk,1)=max(abs(r_128(kk,:)));
    norm_256(kk,1)=max(abs(r_256(kk,:)));
    norm_512(kk,1)=max(abs(r_512(kk,:)));
end

figure()
semilogy(iterations,norm_128,iterations,norm_256,iterations,norm_512)
title('Problem 4');
xlabel('Iterations');
ylabel('Residual Norm');
legend('128 Elements','256 Elements','512 Elements');
end
%% Point Jacobi
function [phi] = point_Jacobi(phi,f,h,n)
M=length(phi);

b=h^2*f;

new_phi=phi;
for kk=1:n
    for ii=2:M-1
        new_phi(ii)=0.5*(phi(ii-1)+phi(ii+1)-b(ii));
    end
end

phi=new_phi;
end
%% Gauss-Seidel 1D
function [phi] = Gauss_Seidel_1D(phi,f,h,n)
M=length(phi);

b=h^2*f;

for kk=1:n
    for ii=2:M-1
        phi(ii)=0.5*(phi(ii-1)+phi(ii+1)-b(ii));
    end
end

end
%% SOR_1D
function [phi] = SOR_1D(phi,f,h,w,n)
M=length(phi);

b=h^2*f;

for kk=1:n
    for ii=2:M-1
        phi(ii)=phi(ii)+w*(0.5*(phi(ii-1)+phi(ii+1)-b(ii))-phi(ii));
    end
end

end
%% Residual 1D
function [r]=residual_1D(phi,f,h)
M=length(phi);

r=zeros(1,M);
for ii=2:M-1
    r(ii)=f(ii)-(phi(ii+1)-2*phi(ii)+phi(ii-1))/h^2;
end

end
%% Gauss-Seidel 2D
function [phi] = Gauss_Seidel_2D(phi,f,h,n)
[M,N]=size(phi);

b=h^2*f;

for kk=1:n
    for jj=2:N-1
        for ii=2:M-1
            phi(ii,jj)=0.25*(phi(ii-1,jj)+phi(ii+1,jj)+phi(ii,jj-1)+phi(ii,jj+1)-b(ii,jj));
        end
    end
end

end
%% SOR 2D
function [phi] = SOR_2D(phi,f,h,w,n)
[M,N]=size(phi);

b=h^2*f;

for kk=1:n
    for jj=2:N-1
        for ii=2:M-1
            phi(ii,jj)=phi(ii,jj)+w*(0.25*(phi(ii-1,jj)+phi(ii+1,jj)+phi(ii,jj-1)+phi(ii,jj+1)-b(ii,jj))-phi(ii,jj));
        end
    end
end
end
%% Residual 2D
function [r]=residual_2D(phi,f,h)
[M,N]=size(phi);

r=zeros(N,M);
for jj=2:N-1
    for ii=2:M-1
        r(ii,jj)=f(ii,jj)-1/h^2*((phi(ii-1,jj)+phi(ii+1,jj)+phi(ii,jj-1)+phi(ii,jj+1)-4*phi(ii,jj)));
    end
end
end