clear; clc;
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
xlabel('Iterations');
ylabel('Residual Norm');
legend('Gauss-Seidel','SOR');