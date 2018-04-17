clear; clc;

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

num_it=1000;
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
xlabel('Iterations');
ylabel('Residual Norm');
legend('128 Elements','256 Elements','512 Elements');