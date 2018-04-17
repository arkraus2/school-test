HW4();

%--------------------------------
[r,c]=size(phi);
M=r-2;
N=c-2;

phi=Gauss_Seidel_2D(phi,f,h,1);
phi=Neumann0(phi);

rh=residual_2D(phi,f,h);
rh=Neumann0(rh);
m=M/2;
n=N/2;

r2h=restr(rh);
e2h=zeros(m+2,n+2);
%==============================
h=2*h;
phi=e2h;
f=r2h;
%--------------------------------
[r,c]=size(phi);
M=r-2;
N=c-2;

phi=Gauss_Seidel_2D(phi,f,h,1);
phi=Neumann0(phi);

rh=residual_2D(phi,f,h);
rh=Neumann0(rh);
m=M/2;
n=N/2;

r2h=restr(rh);
e2h=zeros(m+2,n+2);
%==============================
h=2*h;
phi=e2h;
f=r2h;
%--------------------------------
[r,c]=size(phi);
M=r-2;
N=c-2;

phi=Gauss_Seidel_2D(phi,f,h,1);
phi=Neumann0(phi);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%===============================
e2h=phi;
h=h/2;
%--------------------------------
eh=prolong(e2h);