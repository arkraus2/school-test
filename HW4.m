clc;
M=8;
N=8;

xl=0;
xr=1;
yb=0;
yt=1;

L=xr-xl;
h=L/M;

x=xl-h/2:h:xr+h/2;
y=yb-h/2:h:yt+h/2;

phi=zeros(M+2,N+2);
f=zeros(M+2,N+2);

for jj=1:N+2
    for ii=1:M+2
        f(ii,jj)=-9*pi*x(ii)*cos(2*pi*y(jj)^2)*(2*sin(3*pi*x(ii)^3)+9*pi*x(ii)^3*cos(3*pi*x(ii)^3))-...
            4*pi*cos(3*pi*x(ii)^3)*(sin(2*pi*y(jj)^2)+4*pi*y(jj)^2*cos(2*pi*y(jj)^2));
    end
end

% f=-9*pi.*x.*cos(2*pi*y.^2).*(2*sin(3*pi*x.^3)+9*pi*x.^3.*cos(3*pi*x.^3))-...
%     4*pi*cos(3*pi*x.^3).*(sin(2*pi*y.^2)+4*pi*y.^2.*cos(2*pi*y.^2));

res=residual_2D(phi,f,h);
res=Neumann0(res);
LinfNorm=Linf(res);

[LinfR,phi]=poisson(phi,f,h,2,0);

% % for jj=1:N+2
% %     for ii=1:M+2
% %         fprintf('phi(%1.0f,%1.0f) = %+1.16e\n',ii-1,jj-1,phi(ii,jj));
% %     end
% % end
% 
% res2=residual_2D(phi,f,h);
% res2=Neumann0(res2);
% LinfNorm2=Linf(res2);