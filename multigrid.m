function phi=multigrid(phi,f,h)
[r,c]=size(phi);
M=r-2;
N=c-2;

phi=Gauss_Seidel_2D(phi,f,h,1);
phi=Neumann0(phi);

if M>2 || N>2
    rh=residual_2D(phi,f,h);
    rh=Neumann0(rh);
    m=M/2;
    n=N/2;
    
    r2h=restr(rh);
    e2h=zeros(m+2,n+2);
    
    e2h=multigrid(e2h,r2h,2*h);
    eh=prolong(e2h);
    
    phi=phi+eh;
    
    phi=Gauss_Seidel_2D(phi,f,h,1);
    
    phi=Neumann0(phi);
end

end