function [LinfR,phi]=poisson(phi,f,h,nIterMax,eps)
% [r,c]=size(phi);
% M=r-2;
% N=c-2;

count=1;
curLinfR=1e28;

while count<=nIterMax && curLinfR>eps
    phi=multigrid(phi,f,h);
    r=residual_2D(phi,f,h);
    
    curLinfR=Linf(r);
    LinfR(count)=curLinfR;
    count=count+1;
end

end