function [r]=residual_2D(phi,f,h)
[row,c]=size(phi);

r=zeros(row,c);
for jj=2:c-1
    for ii=2:row-1
        r(ii,jj)=f(ii,jj)-1/h^2*((phi(ii-1,jj)+phi(ii+1,jj)+phi(ii,jj-1)+phi(ii,jj+1)-4*phi(ii,jj)));
    end
end
end