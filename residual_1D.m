function [r]=residual_1D(phi,f,h)
M=length(phi);

r=zeros(1,M);
for ii=2:M-1
    r(ii)=f(ii)-(phi(ii+1)-2*phi(ii)+phi(ii-1))/h^2;
end

end