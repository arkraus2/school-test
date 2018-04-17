function [phi] = SOR_1D(phi,f,h,w,n)
M=length(phi);

b=h^2*f;

for kk=1:n
    for ii=2:M-1
        phi(ii)=phi(ii)+w*(0.5*(phi(ii-1)+phi(ii+1)-b(ii))-phi(ii));
    end
end

end