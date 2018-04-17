function [phi] = Gauss_Seidel_1D(phi,f,h,n)
M=length(phi);

b=h^2*f;

for kk=1:n
    for ii=2:M-1
        phi(ii)=0.5*(phi(ii-1)+phi(ii+1)-b(ii));
    end
end

end