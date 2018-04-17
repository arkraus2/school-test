function [phi] = point_Jacobi(phi,f,h,n)
M=length(phi);

b=h^2*f;

new_phi=phi;
for kk=1:n
    for ii=2:M-1
        new_phi(ii)=0.5*(phi(ii-1)+phi(ii+1)-b(ii));
    end
end

phi=new_phi;
end