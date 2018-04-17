function [phi] = Gauss_Seidel_2D(phi,f,h,n)
[r,c]=size(phi);

b=h^2*f;

for kk=1:n
    for jj=2:c-1
        for ii=2:r-1
            phi(ii,jj)=0.25*(phi(ii-1,jj)+phi(ii+1,jj)+phi(ii,jj-1)+phi(ii,jj+1)-b(ii,jj));
        end
    end
end

end