function [phi] = SOR_2D(phi,f,h,w,n)
[M,N]=size(phi);

b=h^2*f;

for kk=1:n
    for jj=2:N-1
        for ii=2:M-1
            phi(ii,jj)=phi(ii,jj)+w*(0.25*(phi(ii-1,jj)+phi(ii+1,jj)+phi(ii,jj-1)+phi(ii,jj+1)-b(ii,jj))-phi(ii,jj));
        end
    end
end
end