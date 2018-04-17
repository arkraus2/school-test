function rh=prolong(r2h)
[r,c]=size(r2h);

M=r-2;
N=c-2;

m=M*2;
n=N*2;

rh=zeros(m+2,n+2);

for jj=2:N+1
    for ii=2:M+1
        rh(2*ii-2,2*jj-2)=r2h(ii,jj);
        rh(2*ii-1,2*jj-2)=r2h(ii,jj);
        rh(2*ii-2,2*jj-1)=r2h(ii,jj);
        rh(2*ii-1,2*jj-1)=r2h(ii,jj);
    end
end

rh=Neumann0(rh);

end