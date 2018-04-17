function r2h=restr(rh)

[r,c]=size(rh);

M=r-2;
N=c-2;

m=M/2;
n=N/2;
xel=m+2;
yel=n+2;

r2h=zeros(xel,yel);

for jj=2:yel-1
    for ii=2:xel-1
        rh1=rh(2*ii-2,2*jj-2);
        rh2=rh(2*ii-1,2*jj-2);
        rh3=rh(2*ii-2,2*jj-1);
        rh4=rh(2*ii-1,2*jj-1);
        r2h(ii,jj)=0.25*(rh1+rh2+rh3+rh4);
    end
end

r2h=Neumann0(r2h);

end