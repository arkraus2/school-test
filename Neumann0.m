function ary=Neumann0(ary)
[m,n]=size(ary);

for ii=1:m
    if ii==1 || ii==m
        ary(ii,1)=ary(ii,2);
        ary(ii,n)=0;
    else
        ary(ii,1)=ary(ii,2);
        ary(ii,n)=ary(ii,n-1);
    end
end

for jj=1:n
    ary(1,jj)=ary(2,jj);
    ary(m,jj)=ary(m-1,jj);
end

ary(1,1)=ary(2,2);
ary(1,n)=ary(1,n-1);
ary(m,1)=ary(m-1,1);
ary(m,n)=ary(m-1,n-1);

end