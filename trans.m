function t=trans(t)
[M,N]=size(t);

new_t=t;
for ii=1:M
    for jj=1:N
        new_t(ii,jj)=t(jj,ii);
    end
end
t=new_t;
end