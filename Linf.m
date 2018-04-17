function L=Linf(res)
[r,c]=size(res);

for jj=1:c
    for ii=1:r
        if ii==1 && jj==1
            L=abs(res(ii,jj));
        else
            temp=abs(res(ii,jj));
            if temp>L
                L=temp;
            end
        end
    end
end

end