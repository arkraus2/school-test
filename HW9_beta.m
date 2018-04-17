function beta=HW9_beta(u,G)
len=length(u)-1;
beta=zeros(1,len);

eps=1e-12;

for ii=1:len
    if abs(u(ii+1)-u(ii))>=eps
        beta(ii)=(G(ii+1)-G(ii))/(u(ii+1)-u(ii));
    else
        beta(ii)=0;
    end
end

end