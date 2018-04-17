function alpha=HW9_alpha(u,E)
len=length(u)-1;
alpha=zeros(1,len);

eps=1e-12;

for ii=1:len
    if abs(u(ii+1)-u(ii))>=eps
        alpha(ii)=(E(ii+1)-E(ii))/(u(ii+1)-u(ii));
    else
        alpha(ii)=(u(ii+1)+u(ii))/2;
    end
end

end