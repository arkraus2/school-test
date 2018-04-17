function G=HW9_G(sigma,u)
len=length(u);
G=zeros(1,len);

for ii=3:len-2
    up=u(ii+1)-u(ii);
    um=u(ii)-u(ii-1);
    
    S=sign(up);
    
    G(ii)=S*max([0,min([sigma(ii)*abs(up),S*sigma(ii-1)*um])]);
    
end

G=HW9_bc(G);

end