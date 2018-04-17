function sigma=HW9_sigma(a,dt,dx)
len=length(a);
sigma=zeros(1,len);

for ii=1:len
    sigma(ii)=1/2*(HW9_psi(a(ii))-dt/dx*a(ii)^2);
end

end