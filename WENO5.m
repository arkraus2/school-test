function out=WENO5(phi,ind,pos,h)

out=1/(12*h)*(-dP(phi,ind-2)+7*dP(phi,ind-1)+7*dP(phi,ind)-dP(phi,ind+1));

if pos>0
    a=dMdP(phi,ind-2)/h;
    b=dMdP(phi,ind-1)/h;
    c=dMdP(phi,ind)/h;
    d=dMdP(phi,ind+1)/h;
    out=out-psiWENO(a,b,c,d);
else
    a=dMdP(phi,ind+2)/h;
    b=dMdP(phi,ind+1)/h;
    c=dMdP(phi,ind)/h;
    d=dMdP(phi,ind-1)/h;
    out=out+psiWENO(a,b,c,d);
end

end
%% Delta+
function dp=dP(phi,ind)
dp=phi(ind+1)-phi(ind);
end
%% Delta- Delta+
function dmdp=dMdP(phi,ind)
dmdp=phi(ind+1)-2*phi(ind)+phi(ind-1);
end