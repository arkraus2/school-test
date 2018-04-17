function psi=psiWENO(a,b,c,d)

IS_0=13*(a-b)^2+3*(a-3*b)^2;
IS_1=13*(b-c)^2+3*(b+c)^2;
IS_2=13*(c-d)^2+3*(3*c-d)^2;

epsilon=1e-6;

alpha0=1/(epsilon+IS_0)^2;
alpha1=6/(epsilon+IS_1)^2;
alpha2=3/(epsilon+IS_2)^2;

w0=alpha0/(alpha0+alpha1+alpha2);
w2=alpha2/(alpha0+alpha1+alpha2);

psi=1/3*w0*(a-2*b+c)+1/6*(w2-1/2)*(b-2*c+d);

end