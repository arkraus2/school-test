clc;

M_1=256;
phi_1=4.6603455097273294e-01;

M_2=M_1*2;
phi_2=4.6967762130439206e-01;

M_3=M_2*2;
phi_3=4.7138255039699994e-01;

M_4=M_3*2;
phi_4=4.7230196350605191e-01;

M_5=M_4*2;
phi_5=4.7274427087938542e-01;

M_6=M_5*2;
phi_6=4.7296689307598505e-01;


% f1=phi_3; f2=phi_2; f3=phi_1;
% f1=phi_4; f2=phi_3; f3=phi_2;
% f1=phi_5; f2=phi_4; f3=phi_3;
f1=phi_6; f2=phi_5; f3=phi_4;

p=log(abs(f3-f2)/abs(f2-f1))/log(2);
% p=log((f3-f2)/(f2-f1))/log(2);

f=f1+(f1-f2)/(2^p-1);
GCI_12=1.25*abs((f1-f2)/f1)/(2^p-1)*100;
GCI_23=1.25*abs((f2-f3)/f2)/(2^p-1)*100;
approx1=GCI_12/GCI_23*2^p;

fprintf('M = %4.0f\n',M_6)
fprintf('p = %1.8f\n',p)
fprintf('phi at x=0 and t=1.25: %1.16e\n',f)
fprintf('The error band is %0.5f%%\n',GCI_12)