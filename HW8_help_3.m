clc;

M_1=256;
phi_1=4.6946860600322543e-01;

M_2=M_1*2;
phi_2=4.7161198683189803e-01;

M_3=M_2*2;
phi_3=4.7239423763552735e-01;

M_4=M_3*2;
phi_4=4.7279160026414668e-01;

M_5=M_4*2;
phi_5=4.7299038692654882e-01;

M_6=M_5*2;
phi_6=4.7308917843340093e-01;


% f1=phi_3; f2=phi_2; f3=phi_1; outM=M_3;
% f1=phi_4; f2=phi_3; f3=phi_2; outM=M_4;
f1=phi_5; f2=phi_4; f3=phi_3; outM=M_5;
% f1=phi_6; f2=phi_5; f3=phi_4; outM=M_6;

p=log(abs(f3-f2)/abs(f2-f1))/log(2);
% p=log((f3-f2)/(f2-f1))/log(2);

f=f1+(f1-f2)/(2^p-1);
GCI_12=1.25*abs((f1-f2)/f1)/(2^p-1)*100;
GCI_23=1.25*abs((f2-f3)/f2)/(2^p-1)*100;
approx1=GCI_12/GCI_23*2^p;

fprintf('M = %4.0f\n',outM)
fprintf('p = %1.8f\n',p)
fprintf('phi at x=0 and t=1.25: %1.16e\n',f)
fprintf('The error band is %0.5f%%\n',GCI_12)