clc;

load('HW9_vars.mat');

% M=[256,512,1024,2048,4096,8192];
M=[256,512,1024,2048];


xl=3;
xr=7;

outputTime=1;

L=xr-xl;

for ii=1:length(M)
    dx=L/M(ii);
    x=xl-dx*3/2:dx:xr+dx*3/2;

    xpos_ind=find(x>6,1);
    xneg_ind=xpos_ind-1;
    
    switch ii
        case 1
            phi=u256;
        case 2
            phi=u512;
        case 3
            phi=u1024;
        case 4
            phi=u2048;
        case 5
            phi=u4096;
        case 6
            phi=u8192;
    end
    
    u_x6=(phi(xpos_ind)+phi(xneg_ind))/2;

    fprintf('M = %4.0f\n',M(ii))
    fprintf('u @x=6: %1.16e\n\n',u_x6)
end