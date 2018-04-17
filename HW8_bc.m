function phi=HW8_bc(phi,x,xl,t)

phi(1)=phi(6)-((x(6)-x(1))*(phi(6)-Left(t)))/(x(6)-xl);
phi(2)=phi(5)-((x(5)-x(2))*(phi(5)-Left(t)))/(x(5)-xl);
phi(3)=phi(4)-((x(4)-x(3))*(phi(4)-Left(t)))/(x(4)-xl);

% phi(1)=2*Left(t)-phi(6);
% phi(2)=2*Left(t)-phi(5);
% phi(3)=2*Left(t)-phi(4);

phi(end)=phi(end-5);
phi(end-1)=phi(end-4);
phi(end-2)=phi(end-3);
end
%% Left Boundary at t
function phi=Left(t)

if t<=1/4
    phi=1;
elseif t<=1/2
    phi=0;
elseif t<=1
    phi=(1-cos(4*pi*t))/2;
else
    phi=0;
end

end