function out=HW9_psi(y)

eps=0.1;

abs_y=abs(y);
if abs_y>=eps
    out=abs_y;
else
    out=(y^2+eps^2)/(2*eps);
end

end