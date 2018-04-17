function u=HW9_bc(u)

u(1)=u(length(u)-3);
u(2)=u(length(u)-2);
u(length(u)-1)=u(3);
u(length(u))=u(4);

end