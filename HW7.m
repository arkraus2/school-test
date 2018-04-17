%% Homework 7

clc;
M=32;
N=16;
CFL=4;

Re=2;
Sc=0.25;

xl=0;
xr=4;

yb=0;
yt=2;

h=(xr-xl)/M;

alpha_uv=1/Re;

outputTime=[0.1,0.5,1,2];

t_step=5e-3;
time=0;


u=zeros(M+1,N+2);
x_u=(xl:h:xr)';
y_u=(yb-h/2:h:yt+h/2)';

u_avg1=2;
u_avg2=-1;
R1=.25;
R2=.25;

for jj=1:length(y_u)
    for ii=1:length(x_u)
        %Inlet 1
        if x_u(ii)==xl && y_u(jj)>=0.5 && y_u(jj)<=.75
            u(ii,jj)=3/2*u_avg1*(1-(R1-(y_u(jj)-.5))^2/R1^2);
        elseif x_u(ii)==xl && y_u(jj)>0.75 && y_u(jj)<=1
            u(ii,jj)=3/2*u_avg1*(1-(-R1+y_u(jj)-.5)^2/R1^2);
        end
        
        %Inlet 2
        if x_u(ii)==xr && y_u(jj)>=1 && y_u(jj)<=1.25
            u(ii,jj)=3/2*u_avg2*(1-(R2-(y_u(jj)-1))^2/R2^2);
        elseif x_u(ii)==xr && y_u(jj)>1.25 && y_u(jj)<=1.5
            u(ii,jj)=3/2*u_avg2*(1-(-R2+y_u(jj)-1)^2/R2^2);
        end
    end
end


v=zeros(M+2,N+1);
x_v=(xl-h/2:h:xr+h/2)';
y_v=(yb:h:yt)';

v_avg3=-1;
R3=.25;

for ii=1:length(x_v)
    for jj=1:length(y_v)
        %Inlet 3
        if y_v(jj)==yt && x_v(ii)>=0.5 && x_v(ii)<=0.75
            v(ii,jj)=3/2*v_avg3*(1-(R3-(x_v(ii)-.5))^2/R3^2);
        elseif y_v(jj)==yt && x_v(ii)>0.75 && x_v(ii)<=1
            v(ii,jj)=3/2*v_avg3*(1-(-R3+x_v(ii)-.5)^2/R3^2);
        end
    end
end



% while time<outputTime(end)
for a=1%:4
    out=0;
    dt=t_step;
    if time<max(outputTime)
        for ii=1:length(outputTime)
            if time<outputTime(ii) && time+dt>=outputTime(ii)
                dt=outputTime(ii)-time;
                out=1;
            end
        end        
    end

    d_uv=alpha_uv*dt/(h^2)/2;
    
    % Solve for u at t^(n+1/2)
    a_u1=-d_uv*ones((M-1)*N,1);
    b_u1=ones((M-1)*N,1)*(1+2*d_uv);
    c_u1=a_u1;
    d_u1=ones((M-1)*N,1);
    for jj=2:N+1
        for ii=2:M
            d_u1((ii-1)+(M-1)*(jj-2))=d_uv*u(ii,jj+1)+(1-2*d_uv)*u(ii,jj)+d_uv*u(ii,jj-1);
            if ii==2
                a_u1((ii-1)+(M-1)*(jj-2))=0;
                if (y_u(jj)>0.5 && y_u(jj)<1)
                    d_u1((ii-1)+(M-1)*(jj-2))=d_u1((ii-1)+(M-1)*(jj-2))+d_uv*u(ii-1,jj);
                end
            end
            if ii==M
                c_u1((ii-1)+(M-1)*(jj-2))=0;
                if (y_u(jj)>1 && y_u(jj)<1.5)
                    d_u1((ii-1)+(M-1)*(jj-2))=d_u1((ii-1)+(M-1)*(jj-2))+d_uv*u(ii+1,jj);
                end
            end
                
        end
    end
    d_u11=d_u1;
    d_u1=tri_diagonal(a_u1,b_u1,c_u1,d_u1);
    
    for jj=2:N+1
        for ii=2:M
            u(ii,jj)=d_u1((ii-1)+(M-1)*(jj-2));
        end
    end
    
    %Horizontal Walls
    for ii=1:length(x_u)
        jj=1;
        %Bottom wall
        u(ii,jj)=-u(ii,jj+1);
        %Top wall
        jj=length(y_u);
        u(ii,jj)=-u(ii,jj-1);
    end
    
    u1=u;
    
    % Solve for v at t^(n+1/2)
    d_v1=ones(M*(N-1),1);
    a_v1=-d_uv*d_v1;
    b_v1=(1+2*d_uv)*d_v1;
    c_v1=a_v1;
    
    for jj=2:N
        for ii=2:M+1
            d_v1((ii-1)+M*(jj-2))=d_uv*v(ii,jj+1)+(1-2*d_uv)*v(ii,jj)+d_uv*v(ii,jj-1);
            if ii==2
                a_v1((ii-1)+M*(jj-2))=0;
                b_v1((ii-1)+M*(jj-2))=b_v1((ii-1)+M*(jj-2))+d_uv;
            end
            if ii==M+1
                c_v1((ii-1)+M*(jj-2))=0;
                b_v1((ii-1)+M*(jj-2))=b_v1((ii-1)+M*(jj-2))+d_uv;
            end
        end
    end
    d_v11=d_v1;
    d_v1=tri_diagonal(a_v1,b_v1,c_v1,d_v1);
    
    for jj=2:N
        for ii=2:M+1
            v(ii,jj)=d_v1((ii-1)+M*(jj-2));
        end
    end
    
    %Dirichlet BC's
    for jj=1:length(y_v)
        %Left wall
        ii=1;
        v(ii,jj)=-v(ii+1,jj);
        %Right wall
        ii=length(x_v);
        v(ii,jj)=-v(ii-1,jj);
    end
    
    %Neumann at outlet
    for ii=1:length(x_v)
        jj=1;
        if (x_v(ii)>1.5 && x_v(ii)<2.5)
            v(ii,jj)=4/3*v(ii,jj+1)-1/3*v(ii,jj+2);
        end
    end
    v1=v;
    
    
    
    % ^*-------------------------------------------------------------
    
    % Solve for u^*
    a_u2=-d_uv*ones((M-1)*N,1);
    b_u2=ones((M-1)*N,1)*(1+2*d_uv);
    c_u2=a_u2;
    d_u2=ones((M-1)*N,1);
    for ii=2:M
        for jj=2:N+1
            d_u2((ii-2)*(N)+(jj-1))=d_uv*u(ii+1,jj)+(1-2*d_uv)*u(ii,jj)+d_uv*u(ii-1,jj);
            if jj==2
                a_u2((ii-2)*(N)+(jj-1))=0;
                b_u2((ii-2)*(N)+(jj-1))=b_u2((ii-2)*(N)+(jj-1))+d_uv;
%                 d_u2((ii-2)*(N)+(jj-1))=d_u2((ii-2)*(N)+(jj-1))+d_uv*u(ii-1,jj);
            end
            if jj==N+1
                c_u2((ii-2)*(N)+(jj-1))=0;
                b_u2((ii-2)*(N)+(jj-1))=b_u2((ii-2)*(N)+(jj-1))+d_uv;
%                 d_u2((ii-2)*(N)+(jj-1))=d_u2((ii-2)*(N)+(jj-1))+d_uv*u(ii+1,jj);
            end
                
        end
    end
    d_u22=d_u2;
    d_u2=tri_diagonal(a_u2,b_u2,c_u2,d_u2);
    
    for ii=2:M
        for jj=2:N+1
            u(ii,jj)=d_u2((ii-2)*(N)+(jj-1));
        end
    end
%     u_out=u;
    %Horizontal Walls
    for ii=1:length(x_u)
        jj=1;
        %Bottom wall
        u(ii,jj)=-u(ii,jj+1);
        %Top wall
        jj=length(y_u);
        u(ii,jj)=-u(ii,jj-1);
    end
    
    
    % Solve for v^*
    d_v2=ones(M*(N-1),1);
    a_v2=-d_uv*d_v2;
    b_v2=(1+2*d_uv)*d_v2;
    c_v2=a_v2;
    
    for ii=2:M+1
        for jj=2:N
            d_v2((ii-2)*(N-1)+(jj-1))=d_uv*v(ii+1,jj)+(1-2*d_uv)*v(ii,jj)+d_uv*v(ii-1,jj);
            if (ii-2)*(N-1)+(jj-1)==(ii-2)*(N-1)+1
                a_v2((ii-2)*(N-1)+(jj-1))=0;
%                 d_v2((ii-2)*(N-1)+(jj-1))=d_v2((ii-2)*(N-1)+(jj-1))+d_uv*v(ii,jj-1);
                if (x_v(ii)>=1.5 && x_v(ii)<=2.5)
                    b_v2((ii-2)*(N-1)+(jj-1))=b_v2((ii-2)*(N-1)+(jj-1))-4/3*d_uv;
                    c_v2((ii-2)*(N-1)+(jj-1))=c_v2((ii-2)*(N-1)+(jj-1))+1/3*d_uv;
                end
            end
            if (ii-2)*(N-1)+(jj-1)==(ii-1)*(N-1)
                c_v2((ii-2)*(N-1)+(jj-1))=0;
                if (x_v(ii)>0.5 && x_v(ii)<1)
                    d_v2((ii-2)*(N-1)+(jj-1))=d_v2((ii-2)*(N-1)+(jj-1))+d_uv*v(ii,jj+1);
                end
            end
        end
    end
    d_v22=d_v2;
    d_v2=tri_diagonal(a_v2,b_v2,c_v2,d_v2);
    
    for ii=2:M+1
        for jj=2:N
            v(ii,jj)=d_v2((ii-2)*(N-1)+(jj-1));
        end
    end
    
    %Dirichlet BC's
    for jj=1:length(y_v)
        %Left wall and inlet 1
        ii=1;
        v(ii,jj)=-v(ii+1,jj);
        %Right wall and inlet 2
        ii=length(x_v);
        v(ii,jj)=-v(ii-1,jj);
    end
    
    %Neumann at outlet
    for ii=1:length(x_v)
        jj=1;
        if (x_v(ii)>1.5 && x_v(ii)<2.5)
            v(ii,jj)=4/3*v(ii,jj+1)-1/3*v(ii,jj+2);
        end
    end
    
    %Calculate Outlet correction velocity
    ucorr=0;
    for jj=N+1
        for ii=1:M+2
            ucorr=ucorr-v(ii,jj)*h;
        end
    end

    for jj=1:N+2
        for ii=[1,M+1]
            if ii==1
                ucorr=ucorr+u(ii,jj)*h;
            else
                ucorr=ucorr-u(ii,jj)*h;
            end
        end
    end
    
    for jj=1
        for ii=1:M+2
            if (x_v(ii)>1.5 && x_v(ii)<2.5)
                v(ii,jj)=-ucorr;
            end
        end
    end
    
    phi=zeros(M+2,N+2);
    rhs=zeros(M+2,N+2);
    for jj=2:N+1
        for ii=2:M
            rhs(ii,jj)=rhs(ii,jj)+1/(dt*h)*(u(ii,jj)-u(ii-1,jj));
        end
    end
    for jj=2:N
        for ii=2:M+1
            rhs(ii,jj)=rhs(ii,jj)+1/(dt*h)*(v(ii,jj)-v(ii,jj-1));
        end
    end
    [~,phi]=poisson(phi,rhs,1,50,1e-10);
    
    time=time+dt;
end