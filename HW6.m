%% Homework 6

clc;
M=16;
N=8;
CFL=4;

Re=2;
Sc=0.25;

xl=0;
xr=4;

yb=0;
yt=2;

h=(xr-xl)/M;

alpha_uv=1/Re;
alpha_Y=1/(Re*Sc);

outputTime=[0.1,0.5,1,10];

t_uv_stable=CFL/4*h^2/alpha_uv;
t_Y_stable=CFL/4*h^2/alpha_Y;

t_stable=min([t_uv_stable,t_Y_stable]);
time=0;


u=zeros(M+1,N+2);
x_u=(xl:h:xr)';
y_u=(yb-h/2:h:yt+h/2)';

for jj=1:length(y_u)
    for ii=1:length(x_u)
        %Inlet 1
        if x_u(ii)==xl && y_u(jj)>=0.5 && y_u(jj)<=1
            u(ii,jj)=2;
        end
        
        %Inlet 2
        if x_u(ii)==xr && y_u(jj)>=1 && y_u(jj)<=1.5
            u(ii,jj)=-1;
        end
    end
end

v=zeros(M+2,N+1);
x_v=(xl-h/2:h:xr+h/2)';
y_v=(yb:h:yt)';

for ii=1:length(x_v)
    for jj=1:length(y_v)
        %Inlet 3
        if y_v(jj)==yt && x_v(ii)>=0.5 && x_v(ii)<=1
            v(ii,jj)=-1;
        end
    end
end

Y=zeros(M+2,N+2);
x_Y=(xl-h/2:h:xr+h/2)';
y_Y=(yb-h/2:h:yt+h/2)';

for jj=1:length(y_Y)
    for ii=1:length(x_Y)
        %Inlet 1
        if x_Y(ii)==xl-h/2 && y_Y(jj)>=0.5 && y_Y(jj)<=1
            Y(ii,jj)=2;
        end
        
        %Inlet 2
        if x_Y(ii)==xr+h/2 && y_Y(jj)>=1 && y_Y(jj)<=1.5
            Y(ii,jj)=.5;
        end
    end
end

% while time<outputTime(end)
for a=1:4
    out=0;
    dt=t_stable;
    if time<max(outputTime)
        for ii=1:length(outputTime)
            if time<outputTime(ii) && time+dt>=outputTime(ii)
                dt=outputTime(ii)-time;
                out=1;
            end
        end        
    end

    d_uv=alpha_uv*dt/(h^2)/2;
    d_Y=alpha_Y*dt/(h^2)/2;
    
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
    
    
    % Solve for Y at t^(n+1/2)
    d_Y1=ones(M*N,1);
    a_Y1=-d_Y*d_Y1;
    b_Y1=(1+2*d_Y)*d_Y1;
    c_Y1=a_Y1;
    
    for jj=2:N+1
        for ii=2:M+1
            d_Y1((ii-1)+M*(jj-2))=d_Y*Y(ii,jj+1)+(1-2*d_Y)*Y(ii,jj)+d_Y*Y(ii,jj-1);
            if ii==2
                a_Y1((ii-1)+(M)*(jj-2))=0;
                if (y_Y(jj)>0.5 && y_Y(jj)<1)
                    b_Y1((ii-1)+(M)*(jj-2))=b_Y1((ii-1)+(M)*(jj-2))+d_Y;
                    d_Y1((ii-1)+(M)*(jj-2))=d_Y1((ii-1)+(M)*(jj-2))+2*d_Y*1;
%                     d_Y1((ii-1)+(M)*(jj-2))=d_Y1((ii-1)+(M)*(jj-2))+d_Y*Y(ii-1,jj);
                else
                    b_Y1((ii-1)+(M)*(jj-2))=b_Y1((ii-1)+(M)*(jj-2))-d_Y;
                end
            end
            if ii==M+1
                c_Y1((ii-1)+(M)*(jj-2))=0;
                if (y_Y(jj)>1 && y_Y(jj)<1.5)
                    b_Y1((ii-1)+(M)*(jj-2))=b_Y1((ii-1)+(M)*(jj-2))+d_Y;
                    d_Y1((ii-1)+(M)*(jj-2))=d_Y1((ii-1)+(M)*(jj-2))+2*d_Y*.25;
%                     d_Y1((ii-1)+(M)*(jj-2))=d_Y1((ii-1)+(M)*(jj-2))+d_Y*Y(ii+1,jj);
                else
                    b_Y1((ii-1)+(M)*(jj-2))=b_Y1((ii-1)+(M-1)*(jj-2))-d_Y;
                end
            end
        end
    end
    d_Y11=d_Y1;
    d_Y1=tri_diagonal(a_Y1,b_Y1,c_Y1,d_Y1);
    
    for jj=2:N+1
        for ii=2:M+1
            Y(ii,jj)=d_Y1((ii-1)+(M)*(jj-2));
        end
    end
    
    %Horizontal Walls
    for ii=1:length(x_Y)
        %Bottom wall
        jj=1;
        Y(ii,jj)=Y(ii,jj+1);
        %Top wall
        jj=length(y_Y);
        if (x_Y(ii)>0.5 && x_Y(ii)<1)
            Y(ii,jj)=-Y(ii,jj-1);
        else
            Y(ii,jj)=Y(ii,jj-1);
        end
    end
    
    %Vertical Walls
    for jj=1:length(y_Y)
        %Left Wall
        ii=1;
        if (y_Y(jj)>0.5 && y_Y(jj)<1)
            Y(ii,jj)=2*1-Y(ii+1,jj);
        else
            Y(ii,jj)=Y(ii+1,jj);
        end
        %Right Wall
        ii=length(x_Y);
        if (y_Y(jj)>1 && y_Y(jj)<1.5)
            Y(ii,jj)=2*0.25-Y(ii-1,jj);
        else
            Y(ii,jj)=Y(ii-1,jj);
        end
    end
    Y1=Y;
    
    % t^(n+1)-------------------------------------------------------------
    
    % Solve for u at t^(n+1)
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
    
    
    % Solve for v at t^(n+1)
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
    
    
    % Solve for Y at t^(n+1)
    d_Y2=ones(M*N,1);
    a_Y2=-d_Y*d_Y2;
    b_Y2=(1+2*d_Y)*d_Y2;
    c_Y2=a_Y2;
    
    for ii=2:M+1
        for jj=2:N+1
            d_Y2((ii-2)*N+(jj-1))=d_Y*Y(ii+1,jj)+(1-2*d_Y)*Y(ii,jj)+d_Y*Y(ii-1,jj);
            if jj==2
                a_Y2((ii-2)*N+(jj-1))=0;
%                 if (x_Y(ii)>=1.5 && x_Y(ii)<=2.5)
%                     b_Y2((ii-2)*N+(jj-1))=b_Y2((ii-2)*N+(jj-1))+d_Y;
%                 else
%                     b_Y2((ii-2)*N+(jj-1))=b_Y2((ii-2)*N+(jj-1))-d_Y;
%                 end
                b_Y2((ii-2)*N+(jj-1))=b_Y2((ii-2)*N+(jj-1))-d_Y;
%                 d_Y2((ii-2)*N+(jj-1))=d_Y2((ii-2)*N+(jj-1))+d_Y*Y(ii,jj-1);
            end
            if jj==N+1
                c_Y2((ii-2)*N+(jj-1))=0;
                if (x_Y(ii)>0.5 && x_Y(ii)<1)
                    b_Y2((ii-2)*N+(jj-1))=b_Y2((ii-2)*N+(jj-1))+d_Y;
                else
                    b_Y2((ii-2)*N+(jj-1))=b_Y2((ii-2)*N+(jj-1))-d_Y;
                end
%                 d_Y2((ii-2)*N+(jj-1))=d_Y2((ii-2)*N+(jj-1))+d_Y*Y(ii,jj+1);
            end
        end
    end
    d_Y22=d_Y2;
    d_Y2=tri_diagonal(a_Y2,b_Y2,c_Y2,d_Y2);
    
    for ii=2:M+1
        for jj=2:N+1
            Y(ii,jj)=d_Y2((ii-2)*N+(jj-1));
        end
    end
%     Y_out=Y;
    %Horizontal Walls
    for ii=1:length(x_Y)
        %Bottom wall
        jj=1;
        Y(ii,jj)=Y(ii,jj+1);
        %Top wall
        jj=length(y_Y);
        if (x_Y(ii)>0.5 && x_Y(ii)<1)
            Y(ii,jj)=-Y(ii,jj-1);
        else
            Y(ii,jj)=Y(ii,jj-1);
        end
    end
    
    %Vertical Walls
    for jj=1:length(y_Y)
        %Left Wall
        ii=1;
        if (y_Y(jj)>0.5 && y_Y(jj)<1)
            Y(ii,jj)=2*1-Y(ii+1,jj);
        else
            Y(ii,jj)=Y(ii+1,jj);
        end
        %Right Wall
        ii=length(x_Y);
        if (y_Y(jj)>1 && y_Y(jj)<1.5)
            Y(ii,jj)=2*0.25-Y(ii-1,jj);
        else
            Y(ii,jj)=Y(ii-1,jj);
        end
    end
    
    time=time+dt;
end