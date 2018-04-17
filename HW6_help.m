% clear; 
clc;
% 
% syms h;
% 
% A=[1,1,1,0;
%     0,-h,-2*h,0;
%     0,1/2*h^2,2*h^2,-1];
% 
% B=rref(A);
% disp(B)

% for kk=1:N
% %     disp((kk-1)*(M-1)+1)
%     disp(kk*(M-1))
% end

% for jj=2:N+1
%      disp((jj-2)*(M-1)+1)
% %      disp((jj-1)*(M-1))
% end

% for jj=2:N+1
%     for ii=2:M
%         disp((ii-1)+(M-1)*(jj-2))
%     end
% end

% for mm=1:length(d_u)
%     fprintf('d_u(%3.0f) = %1.16e\n',mm,d_u(mm))
% end



% for ii=2:M+1
%     for jj=2:N
%         disp((ii-2)*(N-1)+(jj-1))
%     end
% end

% for jj=2:N
%     for ii=2:M+1
%         disp((ii-1)+M*(jj-2))
%     end
% end

% for jj=2:N
%     for ii=2:M+1
%         disp((jj-2)*M+1)
%     end
% end

% for jj=2:N
%     disp((jj-1)*M)
% end

% for jj=2:N
%     for ii=2:M+1
%         fprintf('ind = %3.0f\tright = %3.0f\n',(ii-1)+M*(jj-2),(jj-1)*M)
% %         disp(((ii-1)+(M-1)*(jj-2))==((jj-1)*M))
%     end
% end


% for ii=1:length(x_u)
%     jj=1;
%     if (x_u(ii)<1.5 || x_u(ii)>2.5)
%         u(ii,jj)=-u(ii,jj+1);
%     end
% end
% 
% for ii=1:length(x_u)
%     jj=length(y_u);
%     if (x_u(ii)<0.5 || x_u(ii)>1)
%         u(ii,jj)=-u(ii,jj-1);
%     end
% end

% for jj=1:N+2
%     for ii=1:M+1
%         fprintf('u(%2.0f,%2.0f) = %1.16e\n',ii,jj,u(ii,jj))
%     end
% end


%     for jj=1:length(y_v)
%         ii=1;
%         if (y_v(jj)<0.5 || y_v(jj)>1)
%             v(ii,jj)=-v(ii+1,jj);
%         end
%     end
% 
%     for jj=1:length(y_v)
%         ii=length(x_v);
%         if (y_v(jj)<1 || y_v(jj)>1.5)
%             v(ii,jj)=-v(ii-1,jj);
%         end
%     end

% for jj=1:N+1
%     for ii=1:M+2
%         fprintf('v(%2.0f,%2.0f) = %1.16e\n',ii,jj,v(ii,jj))
%     end
% end

% for jj=2:length(y_Y)-1
%     for ii=2:length(x_Y)-1
%         disp((jj-2)*M+1)
%     end
% end

% for jj=2:length(y_Y)-1
%     for ii=2:length(x_Y)-1
% 
%         if ii==2
%             fprintf('ii = %3.0f, jj = %3.0f\tind = %3.0f\n',ii,jj,(ii-1)+(M)*(jj-2))
%         end
% 
%     end
% end

% for mm=1:length(d_Y1)
%     fprintf('d_Y1(%3.0f) = %1.16e\n',mm,d_Y1(mm))
% end

% for mm=1:length(b_Y1)
%     fprintf('b_Y1(%3.0f) = %1.16e\n',mm,b_Y1(mm))
% end


% for jj=1:N+2
%     for ii=1:M+2
%         fprintf('Y(%2.0f,%2.0f) = %1.16e\n',ii,jj,Y(ii,jj))
%     end
% end

% for jj=2:N+1
%     for ii=2:M
%         fprintf('ii=%2.0f, jj=%2.0f, ind=%3.0f\n',ii,jj,(ii-2)*(N)+(jj-1))
%     end
% end

% for ii=2:M
% %      disp((ii-2)*(N)+1)
%      disp((ii-1)*(N))
% end

% for jj=2:N
%     for ii=2:M+1
%         fprintf('ii=%2.0f, jj=%2.0f, ind=%3.0f\n',ii,jj,(ii-2)*(N-1)+(jj-1))
%     end
% end

% for ii=2:M+1
% %      disp((ii-2)*(N-1)+1)
%      disp((ii-1)*(N-1))
% end

% for jj=1:N+1
%     for ii=1:M+2
%         fprintf('v(%2.0f,%2.0f) = %1.16e\n',ii,jj,v(ii,jj))
%     end
% end

% for jj=2:N+1
%     for ii=2:M+1
%         fprintf('ii=%2.0f, jj=%2.0f, ind=%3.0f\n',ii,jj,(ii-2)*N+(jj-1))
%     end
% end

% for jj=1:N+2
%     for ii=1:M+2
%         fprintf('Y(%2.0f,%2.0f) = %+1.16e\n',ii,jj,Y(ii,jj))
%     end
% end



% for mm=1:length(a_Y2)
%     fprintf('a_Y2(%3.0f) = %+1.16e\n',mm,a_Y2(mm))
% end

% for mm=1:length(b_Y2)
%     fprintf('b_Y2(%3.0f) = %+1.16e\n',mm,b_Y2(mm))
% end

% for mm=1:length(d_Y2)
%     fprintf('d_Y2(%3.0f) = %1.16e\n',mm,d_Y2(mm))
% end


% for jj=1:N+2
%     for ii=1:M+1
%         fprintf('u(%2.0f,%2.0f) = %+1.16e\n',ii,jj,u(ii,jj))
%     end
% end
% 
% pause(); clc;
% 
% for jj=1:N+1
%     for ii=1:M+2
%         fprintf('v(%2.0f,%2.0f) = %+1.16e\n',ii,jj,v(ii,jj))
%     end
% end
% 
% pause(); clc;
% 
for jj=1:N+2
    for ii=1:M+2
        fprintf('Y(%2.0f,%2.0f) = %+1.16e\n',ii,jj,Y(ii,jj))
    end
end

% for mm=1:length(d_u2)
%     fprintf('d_u2(%3.0f) = %+1.16e\n',mm,d_u2(mm))
% end

% for mm=1:length(d_v1)
%     fprintf('d_v1(%3.0f) = %+1.16e\n',mm,d_v11(mm))
% end

% for mm=1:length(d_Y1)
%     fprintf('d_Y1(%3.0f) = %+1.16e\n',mm,d_Y1(mm))
% end

% for jj=1:N+2
%     for ii=1:M+2
%         fprintf('Y(%2.0f,%2.0f) = %+1.16e\n',ii,jj,Y1(ii,jj))
%     end
% end

% for mm=1:length(d_Y2)
%     fprintf('d_Y2(%3.0f) = %+1.16e\n',mm,d_Y22(mm))
% end

% for mm=1:length(b_v1)
%     fprintf('b_v1(%3.0f) = %+1.16e\n',mm,b_v1(mm))
% end

% for mm=1:length(d_v22)
%     fprintf('d_v22(%3.0f) = %+1.16e\n',mm,d_v22(mm))
% end