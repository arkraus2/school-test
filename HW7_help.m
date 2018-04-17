clc;

% for jj=1:N+2
%     for ii=1:M+1
%         fprintf('u(%2.0f,%2.0f) = %+1.16e\n',ii,jj,u(ii,jj))
%     end
% end

% for jj=1:N+1
%     for ii=1:M+2
%         fprintf('v(%2.0f,%2.0f) = %+1.16e\n',ii,jj,v(ii,jj))
%     end
% end

% ucorr=0;
% for jj=N+1
%     for ii=1:M+2
%         ucorr=ucorr-v(ii,jj)*h;
%     end
% end
% 
% for jj=1:N+2
%     for ii=[1,M+1]
%         if ii==1
%             ucorr=ucorr+u(ii,jj)*h;
%         else
%             ucorr=ucorr-u(ii,jj)*h;
%         end
%     end
% end

% phi=zeros(M+2,N+2);
% for jj=2:N+1
%     for ii=2:M
%         phi(ii,jj)=phi(ii,jj)+1/(dt*h)*(u(ii,jj)-u(ii-1,jj));
%     end
% end
% for jj=2:N
%     for ii=2:M+1
%         phi(ii,jj)=phi(ii,jj)+1/(dt*h)*(v(ii,jj)-v(ii,jj-1));
%     end
% end

for jj=1:N+2
    for ii=1:M+2
        fprintf('phi(%2.0f,%2.0f) = %+1.16e\n',ii,jj,phi(ii,jj))
    end
end