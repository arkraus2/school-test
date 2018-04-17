clc;

% for jj=1:N+2
%     for ii=1:M+2
%         fprintf('res(%2.0f,%2.0f) = %+1.16e\n',ii-1,jj-1,res(ii,jj))
%     end
% end

% for jj=1:N+2
%     for ii=1:M+2
%         fprintf('phi(%2.0f,%2.0f) = %+1.16e\n',ii-1,jj-1,phi(ii,jj))
%     end
% end

% for jj=1:N+2
%     for ii=1:M+2
%         fprintf('f(%2.0f,%2.0f) = %+1.16e\n',ii-1,jj-1,f(ii,jj))
%     end
% end

% for jj=1:N+2
%     for ii=1:M+2
%         fprintf('rh(%2.0f,%2.0f) = %+1.16e\n',ii-1,jj-1,rh(ii,jj))
%     end
% end

% h=2*h;
% phi=e2h;
% f=r2h;

% for jj=1:n+2
%     for ii=1:m+2
%         fprintf('r2h(%2.0f,%2.0f) = %+1.16e\n',ii-1,jj-1,r2h(ii,jj))
%     end
% end

% e2h=phi;

[s,t]=size(eh);
for jj=1:t
    for ii=1:s
        fprintf('eh(%2.0f,%2.0f) = %+1.16e\n',ii-1,jj-1,eh(ii,jj))
    end
end

% [r,c]=size(e2h);
% for jj=1:c
%     for ii=1:r
%         fprintf('e2h(%2.0f,%2.0f) = %+1.16e\n',ii-1,jj-1,e2h(ii,jj))
%     end
% end