function d=tri_diagonal(a,b,c,d)
P=length(d);

for ii = 2:P
    b(ii)=b(ii) - c(ii-1)*a(ii)/b(ii-1);
    d(ii)=d(ii) - d(ii-1)*a(ii)/b(ii-1);
end

d(P) = d(P)/b(P);
for jj = P-1:-1:1
    d(jj) = (d(jj) - c(jj)*d(jj+1))/b(jj);
end

end