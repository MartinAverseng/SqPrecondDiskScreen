function [ v ] = TrefethenSqrt2(A,N,u,M,min_m,max_M)


k2 = min_m/max_M; % elliptic functions parameter k^2
Kp = ellipke(1-k2);
t = 1i*(.5:N)*Kp/N;
[sn, cn, dn] = ellipj(imag(t),1-k2);
cn = 1./cn; dn = dn.*cn; sn = 1i*sn.*cn;
w = sqrt(min_m)*sn;
dzdt = cn.*dn;
lambda = M*u;

n = length(lambda);
vj = zeros(n,N);

parfor (j = 1:N,4)
    C = A-w(j)^2*M;
    vj(:,j) = - dzdt(j)*(C\lambda);
end
v = 0;
for j = 1:N
    v = v + vj(:,j);
end
v = (-2*Kp*sqrt(min_m)/(pi*N))*(M*v);

clearvars -except v;

end

