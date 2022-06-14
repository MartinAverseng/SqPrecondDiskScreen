function[ C0,Aj,Bj ] = rationalCoeffSqrt(Np,theta,k)

% Computes rational coefficients C0, Aj, Bj, for an approximation
% sqrt(X - k^2) approx C0 + sum_j Aj X / (1 + Bj X)

c0 = 1;
j = (1:Np)';
aj = 2/(2*Np + 1)*sin(j*pi/(2*Np+1)).^2;
bj = cos(j*pi/(2*Np+1)).^2;

Rp = @(x)(c0 + sum(aj*x./(1 + bj*x)));

C0 = exp(1i*theta/2)*Rp(exp(-1i*theta)-1);
Aj = exp(-1i*theta/2)*aj./(1 + bj*(exp(-1i*theta)-1)).^2;
Bj = exp(-1i*theta)*bj./(1 + bj*(exp(-1i*theta)-1));

C0 = 1i*k*C0;
Aj = -1i*k*Aj/k^2;
Bj = -1i*k*Bj/k^2;

end

