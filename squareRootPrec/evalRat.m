function [ phi ] = evalRat(Np,C0,Aj,Bj,I,D,l)

lambda = I*l;
N = size(lambda,1);
phi = C0*lambda;
phi_j = zeros(N,Np);
parfor (j = 1:Np,4)
    Gj =  I + Bj(j)*D;
    phi_j(:,j) = Gj\lambda; 
end

for j = 1:Np
    phi = phi + Aj(j)*D*phi_j(:,j);
end

clearvars -except phi

end

