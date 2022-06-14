function [X,W] = intType3(A,B,C,N)

% Computes a quadratue for \int_{U} f(x)dx
% where U is a region of the disk delimited by two segments [A,C] and [B,C]
% and by the arc of circle joining A and B. (A and B are on the boundary of
% the disk, while C is inside the disk).
% f is assumed to have a singularity of the form g/omega where 
% omega(x) = sqrt(1 - ||x||^2)
% Method: combine two changes of variables.
% First one : vertical projection of U on a region V of the uper half sphere
% Second one : gnomonic projection from C to the plane P tangent to the
% sphere at C'. (C' is the projection of C on the upper half sphere)

a = C(1);
b = C(2);
omegaC = sqrt(1 - a^2 - b^2);
Cp = [a b omegaC]; n = Cp;

P_A = A - C;
P_B = B - C;
ndotPA = P_A(1)*n(1) + P_A(2)*n(2) + P_A(3)*n(3);
ndotPB = P_B(1)*n(1) + P_B(2)*n(2) + P_B(3)*n(3);
phiA = C + omegaC^2/(ndotPA)*P_A;
phiB = C + omegaC^2/(ndotPB)*P_B;
phiC = Cp;

[Y,W] = TriGaussABC_3d(phiA,phiB,phiC,N);

Q = Y-C;
CdotQ = C(1)*Q(:,1) + C(2)*Q(:,2);
normQ2 = Q(:,1).^2 + Q(:,2).^2 + Q(:,3).^2;
t = (-CdotQ + sqrt(CdotQ.^2 + omegaC^2*normQ2))./normQ2;
x = C + t.*Q;
P = x - C;
Pdotx = x(:,1).*P(:,1) + x(:,2).*P(:,2) + x(:,3).*P(:,3);
Pdotn = n(1).*P(:,1) + n(2).*P(:,2) + n(3).*P(:,3);
Pix = [x(:,1), x(:,2) 0*x(:,3)];

Jphi = omegaC^4*Pdotx./Pdotn.^3;
W = W./Jphi.*omega(Pix);
X = Pix;



end

