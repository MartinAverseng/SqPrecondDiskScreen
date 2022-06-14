function [x] = barycentricCoord(A,B,C,X)

AB = B - A;
AB2 = AB(1)^2 + AB(2)^2 + AB(3)^2;
AC = C - A;
AC2 = AC(1)^2 + AC(2)^2 + AC(3)^2;
ABAC = AB(1)*AC(1) + AB(2)*AC(2) + AB(3)*AC(3);
detK = AB2*AC2 - ABAC^2;
Kinv = 1/detK*[AC2 -ABAC;-ABAC AB2];
AX = X - A;
L1 = AB(1)*AX(:,1) + AB(2)*AX(:,2) + AB(3)*AX(:,3);
L2 = AC(1)*AX(:,1) + AC(2)*AX(:,2) + AC(3)*AX(:,3);

L = [L1 L2];


x = L*Kinv';



end

