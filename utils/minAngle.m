function [thetaMin] = minAngle(m)


A = m.vtx(m.elt(:,1),:);
B = m.vtx(m.elt(:,2),:);
C = m.vtx(m.elt(:,3),:);

norm2 = @(X)(sqrt(X(:,1).^2 + X(:,2).^2));
scal = @(X,Y)(abs(X(:,1).*Y(:,1) + X(:,2).*Y(:,2)));
c = norm2(B-A);
b = norm2(C-A);
a = norm2(C-B);
ab = scal(C-B,C-A);
ac = scal(C-B,B-A);
bc = scal(C-A,B-A);


theta1 = acos(ab./(a.*b));
theta2 = acos(ac./(a.*c));
theta3 = acos(bc./(b.*c));
thetaMin = min([theta1;theta2;theta3])*160/(2*pi);

end

