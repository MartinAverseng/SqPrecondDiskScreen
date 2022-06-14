function [X,W] = TriGaussABC_3d(A,B,C,N)
% Computes the Gaussian quadrature of order N where A B C is a triange in
% R3.


xwref = TriGaussPoints(N);
x = xwref(:,1);
y = xwref(:,2);
w = xwref(:,3);
AB = B - A;
AC = C - A;
Jac = norm(cross(AB,AC));

% find number of Gauss points

X = A + x*AB + y*AC;
W = Jac*w/2;

end