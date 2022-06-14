function [X,W] = intType2(A,B,C)
% Computes quadrature for int_T(f)
% where T = triangle ABC
% f assumed to have a singularity of order 1/2 at the point A.
% Method : Duffy trick (y -> xy) combined with var change (x -> u^2)
% Evaluate the integral by sum(f(X).*W)


AB = B - A;
AC = C - A;
K = [AB(1),AC(1)-AB(1);AB(2),AC(2)-AB(2)];
JacK = abs(det(K));


x = [0.013046735741414
0.067468316655508
0.160295215850488
0.283302302935376
0.425562830509184
0.574437169490816
0.716697697064624
0.839704784149512
0.932531683344492
0.986953264258586];

w = [0.033335672154344
0.074725674575290
0.109543181257991
0.134633359654998
0.147762112357376
0.147762112357376
0.134633359654998
0.109543181257991
0.074725674575290
0.033335672154344];
y = [repmat(x,10,1),repelem(x,10)];
w = w.*w'; w = w(:);
% Tensorized quadrature of the square

u = y(:,1);
v = y(:,2);

x1 = u.^2;
x2 = u.^2.*v;
Kx1 = K(1,1)*x1 + K(1,2)*x2;
Kx2 = K(2,1)*x1 + K(2,2)*x2;

X = [Kx1 + A(1),Kx2 + A(2),0*Kx1];
W = JacK*(x1.^2 + x2.^2).^(1/4)./(1 + v.^2).^(1/4).*2.*u.^2.*w;






end

