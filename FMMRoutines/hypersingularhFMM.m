function [U] = hypersingularFMM(Gamma,Vh,k,uh)

Wh = nxgrad(Vh);
M1 = Wh.uqm(Gamma);
[X,W] = Gamma.qud;

U1= M1{1}.'*(W.*evalhFmm(X,k,W.*(M1{1}*uh)));
U2= M1{2}.'*(W.*evalhFmm(X,k,W.*(M1{2}*uh)));

M2 = Vh.uqm(Gamma);
U3 = M2.'*(W.*evalhFmm(X,k,W.*(M2*uh)));

U = U1 + U2 - k^2*U3;


end

