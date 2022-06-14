function [U] = hypersingularlFMM(Gamma,Vh,uh)

Wh = nxgrad(Vh);
M1 = Wh.uqm(Gamma);
[X,W] = Gamma.qud;

U1= M1{1}.'*(W.*evallFmm(X,W.*(M1{1}*uh)));
U2= M1{2}.'*(W.*evallFmm(X,W.*(M1{2}*uh)));

U = U1 + U2;


end

