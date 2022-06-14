function [U] = singleLayerlFMM(Gamma,Vh,uh)

M = Vh.uqm(Gamma);
[X,W] = Gamma.qud;
U = M.'*(W.*evallFmm(X,W.*(M*uh)));



end

