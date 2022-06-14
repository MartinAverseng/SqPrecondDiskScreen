function [U] = singleLayerhFMM(Gamma,Vh,k,uh)

M = Vh.uqm(Gamma);
[X,W] = Gamma.qud;
U = M.'*(W.*evalhFmm(X,k,W.*(M*uh)));



end

