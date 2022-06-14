function [y] = evallFmm(X,u)


source = X.';
iprec = 1;
nsource = size(X,1);
ifcharge = 1;
charge = u.';
ifdipole = 0;
dipstr = 0;
dipvec = zeros(3,nsource);
ifpot = 1;
iffld = 0;
[U]=lfmm3dpart(iprec,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld);
y = U.pot.';


end

