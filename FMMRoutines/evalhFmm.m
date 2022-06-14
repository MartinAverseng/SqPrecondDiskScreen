function [y] = evalhFmm(X,k,u)


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
[U]=hfmm3dpart(iprec,k,nsource,source,ifcharge,charge,ifdipole,dipstr,dipvec,ifpot,iffld);
y = (U.pot + 1i*k*charge).';


end

