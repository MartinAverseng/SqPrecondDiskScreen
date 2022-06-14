function [res] = omega(X)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

r2 = X(:,1).^2 + X(:,2).^2;
res = sqrt(1 - r2);
res(r2>1) = 0;
res(res>1) =1;



end

