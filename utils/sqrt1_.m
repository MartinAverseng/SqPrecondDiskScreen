function [res] = sqrt1_(X)

res = real(sqrt(1 - X));
res(res>1) = 1;
res(X > 1) = 0;

end

