function [X] = triCoord(A,B,C,x)

X = 0*x;
X(:,1) = A(1) + x(:,1)*(B(1) - A(1)) + x(:,2)*(C(1) - A(1));
X(:,2) = A(2) + x(:,1)*(B(2) - A(2)) + x(:,2)*(C(2) - A(2));


% AX = x1 AB + x2 AC
% AX.AB = x1 ||AB||2 + x2 AB.AC
% AX.AC = x1 AB.AC + x2 ||AC||2
% L = K x
% x = K\L

end

