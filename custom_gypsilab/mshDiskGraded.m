function [mesh] = mshDiskGraded(N,rad,beta)

e = exp(1);
if beta == 2
    g = @(x)(x/lambertw(x));
    Y = g(beta*N/(2*pi*e));
    M = ceil(sqrt(e*Y));
elseif beta > 2
    M = ceil(((beta*N)/(2*pi*zeta(beta-1)))^(1/beta));
else
    M = ceil(((2-beta)*N/(pi*beta))^(1/2));
end

r = fliplr(1 - linspace(0,1,M).^beta);
r = r(2:end);
rs = [0];
ts = [0];
for i = 1:M-1
    if i == 1
        dtheta = 1;
    else
        dtheta = (r(i) - r(i-1))/r(i);
    end
    L = ceil(2*pi/dtheta);
    thet = linspace(0,2*pi,L);
    thet = thet(1:end-1);
    ts = [ts, thet];
    rs = [rs,0*thet + r(i)];
end

DT = delaunayTriangulation(rs(:).*cos(ts(:)),rs(:).*sin(ts(:)));

elt  = DT.ConnectivityList;
vtx  = [DT.Points,zeros(size(DT.Points,1),1)];
mesh = msh(rad*vtx,elt);


end

