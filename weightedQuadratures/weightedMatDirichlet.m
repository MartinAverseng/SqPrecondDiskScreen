function [Iw,Dw,Iw2] = weightedMatDirichlet(Vh)

m = Vh.msh;
dm = m.bnd;
[~,~,indBound] = intersect(dm.vtx,m.vtx,'rows','stable');

Nelt = size(m.elt,1);
Nvtx = size(m.vtx,1);
[Xel,T] = Vh.dof;
[x,w] = domReference(dom(m,12));
V = m.ndv;

eta{1} = @(x)(1 - x(:,1) - x(:,2));
dxeta(1) = -1;
dyeta(1) = -1;
eta{2} = @(x)(x(:,1));
dxeta(2) = 1;
dyeta(2) = 0;
eta{3} = @(x)(x(:,2));
dxeta(3) = 0;
dyeta(3) = 1;


I = [];
J = [];
valIw = [];
valDw = [];
valIw2 = [];
h = m.stp; h = h(1);

for el = 1:Nelt
    mask = ismember(m.elt(el,:),indBound);
    testBound = sum(mask);
    A = m.vtx(m.elt(el,1),:);
    B = m.vtx(m.elt(el,2),:);
    C = m.vtx(m.elt(el,3),:);
    switch testBound
        case 0            
            xel = x;
            Xel = triCoord(A,B,C,x);
            Wel = V(el)*w;
        case 1
            i1 = find(mask==1);
            i2 = mod(i1,3)+ 1;
            i3 = mod(i2,3)+ 1;
            a = m.vtx(m.elt(el,i1),:);
            b = m.vtx(m.elt(el,i2),:);
            c = m.vtx(m.elt(el,i3),:);
            [Xel,Wel] = intType2(a,b,c);
            xel = barycentricCoord(A,B,C,Xel);
        case 2
            i1 = find(mask~=1);
            i2 = mod(i1,3)+ 1;
            i3 = mod(i2,3)+ 1;
            a = m.vtx(m.elt(el,i2),:);
            b = m.vtx(m.elt(el,i3),:);
            c = m.vtx(m.elt(el,i1),:);
            [Xel,Wel] = intType3(a,b,c,6);
            xel = barycentricCoord(A,B,C,Xel);
    end
    E1 = B -A;
    E2 = C -A;
    a = sum(E1.^2);
    b = sum(E1.*E2);
    c = b;
    d = sum(E2.^2);
    
    % Determinant
    detG = a.*d - b.*c;
    
    % Inverse by cofactors
    Dx1 = d  ./ detG;
    Dx2 = -b ./ detG;
    Dy1 = -c ./ detG;
    Dy2 = a  ./ detG;
    
    
    for p = 1:3
        for q = 1:3
            nabp = (Dx1*E1 + Dx2*E2)*dxeta(p) + (Dy1*E1 + Dy2*E2)*dyeta(p);
            nabq = (Dx1*E1 + Dx2*E2)*dxeta(q) + (Dy1*E1 + Dy2*E2)*dyeta(q);
            nabpq = sum(nabp.*nabq);
            
            
            Apq = sum(Wel.*eta{p}(xel).*eta{q}(xel)./omega(Xel));
            Bpq = sum(Wel*nabpq.*omega(Xel));
            Cpq = sum(Wel.*eta{p}(xel).*eta{q}(xel).*omega(Xel));
            I = [I,T(el,p)]; %#ok
            J = [J,T(el,q)]; %#ok
            valIw = [valIw,Apq]; %#ok
            valDw = [valDw,Bpq]; %#ok            
            valIw2 = [valIw2,Cpq]; %#ok
        end
    end
    
end

Iw = sparse(I,J,valIw,Nvtx,Nvtx);
Dw = sparse(I,J,valDw,Nvtx,Nvtx);
Iw2 = sparse(I,J,valIw2,Nvtx,Nvtx);





% 
% 
% 
% dxF = zeros(dim,Ngss);
% dyF = zeros(dim,Ngss);
% 
% dxF(1,:) = - 1;
% dxF(2,:) = 1;
% 
% dyF(1,:) = - 1;
% dyF(3,:) = 1;
% 
% % Gradient projection to integration points
% dbas = cell(1,3);
% for n = 1:3
%     dbas{n} = zeros(Nelt,Nbas,Ngss);
%     DCVx    = Dx1.*E1(:,n) + Dx2.*E2(:,n);
%     DCVy    = Dy1.*E1(:,n) + Dy2.*E2(:,n);
%     for i = 1:Nbas
%         for j = 1:Ngss
%             dbas{n}(:,i,j) = DCVx(:)*dxF(i,j) + DCVy(:)*dyF(i,j);
%         end
%     end
% end
% 


