function [] = job(levels,Lap,Dir,Unif)


if Lap && Dir && Unif
    testName = "LaplaceS0Unif"
elseif Lap && Dir && ~Unif
    testName = "LaplaceS0NonUnif"
elseif Lap && ~Dir && Unif
    testName = "LaplaceH0Unif"
elseif Lap && ~Dir && ~Unif
    testName = "LaplaceH0NonUnif"
elseif ~Lap && Dir && Unif
    testName = "HelmholtzSkUnif"
elseif ~Lap && Dir && ~Unif
    testName = "HelmholtzSkNonUnif"
elseif ~Lap && ~Dir && Unif
    testName = "HelmholtzHkUnif"
elseif ~Lap && ~Dir && ~Unif
    testName = "HelmholtzHkNonUnif"
end



if Unif
    Ns = 5*4.^(levels);
else
    Ns = 2*4.^(levels+1);
end
c = 0.5;
RESTART = 20;
MAXIT = 10;
tolGMRES = 1e-6;
n1 = 0*levels;
t1 = 0*levels;
n2 = 0*levels;
t2 = 0*levels;
if ~Lap
    waveNum = 0*levels;
end
if ~Dir
    t3 = 0*levels;
    n3 = 0*levels;
end



for i = 1:length(Ns)
    
    'Clearing memory'
    
    clearvars -except i Ns c RESTART MAXIT tolGMRES...
        n1 t1 n2 t2 n3 t3 ...
        levels Lap Dir Unif waveNum
    
    
    'Assembling'
    
    tic;
    N = Ns(i);
    
    rad = 1;
    if Unif
        m = mshDisk(N,rad);
    else
        beta = 2;
        m = mshDiskGraded(N,rad,beta);
    end
    h = m.stp;
    if ~Lap
        if Unif
            waveNum(i) = 1/(2*h(2))
            k = waveNum(i);
        else
            waveNum(i) = 1/h(2)
            k = waveNum(i);
        end
    end
    
    Ngss = 3;
    D = dom(m,Ngss);
    
    if Dir
        Vh = fem(m,'P1');
    else
        Vh = dirichlet(fem(m,'P1'),m.bnd);
    end
    
    
    if Lap && Dir
        reg = regularize(D,D,Vh,'[1/r]',Vh);
        S = @(uh)(singleLayerlFMM(D,Vh,uh) + reg*uh); % Laplace Single layer
        u = integral(D,Vh,@(X)(1./omega(X))); % rhs
    elseif Lap && ~Dir
        reg1 = regularize(D,D,Vh,'[1/r]',Vh);
        reg2 = regularize(D,D,nxgrad(Vh),'[1/r]',nxgrad(Vh));
        S = @(uh)(singleLayerlFMM(D,Vh,uh) + reg1*uh); % Laplace Single layer
        H = @(uh)(hypersingularlFMM(D,Vh,uh) + reg2*uh);  %Laplace Hypersingular
        u = integral(D,Vh,@(X)(omega(X))); % rhs
    elseif ~Lap && Dir
        reg = regularize(D,D,Vh,'[1/r]',Vh);
        S = @(uh)(singleLayerhFMM(D,Vh,k,uh) + reg*uh); % Helmholtz Single layer
        PW = @(X)(exp(1i*k/sqrt(2)*(X(:,1) + X(:,3))));
        u = integral(D,Vh,PW); % rhs
    elseif ~Lap && ~Dir
        reg1 = regularize(D,D,Vh,'[1/r]',Vh);
        reg2 = regularize(D,D,nxgrad(Vh),'[1/r]',nxgrad(Vh)) ...
            - k^2*reg1;
        S = @(uh)(singleLayerhFMM(D,Vh,k,uh) + reg1*uh); % Helmholtz Single layer
        H = @(uh)(hypersingularhFMM(D,Vh,k,uh) + reg2*uh); % Helmtholtz hypersingular
        PW = @(X)(exp(1i*k/sqrt(2)*(X(:,1) + X(:,3))));
        gradPW{1} = @(X)(1i*k*sqrt(2)*PW(X));
        gradPW{2} = @(X)(0*X(:,2));
        gradPW{3} = @(X)(1i*k*sqrt(2)*PW(X));
        u = integral(D,ntimes(Vh),gradPW); % rhs
    end
    
    M = integral(D,Vh,Vh); % Mass matrix
    toc;
    
    'Weigthed quadratures'
    tic;
    if Dir
        [Mw_1,wDelta,Mw] = weightedMatDirichlet(Vh);
    else
        [Mw,Deltaw,Momega2w] = weightedMatNeumann(Vh);
    end
    toc;
    
    
    'Solving'
    
    % Preconditioners
    
    if Lap && Dir
        trefSq = @(x)(TrefethenSqrt1(wDelta+c*Mw_1,5,x,Mw_1,c,N^2));
        precSqrt = @(x)(M\(trefSq(M\x)));
    elseif Lap && ~Dir
        trefSq = @(x)(TrefethenSqrt2(Deltaw,5,x,Mw,1,N^2));
        precSqrt = @(x)(M\(trefSq(M\x)));
        precCald = @(x)(M\S(M\x));
    elseif ~Lap && Dir
        damping = -0.45i*k;
        K = wDelta - k^2*Mw + damping*Mw_1;
        K1 = K + k^2*Mw_1;
        Np = 15;
        theta = pi/3;
        [ C0,Aj,Bj ] = rationalCoeffSqrt(Np,theta,k);
        ratSq = @(x)(evalRat(Np,C0,Aj,Bj,Mw_1,K1,x));
        precSqrt = @(x)(M\(ratSq(M\x)));
    elseif ~Lap && ~Dir
        damping = -0.55i*k;
        Np = 15;
        theta = pi/3;
        [ C0,Aj,Bj ] = rationalCoeffSqrt(Np,theta,k);
        K = Deltaw - k^2*Momega2w + damping*Mw;
        K1 = K + k^2*Mw;
        ratSq = @(x)(evalRat(Np,C0,Aj,Bj,Mw,K1,x));
        precSqrt = @(x)(M\(Mw*(K\(ratSq(M\x)))));
        precCald = @(x)(M\S(M\x));
    end
    
    'No preconditioner'
    t11 = tic;
    if Dir
        [~,~,~,~,ITER1] = gmres(@(x)(S(x)),u,...
            RESTART,tolGMRES,MAXIT);
    else
        [~,~,~,~,ITER1] = gmres(@(x)(H(x)),u,...
            RESTART,tolGMRES,MAXIT);
    end
    t1(i) = toc(t11)
    n1(i) = length(ITER1) - 1
    
    
    parpool(4);
    
    'Square-root preconditioner'
    t22 = tic;
    v = precSqrt(u);
    if Dir
        [~,~,~,~,ITER2] = gmres(@(x)(precSqrt(S(x))),v,...
            RESTART,tolGMRES,MAXIT);
    else
        [~,~,~,~,ITER2] = gmres(@(x)(precSqrt(H(x))),v,...
            RESTART,tolGMRES,MAXIT);
    end
    t2(i) = toc(t22)
    n2(i) = length(ITER2)-1
    
    delete(gcp('nocreate'));
    
    if ~Dir
        'Calderon preconditioner'
        t33 = tic;
        [~,~,~,~,ITER3] = gmres(@(x)(precCald(H(x))),v,...
            RESTART,tolGMRES,MAXIT);
        t3(i) = toc(t33)
        n3(i) = length(ITER3)-1
    end
    
    
    
end


cd tables

if Lap && Dir && Unif
    testName = "LaplaceS0Unif";
elseif Lap && Dir && ~Unif
    testName = "LaplaceS0NonUnif";
elseif Lap && ~Dir && Unif
    testName = "LaplaceH0Unif";
elseif Lap && ~Dir && ~Unif
    testName = "LaplaceH0NonUnif";
elseif ~Lap && Dir && Unif
    testName = "HelmholtzSkUnif";
elseif ~Lap && Dir && ~Unif
    testName = "HelmholtzSkNonUnif";
elseif ~Lap && ~Dir && Unif
    testName = "HelmholtzHkUnif";
elseif ~Lap && ~Dir && ~Unif
    testName = "HelmholtzHkNonUnif";
end

fileID = fopen(strcat(testName,".tex"),"w");



if Lap && Dir
    fprintf(fileID,"\\begin{tabular}{c c c c c}\n");
    fprintf(fileID,"Refinement level  & $n_1$ & $t_1$ & $n_2$ & $t_2$\\\\\\hline\n");
elseif Lap && ~Dir
    fprintf(fileID,"\\begin{tabular}{c c c c c c c}\n");
    fprintf(fileID,"Refinement level  & $n_1$ & $t_1$ & $n_2$ & $t_2$ & $n_3$ & $t_3$ \\\\\\hline\n");
elseif ~Lap && Dir
    fprintf(fileID,"\\begin{tabular}{c c c c c c}\n");
    fprintf(fileID,"Refinement level & $k$ & $n_1$ & $t_1$ & $n_2$ & $t_2$\\\\\\hline\n");
elseif ~Lap && ~Dir
    fprintf(fileID,"\\begin{tabular}{c c c c c c c c}\n");
    fprintf(fileID,"Refinement level & $k$ & $n_1$ & $t_1$ & $n_2$ & $t_2$ & $n3$ & $t_3$ \\\\\\hline\n");
end



for i = 1:length(Ns)
    fprintf(fileID,'%s ',num2str(levels(i)));
    if ~Lap
        fprintf(fileID,'& %s',num2str(waveNum(i)));
    end
    fprintf(fileID,'& %s',num2str(n1(i)));
    fprintf(fileID,'& %.2f',t1(i));
    fprintf(fileID,'& %s',num2str(n2(i)));
    fprintf(fileID,'& %.2f',t2(i));
    if ~Dir
        fprintf(fileID,'& %s',num2str(n3(i)));
        fprintf(fileID,'& %.2f',t3(i));
    end
    if i ~= length(Ns)
        fprintf(fileID,"\\\\\n");
    end
end

fprintf(fileID,"\n");
fprintf(fileID,"\\end{tabular}");
cd ..



end

