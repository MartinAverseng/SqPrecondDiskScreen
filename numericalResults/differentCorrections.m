


close all;
clear all;
delete(gcp('nocreate'));

% Switch to choose type of mesh
Unif = 1; %1 

levels = 1:3;
if Unif 
    Ns = 5*4.^(levels);
else
    Ns = 2*4.^(levels+1);
end


waveNum = 0*Ns;
n2 = 0*Ns;
n3 = 0*Ns;
n4 = 0*Ns;


for i = 1:length(Ns)
    
    'Cleaning memory'
    clearvars -except n1 n2 n3 n4...
        Nvt c i...
        levels Ns waveNum Unif...
        
    
    
    N = Ns(i);
    rad = 1;
    if Unif
        m = mshDisk(N,rad);
    else
        beta = 2;
        m = mshDiskGraded(N,rad,beta);
    end
    h = m.stp;
    
    if Unif        
        waveNum(i) = 1/(2*h(2));
    else
        waveNum(i) = 1/(h(2));
    end
    k = waveNum(i);
    
    D = dom(m,3);
    Vh = fem(m,'P1');
    clear m
    
    
    'Assembling'
    tic;
    reg = regularize(D,D,Vh,'[1/r]',Vh);
    M = integral(D,Vh,Vh);
    PW = @(X)(exp(1i*k/sqrt(2)*(X(:,1) + X(:,3))));
    u = integral(D,Vh,PW);
    toc;
    
    'Weighted matrices'
    
    tic;
    [Mw_1,wDelta,Mw] = weightedMatDirichlet(Vh);
    toc;
    
    G = @(X,Y)femGreenKernel(X,Y,'[exp(ikr)/r]',k);
    S =@(uh)(singleLayerhFMM(D,Vh,k,uh) + reg*uh);
    
    clear reg
    
    Np = 15;
    theta = pi/3;
    
    'solving'
    
    [ C0,Aj,Bj ] = rationalCoeffSqrt(Np,theta,k);
    
    
    
    RESTART = 20;
    MAXIT = 10;
    
    parpool(4);
    
    'square root prec for k = 0'
    K = wDelta + Mw;
    trefSq = @(x)(TrefethenSqrt1(K,5,x,Mw_1,0.5,N^2));
    Prec = @(x)(M\(trefSq(M\x)));
    v = Prec(u);
    [~,~,~,~,ITER] = gmres(@(x)Prec(S(x)),v,RESTART,1e-6,MAXIT);
    n2(i) = length(ITER) - 1
    clear K v
    
    'Naive k depedency'
    damping = -0.45*1i*k;
    K = wDelta + damping*Mw_1;
    ratSq = @(x)(evalRat(Np,C0,Aj,Bj,Mw_1,K,x));
    Prec = @(x)(M\(ratSq(M\x)));
    v = Prec(u);
    [~,~,~,~,ITER] = gmres(@(x)(Prec(S(x))),v,RESTART,1e-6,MAXIT);
    n3(i) = length(ITER) - 1
    clear wDelta Mw Mw_1 K v
    
    'pseudo-diff DtN approx'
    dM = integral(D,grad(Vh),grad(Vh));
    K = dM*k^2/(k+1i*0.4*k^1/3)^2;
    ratSq = @(x)(evalRat(Np,C0,Aj,Bj,M,K,x));
    Prec = @(x)(M\(ratSq(M\x)));
    v = Prec(u);
    [~,~,~,~,ITER4] = gmres(@(x)Prec(S(x)),v,RESTART,1e-6,MAXIT);
    n4(i) = length(ITER4) - 1
    clear K dM v
    
    
    delete(gcp('nocreate'));
end
