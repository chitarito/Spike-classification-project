function [W0,H0] = initialbasis(X,K,L,lambda,numinit) % give number of initial conditions to test
    Er = {};
    WsE = {};
    HsE = {};
    for i = 1:numinit
        [Ws,Hs] = seqNMF(X,'K',K, 'L', L,'lambda', lambda,'maxiter',5,'showplot',0);
        X_hat = helper.reconstruct(Ws,Hs);
        Er{i} = (sumsqr(X)-sumsqr(X_hat))/sumsqr(X);
        WsE{i} = Ws;
        HsE{i} = Hs;
    end
    ERM = cell2mat(Er);
    
    [RMS,I] = min(ERM);
    W0 = WsE{I};
    H0 = HsE{I};
    RMS
end