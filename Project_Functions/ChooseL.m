function Lcosts = ChooseL(X,K,L,nL)
    nLambdas = nL; % increase if you're patient
    lambdas = sort([logspace(-3,-6,nLambdas)], 'ascend'); %% Change range as needed %%
    % Would it not be more accurate to run multiple iterations per lambda
    % value?
    loadings = [];
    regularization = [];
    cost = [];
    
    parfor li = 1:length(lambdas)
        display(['Testing lambda ' num2str(li) '/' num2str(length(lambdas))])
        [N,T] = size(X);
        [W, H, ~,loadings(li,:),power]= seqNMF(X,'K',K,'L',L,'lambda', lambdas(li), 'maxiter', 500, 'showPlot', 0);
        [cost(li),regularization(li),~] = helper.get_seqNMF_cost(X,W,H);
    end
    
    Lcosts = figure;
    
    windowSize = 3;
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    Rs = filtfilt(b,a,regularization);
    minRs = prctile(regularization,10); maxRs= prctile(regularization,90);
    Rs = (Rs-minRs)/(maxRs-minRs);
    R = (regularization-minRs)/(maxRs-minRs);
    Cs = filtfilt(b,a,cost);
    minCs =  prctile(cost,10); maxCs =  prctile(cost,90);
    Cs = (Cs -minCs)/(maxCs-minCs);
    C = (cost -minCs)/(maxCs-minCs);
    
    clf; hold on
    plot(lambdas,Rs, 'b')
    plot(lambdas,Cs,'r')
    scatter(lambdas, R, 'b', 'markerfacecolor', 'flat');
    scatter(lambdas, C, 'r', 'markerfacecolor', 'flat');
    xlabel('Lambda'); ylabel('Cost (au)')
    set(legend('Correlation cost', 'Reconstruction cost'), 'Box', 'on')
    set(gca, 'xscale', 'log', 'ytick', [], 'color', 'none')
    set(gca,'color','none','tickdir','out','ticklength', [0.025, 0.025])

end