function CVKplot = CrossValidationK(X,K,L,nReps)
    Ks = [8:2:K]; %% Change range as needed %%
    [N,T] = size(X);
    RmseTrain = zeros(length(Ks), nReps); 
    RmseTest = zeros(length(Ks), nReps);
    
    CVKplot = figure;
    [~,Lplot] = meshgrid(1:nReps,Ks);
    WS = cell(length(Ks),nReps); HS = cell(length(Ks),nReps);
    
    parfor i = 1:length(Ks)
        for repi = 1:nReps
            display(['Cross validation on masked test set; Testing K = ' num2str(Ks(i)) ', rep ' num2str(repi)])
            rng('shuffle')
            M = rand(N,T)>0.30; % create masking matrix (0's are test set, not used for fit)
            [W,H] = seqNMF(X,'K',Ks(i),'L',L,'lambda',0,'maxiter',500,'showPlot',0, 'M', M);
            WS{i,repi} = W; HS{i,repi} = H;
            Xhat = helper.reconstruct(W,H);
            RmseTrain(i,repi) = sqrt(sum(M(:).*(X(:)-Xhat(:)).^2)./sum(M(:)));
            RmseTest(i,repi) = sqrt(sum((~M(:)).*(X(:)-Xhat(:)).^2)./sum(~M(:)));
        end
    end
    clf; scatter(Lplot(:),RmseTrain(:), 'r', 'markerfacecolor', 'flat');
    hold on; scatter(Lplot(:),RmseTest(:), 'b', 'markerfacecolor', 'flat');
    plot(Ks,mean(RmseTrain,2), 'r'); plot(Ks,mean(RmseTest,2), 'b')
    xlabel('K'); ylabel('RMSE'); title('RMSE vs K')
    legend('Train','Test','location','northwest'); drawnow; shg
    
    %[m,I] = min(mean(RmseTrain,2));[n,J] = min(mean(RmseTest,2))
    %[l,O] = min((mean(RmseTrain,2)+mean(RmseTest,2))/2)
end