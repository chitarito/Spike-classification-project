function CVplot = CrossValidation(X,K,L,nReps,nLambdas,lU,lL)
    lambdas = sort([logspace(lU,lL,nLambdas)], 'ascend'); %% Change range as needed %%

    RmseTrain = zeros(length(lambdas), nReps); 
    RmseTest = zeros(length(lambdas), nReps);
    Error = zeros(length(lambdas),nReps);
    
    CVplot = figure
    [~,Lplot] = meshgrid(1:nReps,lambdas);
    WS = cell(length(lambdas),nReps); HS = cell(length(lambdas),nReps);
    
    parfor li = 1:length(lambdas)
        for repi = 1:nReps
            display(['Cross validation on masked test set; Testing li = ' num2str(lambdas(li)) ', rep ' num2str(repi)])
            rng('shuffle')
            M = rand(N,T)>0.30; % create masking matrix (0's are test set, not used for fit)
            [W,H] = seqNMF(X,'K',K,'L',L,'lambda',lambdas(li),'maxiter',500,'showPlot',0, 'M', M);
            WS{li,repi} = W; HS{li,repi} = H;
            Xhat = helper.reconstruct(W,H);
            RmseTrain(li,repi) = sqrt(sum(M(:).*(X(:)-Xhat(:)).^2)./sum(M(:)));
            RmseTest(li,repi) = sqrt(sum((~M(:)).*(X(:)-Xhat(:)).^2)./sum(~M(:)));
        end
    end
    clf; scatter(Lplot(:),RmseTrain(:), 'r', 'markerfacecolor', 'flat');
    hold on; scatter(Lplot(:),RmseTest(:), 'b', 'markerfacecolor', 'flat');
    plot(lambdas,mean(RmseTrain,2), 'r'); plot(lambdas,mean(RmseTest,2), 'b')
    xlabel('Lambda'); ylabel('RMSE'); title('RMSE vs lambda')
    legend('Train','Test','location','northwest'); drawnow; shg
    
    %[m,I] = min(mean(RmseTrain,2));[n,J] = min(mean(RmseTest,2))
    %[l,O] = min((mean(RmseTrain,2)+mean(RmseTest,2))/2)
end
