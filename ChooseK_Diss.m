function [DK,Ws,Hs,Diss] = ChooseK_Diss(data,K,nReps)
    % Dissimilarity score
    Data = data;
    Ws = {};
    Hs = {};
    numfits = nReps; %number of fits to compare
    tic
    for k = 1:K
        display(sprintf('running NMF with K = %i',k))
        parfor ii = 1:numfits
            opt = statset('MaxIter',5); %,'Display','final');
            [W0,H0] = nnmf(Data,k,'replicates',30,'options',opt,'algorithm','mult');
            [Ws{ii,k},Hs{ii,k},D] = nnmf(Data,k,'algorithm','als','w0',W0,'h0',H0);
        end
        
        inds = nchoosek(1:numfits,2);
        parfor i = 1:size(inds,1) % consider using parfor for larger numfits
            Diss(i,k) = nmfDISSX(Ws{inds(i,1),k},Ws{inds(i,2),k});
        end
    end
    toc
    % It took 3.3 hrs to run this code
    % Plot Diss and choose K with the minimum average diss.
    DK = figure
    plot(1:K,Diss,'ko'); title('Median Diss'); hold on
    h1 = plot(1:K,median(Diss,1),'k-','linewidth',2);
    legend(h1, {'mean Diss'}); xlabel('K'); ylabel('Diss')
end