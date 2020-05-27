function [DissK,Wss,Hss] = dissk(K,numfits)
    % Give range of K and number of fits to compare
    Wss = {}; Hss = {};
    for k = 1:K
        display(sprintf('running seqNMF with K = %i',k))
        parfor ii = 1:numfits
            [Wss{ii,k},Hss{ii,k}] = seqNMF(WWE,'K',k, 'L', L,'lambda', 0,'maxiter',500,'showplot',0);
            % note that max iter set low (30iter) for speed in demo (not recommended in practice)
        end
        inds = nchoosek(1:numfits,2);
        for i = 1:size(inds,1) % consider using parfor for larger numfits
            Diss(i,k) = helper.DISSX(Hss{inds(i,1),k},Wss{inds(i,1),k},Hss{inds(i,2),k},Wss{inds(i,2),k});
        end
    end
    DissK = Diss;
end