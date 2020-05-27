function [Optimal_Rank, PoV] = optrank(Data,fs,BFs,Bfs_of_Interest,Rank)
    LLA = Data;    
    % Percent of Variance explained for each rank
    PoVE = zeros(1,Rank);
    for j=1:Rank
        V = BFs(j).BWs*BFs(j).BHs;
        R = LLA-V;
        PVE = 1-sumsqr(R)/sumsqr(LLA);
        PoVE(j)=PVE;
    end

    % Partial explained variance using the activation vector of interest for
    % each rank
    pPoVE = zeros(1,Rank);
    for j=1:Rank
        voi = Bfs_of_Interest(j); % Get the index of the activation vector of interest for
        % the Basis vectors for every rank
        pHoI = BFs(j).BHs(voi,:); % Select the Activation vector of interest
        pWoI = BFs(j).BWs(:,voi); % Select corresponding Basis vector of interest
        pV = pWoI*pHoI;
        pR = LLA-pV;
        pPVE = 1-sumsqr(pR)/sumsqr(LLA);
        pPoVE(j) = pPVE;
    end

    % The final rank was determined by taking the first rank (in increasing order) 
    % that yielded an R2 PVE within 50% of the minimal R2 PVE for partial reconstructions.
    for i=1:Rank
        if abs(pPoVE(i)) < 1.5*min(pPoVE)
            fprintf('Optimal Rank: %i\n', i) 
            Optimal_Rank = i;
            break
        end
    end
    PoV = figure
    plot([1:Rank],PoVE);ylim([0,1]); title('Percent of Variance explained')
    hold on
    plot([1:Rank],pPoVE); title('Partial Percent of Variance explained')
    legend('Percent of Variance explained','Partial Explained Variance')
end