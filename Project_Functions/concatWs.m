function [CWs,PPOVE] = concatWs(X,WS,HS) % Give it the cell of basis functions
    WWW = cell(size(WS));
    PPoVE = cell(size(WS)); % Get the partial percent of variance explained by each basis pair
    
    parfor i = 1:length(WS)
        K = size(HS{1},1);
        WW = cell(1,K);
        W = cell2mat(WS(i)); H = cell2mat(HS(i));
        pPoVE = zeros(size(H,1),1);
        for j = 1:K
            % Multiply Hs not of interest by zeros to calculate pPoVE
            Mask = zeros(size(H)); Mask(j,:) = ones(size(H(j,:)));
            pV = helper.reconstruct(W,(H.*Mask));
            pPVE = 1-sumsqr(X-pV)/sumsqr(X);
            pPoVE(j) = pPVE; % calculate variance explained by basis
            % if the variance is great it is not a spiking basis so set to
            % zero because it effects the correlations
            if pPVE > 0.15 % set threshold to explain less than 15%
                W(:,j,:) = 0;
            end
            WW{j} = squeeze(W(:,j,:));
        end
        PPoVE{i} = pPoVE;
        WWW{i} = cell2mat(WW);
    end
    PPOVE = cell2mat(PPoVE);
    CWs = cell2mat(WWW);
end