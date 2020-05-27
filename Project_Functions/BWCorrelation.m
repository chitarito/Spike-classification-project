function [BWCorr, I] = BWCorrelation(W,WS,PPOVE) % Returns Correlation of each...
...BW to the extracted significant BWs and the BW pair that most correlates with them...
    ...as well as the percent of variance exlained by each basis in each set.

    % first set all of the basis with high PPOVE equal to zero
    for i = 1:length(WS)
        J = find(PPOVE(:,i) > 0.15);
        for j = 1:length(J)
            WS{i}(:,j,:) = 0;
        end
    end

    BWCorr = zeros(size(WS{1},2),size(W,2),length(WS)); % Correlation of each BW with the other BW from the basis sets
    for ii = 1:length(WS)
        Ws = cell2mat(WS(ii));
        for i = 1:size(W,2)
            for j = 1:size(Ws,2)
                w1 = squeeze(W(:,i,:)); w2 = squeeze(Ws(:,j,:));
                if sum(w2,'all') == 0
                    BWCorr(j,i,ii) = 0;
                else
                    C = xcorr2(w1,w2);
                    % only care about horizontal shift in correlation so,
                    BWCorr(j,i,ii) = max(C(size(Ws,1),:),[],'all'); 
                end
            end
        end
    end
    
    % find the basis that correlates most with the significant basis. This
    % is the most usefull decomposition.
    BWMaxCorr = zeros(length(WS),1);
    for i = 1:length(WS)
        BWMaxCorr(i) = sum(squeeze(BWCorr(:,:,i)),'all');
    end
    [m I] = max(BWMaxCorr);
    sprintf('Basis set %d correlates most with extracted factors.',I)
end