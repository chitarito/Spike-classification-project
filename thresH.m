function [mask,HsTH,OM] = thresH(dLLA,Th1,Th2,H,WW)
    K = size(H,1);
    % Threshold  activation vectors based on their interspike mean
    % Modify Z-score the data during timepoints outside of events
    Med = median(dLLA,2); STD = std(dLLA,0,2); LLAMZ = zeros(size(dLLA));
    %
    for i = 1:min(size(dLLA)) % go across channels
        if STD(i) == 0
            LLAMZ(i,:) = 0;
        else
            LLAMZ(i,:) = (dLLA(i,:)-Med(i))/STD(i); % Z-score linelength data
        end
    end
    % Threshold modified Z-score channel
    LLAMZTH = LLAMZ; LLAMZTH(LLAMZTH > Th1) = 0;
    HiAmp = double(logical(sum(LLAMZTH,1)));
    %
    HThresh = zeros(1,K);
    Hs = zeros(size(H));
    for j = 1:K
        % Set high amp event times to zero
        X = H(j,:).*HiAmp;
        % Find the mean of the remaining non zero values for each basis
%         I = find(X); HThresh(j) = mean(X(I));
%         X(X < mean(X(I))) = 0; %% maybe we should also use the standard deviation
        Hs(j,:) = X;
    end
    % Now we have the thresholded basis Hs
    HsTH = Hs;
    OM = double(logical(WW*HsTH));
  
    % If a high amp event only occurs across a single channel, unlabel event
    LlamZTH = LLAMZ; LlamZTH(LlamZTH < Th2) = 0;
    
    Mask = double(logical(LlamZTH));
    mask = OM.*Mask;
end