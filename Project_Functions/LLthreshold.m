function LLAN = LLthreshold(X,ZTH) % Give line length ratio array and Z_score threshold
    LLA = X;
    % Modify Z-score the data
    Med = zeros(min(size(LLA)),1); STd = zeros(min(size(LLA)),1);
    
    % Find the mean and standard deviation of each channel
    for i = 1:min(size(LLA))
        med = median(LLA(i,:)); Std = std(LLA(i,:));
        Med(i) = med; STd(i) = Std;
    end
    
    % Modify Z-score data
    for i = 1:min(size(LLA)) % go across channels
        if STd(i)==0
            LLAMZ(i,:) = zeros(1,length(LLA));
        else
            for j = 1:length(LLA) % go through IsZ for indicies in data
                LLAMZ(i,j) = (LLA(i,j)-Med(i))/STd(i); % Z-score linelength data
            end
        end
    end
    
    % Threshold modified Z-score channel
    LLAMZTH = LLAMZ; LLAMZTH(LLAMZTH < ZTH) = 0; % Important
    
    % If during a contiguous event, only one channel is involved, unlabel event
    LLAMZTHNew = LLAMZTH;
    for i = 1:min(size(LLA)) % go across channels
        % Find time points of events (the nonzero elements that remain)
        IsnZ = find(LLAMZTHNew(i,:));
        for j = 1:length(IsnZ) % go through IsZ
            if length(find(LLAMZTH(:,IsnZ(j)))) == 1
                % if only one channel is involved, declare no event
                LLAMZTHNew(i,IsnZ(j)) = 0;
            end
        end
    end
    
    mask = double(LLAMZTHNew > 0);
    LLAN = mask.*LLA;
end