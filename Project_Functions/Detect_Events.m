function [Events_Detected,Event_indices,Mask] = Detect_Events(Data,Hs,Optimal_Rank) %His_PDF,
% This function requires the basis corresponding to the spike event
    LLA = Data;
    HoI = Hs;%(Optimal_Rank,:); % This is the basis that contains the spike

    HTH = zeros(size(HoI));
    for j = 1:size(Hs,1)
    % Use Kernel density estimation for PDF
        [f,xi,bw] = ksdensity(HoI(j,:)','Kernel','normal');
        % Smooth the PDF
        fsmooth = f;
        for i = 1%10
            fss = smoothdata(fsmooth,'movmean',3);
            fsmooth = fss;
        end
    % Find line that crosses through max point and max rate of change
        df = [0 diff(fsmooth)];[mdr, I] = min(df);[mfs, J] = max(fsmooth);
        % J = 14; % If there are miltiple peaks, set J at the last peak
        p = polyfit(xi(J:I),fsmooth(J:I),1); yp = p(1)*xi+p(2);
        % Threshold is where line crosses abcissa
        TH = find(yp < 0); Threshold = xi(TH(1));
                          
        ABC = HoI(j,:);
        ABC(ABC < 1*Threshold) = 0; % Threshold our activation vector
        HTH(j,:) = ABC;
        
        His_PDF = figure
        subplot(1,2,1)
        h = histogram(HoI(j,:),50); title('Histogram of H');
        counts = h.Values; xlabel('Amplitude of H'); ylabel('Count')
        subplot(1,2,2)
        plot(xi,f,'linewidth', 2); hold on; plot(xi,fsmooth,'linewidth', 2);
        plot(xi,yp,'linewidth', 2); title(sprintf('PDF of H%d',j));xlabel(sprintf('Amplitude of H%d',j));
        ylim([0 max(f)+100]); legend('Original Density','Smoothed Density', 'Line Fit')
    end
    HTH = sum(HTH,1);
                          
% Modify Z-score the data during timepoints outside of events
    NonZ = find(~HTH); Med = zeros(1,min(size(LLA))); STd = zeros(1,min(size(LLA)));
    
    for i = 1:min(size(LLA))
        VVV = LLA(i,NonZ);
        med = median(VVV);
        Std = std(VVV);
        Med(i) = med;
        STd(i) = Std;
    end
    % Find time points of events (the nonzero elements that remain)
    IsZ = find(HTH); LLAMZ = zeros(size(LLA)); 
% Maybe if the remaining parts of HTH are very near to each other, say like 10ms? combine into one event?
    for i = 1:min(size(LLA)) % go across channels
        for j = 1:length(IsZ) % go through IsZ for indicies in data
            if STd(i)==0
                LLAMZ(i,IsZ(j)) = 0;
            else
                LLAMZ(i,IsZ(j)) = (LLA(i,IsZ(j))-Med(i))/STd(i); % Z-score linelength data
            end
        end
    end
% Threshold modified Z-score channel
    LLAMZTH = LLAMZ; LLAMZTH(LLAMZTH < 6) = 0; % Important
    
    % Find time points of events
    k = [true;diff(IsZ(:))~=1]; % find the difference in IsZ where its not a difference of 1...
    ...indicating a new event. Start it with a true because the first part is also and event.
        s = cumsum(k); % counts the number of events
    x =  histc(s,1:s(end)); %counts the number of indices in each event
    idx = find(k); % Find the index of the beginning of each event in IsZ
    NumberofEvents = length(x);
    IVEvents = cell(1,length(x));
    for i = 1:length(x)
        index = find(s == i); % find the indices for each event in s
        IVEvents{i} = IsZ(index); % find actual indices in data from IsZ...
        ...now these indices are the indices for each event. Put them in cell array
    end
    % If during a contiguous detection, only one channel is ivolved, unlabel event
    HTHNew = HTH;
    for i = 1:length(IVEvents) % now go through each event
        Count = 0;
        for j = 1:min(size(LLA)) % go through the channels
            if sum(LLAMZTH(j,IVEvents{i})) > 0
                Count = Count+1; % if more than one channel is involved during the event increase the count
            end
        end
        if Count <= 1 % if only one channel is involved, declare no event
            HTHNew(:,IVEvents{i}) = 0;
        end
    end
%--------------------------------------------------------------------------
    % count the number of remaining events and return mask if nonzero
    ISZ = find(HTHNew); % Find indices in HTH where not equal to zero

    if length(ISZ) == 0
        Events_Detected = 0;
        Event_indices = 'N/A';
        Mask = 'N/A';
        disp('no events detected')
    else
        kk = [true;diff(ISZ(:))~=1 ]; S = cumsum(kk);
        xx =  histc(S,1:S(end)); IDX = find(kk);
    
        % Now index IsZ and that value is the beginning of the event in the data
        NoSE = length(xx); % Number of spike events
    
        % The indices for each event is as follows
        NOSE = cell(1,length(xx)); % Index vector for events
        for i = 1:length(xx)
            index = find(S == i); % Find indices where value is event i
            NOSE{i} = ISZ(index); % These are the indices for IsZ that are in event i
        end
    
        % Now we have the number of events detected NoSE, and the event indices NOSE
        Events_Detected = NoSE;
        Event_indices = NOSE;
        disp(sprintf('%d events detected',NoSE))
        
        % For each event wherever LLAMZTH is not zero, this is where an event
        % is taking place, label this with red. Create a mask.
        BLLAMZTH = logical(LLAMZTH);
        IIS = find(~HTHNew);
        for i = 1:min(size(LLA))
            for j = 1:length(Event_indices)
                BLLAMZTH(i,IIS) = 0;
            end
        end
        Mask = BLLAMZTH;
    end
end