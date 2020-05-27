function LLA = linelength(Signal,fs,T) % W length in seconds
    LLA = abs(Signal);
    LFP = Signal;
    % Smooth the amplitude of the absolute value of the signal
    LLA1 = LLA;
    for i = 1:min(size(LLA1)) % go across channels
        fsmooth = LLA1(i,:);
        for j = 1:20
            fss = smoothdata(fsmooth,'movmean',3); fsmooth = fss;
        end
        LLA1(i,:) = fsmooth;
    end
%  %  
%     LLA = LLA1;
%     % Modify Z-score the data
%     Med = zeros(min(size(LLA)),1); STd = zeros(min(size(LLA)),1);
%     % Find the mean and standard deviation of each channel
%     for i = 1:min(size(LLA))
%         med = mean(LLA(i,:)); Std = std(LLA(i,:));
%         Med(i) = med; STd(i) = Std;
%     end
%     % Modify Z-score data
%     for i = 1:min(size(LLA)) % go across channels
%         if STd(i)==0
%             LLAMZ(i,:) = zeros(1,length(LLA));
%         else
%             for j = 1:length(LLA) % go through IsZ for indicies in data
%                 LLAMZ(i,j) = (LLA(i,j)-Med(i))/STd(i); % Z-score linelength data
%             end
%         end
%     end
% %  
    LLArray = {};
    % Smooth the amplitude of the absolute value of the difference of the signal
    for i = 1:size(LFP,1) % for each channel
        % Compute the linelength transform of the channel signals
        f = abs([0 diff(LFP(i,:))]);
        LL = f;
        LLArray{i} = LL;
    end
    LFD = cell2mat(LLArray'); LFD = LFD(:,1:size(LFP,2));

    % Smooth
    LLAd = LFD;
    % Smooth the amplitude of the absolute value of the signal
    for i = 1:min(size(LLAd)) % go across channels
        fsmooth = LLAd(i,:);
        for j = 1:20
            fss = smoothdata(fsmooth,'movmean',3); fsmooth = fss;
        end
        LLAd(i,:) = fsmooth;
    end
    
 %
    LLArray = {};
    for i = 1:size(Signal,1) % for each channel
        % Compute the linelength transform of the channel signals
        W = round(T*fs); % 40ms (0.040s) window
        
        f = LLAd(i,:).*LLA1(i,:);
        
        LL = conv(f,rectwin(W)/W);
        LLArray{i} = LL;
    end
    LLA = cell2mat(LLArray'); LLA = LLA(:,1:size(Signal,2));
    display(sprintf('LLA shape: %dx%d' ,size(LLA,1),size(LLA,2)))
end