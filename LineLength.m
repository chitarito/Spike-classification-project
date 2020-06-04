function [XX, LLR, LFPsmooth] = LineLength(LFP,fs)
% Compute linelength ratio of data (play around with long window)
% Smooth the amplitude of the absolute difference of the signal
    LFPsmooth = LFP;
    % Smooth the amplitude of the absolute value of the signal
    for i = 1:min(size(LFPsmooth)) % go across channels
        fsmooth = LFPsmooth(i,:);
        for j = 1:3
            fss = smoothdata(fsmooth,'movmean',3); fsmooth = fss;
        end
        LFPsmooth(i,:) = fsmooth;
    end

% Compute amplitude and difference of LFPsmooth
% Amplitude of the absolute value of the signal
    LLAa = abs(LFPsmooth);
    % Absolute value of the difference of the signal
    for i = 1:size(LFPsmooth,1) % for each channel
    % Compute the linelength transform of the channel signals
        f = abs([0 diff(LFPsmooth(i,:))]);
        LL = f;
        LLArray{i} = LL;
    end
    LFD = cell2mat(LLArray'); LFD = LFD(:,1:size(LFPsmooth,2));
% Smooth LFD (to get rid of zero in middle for when multiplying with amp)
    LFDsmooth = LFD;
    % Smooth the amplitude of the absolute value of the signal
    for i = 1:min(size(LFDsmooth)) % go across channels
        fsmooth = LFDsmooth(i,:);
        for j = 1:10
            fss = smoothdata(fsmooth,'movmean',3); fsmooth = fss;
        end
        LFDsmooth(i,:) = fsmooth;
    end
% Compute linelength
    XX = (1000*LLAa).*(1000*LFDsmooth);
    for i = 1:size(LFP,1) % for each channel
        % Compute the linelength transform of the channel signals
        w = round(0.02*fs); % 20ms (0.02s) window
        LL = conv(XX(i,:),rectwin(w)/w); 
        LLArray{i} = LL;
    end
    LLAss = cell2mat(LLArray'); LLAss = LLAss(:,1:size(LFP,2));

    % Compute long linelength
    for i = 1:size(LFP,1) % for each channel
        % Compute the linelength transform of the channel signals
        W = round(30*fs); % 30s window
        LL = conv(XX(i,:),rectwin(W)/W);  
        LLArray{i} = LL;
    end
    LLAl = cell2mat(LLArray'); LLAl = LLAl(:,1:size(LFP,2)); MM = mean(LLAl,2);
    % Mean the first W seconds
    LLAl(:,1:round(W)) = MM.*ones(size(LLAl,1),length(1:round(W)));
    LLR = LLAss./LLAl;
end