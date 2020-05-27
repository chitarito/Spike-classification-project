function LLR = linlenR(LFP,fs,w,W) % (w and W in seconds)
    LLAshort = linelength(LFP,fs,w); % w length in seconds
    LLAlong = linelength(LFP,fs,W); % W length in seconds
    LLAlong(:,1:round(W*fs)) = 1; 
    % if the LLR is less zero, make the LLr zero
    LLr = (LLAshort)./LLAlong; MM = mean(LLr,2); LLR = LLr;
    LLR(:,1:round(W*fs)) = MM.*ones(size(LLr,1),length(1:round(W*fs))); % Mean the first W seconds
end

