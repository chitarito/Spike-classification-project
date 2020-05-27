load EC202_B9andB11_Alldata.mat % load time series data
%% Compute linelength ratio of data (play around with long window)
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
%
mydata = [1:size(LFP,1)]' + LFP(:,71000:78000)*2000;
t = 1:length(mydata);
figure
plot(t,mydata,'b'); title('LFP')
xlabl1 = floor(length(mydata)/fs);X1 = [0:xlabl1];string(X1);
xticks([0:fs:length(mydata)]); xticklabels({X1}); xlabel('time (s)')
grid on

mydata = [1:size(LFP,1)]' + LFPsmooth(:,71000:78000)*2000;
t = 1:length(mydata);
figure
plot(t,mydata,'b'); title('LFPsmooth')
xlabl1 = floor(length(mydata)/fs);X1 = [0:xlabl1];string(X1);
xticks([0:fs:length(mydata)]); xticklabels({X1}); xlabel('time (s)')
grid on
 
%% Compute amplitude and difference of LFPsmooth
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
%%
mydata = [1:size(LFP,1)]' + LFPsmooth(:,71000:78000)*4000;
t = 1:length(mydata);
figure
plot(t,mydata,'b')
xlabl1 = floor(length(mydata)/fs);X1 = [0:xlabl1];string(X1);
xticks([0:(fs):length(mydata)]); xticklabels({X1}); xlabel('time (s)')
grid on

mydata = [1:size(LFP,1)]' + LFD(:,71000:78000)*2*10^4;
t = 1:length(mydata);
figure
plot(t,mydata,'b');title('Difference')
xlabl1 = floor(length(mydata)/fs);X1 = [0:xlabl1];string(X1);
xticks([0:(fs):length(mydata)]); xticklabels({X1}); xlabel('time (s)')
grid on

mydata = [1:size(LFP,1)]' + LLAa(:,71000:78000)*0.3*10^4;
t = 1:length(mydata);
figure
plot(t,mydata,'b');title('Amplitude')
xlabl1 = floor(length(mydata)/fs);X1 = [0:xlabl1];string(X1);
xticks([0:(fs):length(mydata)]); xticklabels({X1}); xlabel('time (s)')
grid on
%% Smooth LFD (to get rid of zero in middle for when multiplying with amp)
LFDsmooth = LFD;
 % Smooth the amplitude of the absolute value of the signal
 for i = 1:min(size(LFDsmooth)) % go across channels
    fsmooth = LFDsmooth(i,:);
    for j = 1:10
        fss = smoothdata(fsmooth,'movmean',3); fsmooth = fss;
    end
    LFDsmooth(i,:) = fsmooth;
 end
%%
mydata = [1:size(LFP,1)]' + LFD(:,71000:78000)*2*10^4;
t = 1:length(mydata);
figure
plot(t,mydata,'b')
xlabl1 = floor(length(mydata)/fs);X1 = [0:xlabl1];string(X1);
xticks([0:(fs):length(mydata)]); xticklabels({X1}); xlabel('time (s)')
grid on

mydata = [1:size(LFP,1)]' + LFDsmooth(:,71000:78000)*2*10^4;
t = 1:length(mydata);
figure
plot(t,mydata,'b');title('Difference smoothed')
xlabl1 = floor(length(mydata)/fs);X1 = [0:xlabl1];string(X1);
xticks([0:(fs):length(mydata)]); xticklabels({X1}); xlabel('time (s)')
grid on
%% Multiply both together
XX = (LLAa*10^6).*(LFDsmooth*10^6);
%%
mydata = [1:size(LFP,1)]' + LFPsmooth(:,71000:78000)*3000;
t = 1:length(mydata);
figure
plot(t,mydata,'b')
xlabl1 = floor(length(mydata)/fs);X1 = [0:xlabl1];string(X1);
xticks([0:(fs):length(mydata)]); xticklabels({X1}); xlabel('time (s)')
grid on

mydata = [1:size(LFP,1)]' + XX(:,71000:78000)*0.00004;
t = 1:length(mydata);
figure
plot(t,mydata,'b');title('Difference.*Amplitude')
xlabl1 = floor(length(mydata)/fs);X1 = [0:xlabl1];string(X1);
xticks([0:(fs):length(mydata)]); xticklabels({X1}); xlabel('time (s)')
grid on
%% Compute linelength
for i = 1:size(LFP,1) % for each channel
    % Compute the linelength transform of the channel signals
    w = round(0.025*fs); % 25ms (0.025s) window
    LL = conv(XX(i,:),rectwin(w)/w); 
    LLArray{i} = LL;
end
LLAss = cell2mat(LLArray'); LLAss = LLAss(:,1:size(LFP,2));
%%
mydata = [1:size(LFP,1)]' + XX(:,71000:78000)*0.00004;
t = 1:length(mydata);
figure
plot(t,mydata,'b')
xlabl1 = floor(length(mydata)/fs);X1 = [0:xlabl1];string(X1);
xticks([0:(fs):length(mydata)]); xticklabels({X1}); xlabel('time (s)')
grid on

mydata = [1:size(LFP,1)]' + LLAss(:,71000:78000)*0.00006;
t = 1:length(mydata);
figure
plot(t,mydata,'b');title('25ms Line Length')
xlabl1 = floor(length(mydata)/fs);X1 = [0:xlabl1];string(X1);
xticks([0:(fs):length(mydata)]); xticklabels({X1}); xlabel('time (s)')
grid on
%%
% Compute long linelength
for i = 1:size(LFP,1) % for each channel
    % Compute the linelength transform of the channel signals
    W = round(30*fs); % 40ms (0.040s) window
    LL = conv(XX(i,:),rectwin(W)/W);  
    LLArray{i} = LL;
end
LLAl = cell2mat(LLArray'); LLAl = LLAl(:,1:size(LFP,2)); MM = mean(LLAl,2);
% Mean the first W seconds
LLAl(:,1:round(W)) = MM.*ones(size(LLAl,1),length(1:round(W)));
LLR = LLAss./LLAl;
%%
mydata = [1:size(LFP,1)]' + LLR(:,71000:78000)*0.04;
t = 1:length(mydata);
figure
plot(t,mydata,'b');title('30s Line Length Ratio')
xlabl1 = floor(length(mydata)/fs);X1 = [0:xlabl1];string(X1);
xticks([0:(fs):length(mydata)]); xticklabels({X1}); xlabel('time (s)')
grid on

mydata = [1:size(LFP,1)]' + LLAss(:,71000:78000)*0.00006;
t = 1:length(mydata);
figure
plot(t,mydata,'b');title('25ms Line Length')
xlabl1 = floor(length(mydata)/fs);X1 = [0:xlabl1];string(X1);
xticks([0:(fs):length(mydata)]); xticklabels({X1}); xlabel('time (s)')
grid on
%%
mydata = [1:size(LFP,1)]' + LFP(:,71000:78000)*3000;
t = 1:length(mydata);
figure
plot(t,mydata,'b');title('LFP')
xlabl1 = floor(length(mydata)/fs);X1 = [0:xlabl1];string(X1);
xticks([0:(fs):length(mydata)]); xticklabels({X1}); xlabel('time (s)')
grid on