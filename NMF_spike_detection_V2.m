% Ankit already processed data

%load EC206_Baseline_LLX.mat
load EC202_B9andB11_Alldata.mat % load time series data
%% Compute the linelength transform of data
LLA = LLR;
% Downsample it
d = 3;
dLLA = dec(LLA,d);



%% ----------------------- Procedure for choosing K -----------------------
%% Cross validation
K = 15;
nReps = 4;
CrossVal(dLLA,K,nReps)
%%
savefig('CVK_fullZ_K15_EC202_B9B11.fig')

%% Dissimilarity score
Data = dLLA;
Ws = {};
Hs = {};
numfits = 16; %number of fits to compare
tic
for k = 1:20
    display(sprintf('running seqNMF with K = %i',k))
    parfor ii = 1:numfits
        opt = statset('MaxIter',5); %,'Display','final');
        [W0,H0] = nnmf(Data,k,'replicates',30,'options',opt,'algorithm','mult');
        [Ws{ii,k},Hs{ii,k},D] = nnmf(Data,k,'algorithm','als','w0',W0,'h0',H0);
    end
    
    inds = nchoosek(1:numfits,2);
    parfor i = 1:size(inds,1) % consider using parfor for larger numfits
        Diss(i,k) = nmfDISSX(Ws{inds(i,1),k},Ws{inds(i,2),k});
    end
end
toc
% It took 3.3 hrs to run this code
%% Plot Diss and choose K with the minimum average diss.
figure
plot(1:20,DISS,'ko'); title('Mean Diss(corr) vs K: 16numfits K = 20'); hold on
h1 = plot(1:20,mean(DISS,1),'k-','linewidth',2);
legend(h1, {'mean Diss'})
xlabel('K')
ylabel('Diss')




%% Select K with lowest Diss score or lowest Cross Validation and run NMF

%takes a couple minutes to run on an hour of data sampled at 198Hz data
Data = dLLA;
K = 10; % In this case K = 4 seems to be the best option
% First, we obtained RMS of the reconstructions for 30 sets of randomly...
...initialized BFs and five iterations with multiplicative updates.  
opt = statset('MaxIter',5); %,'Display','final');
[W0,H0] = nnmf(Data,K,'replicates',30,'options',opt,'algorithm','mult');
% This will give us output the basis with the smallest RMS to use as...
...initial basis for the full decomposition
    
% Second, we obtained a full decomposition using alternating least squares...
...after initialization with the best set of bases found in the first step.
    OP = statset('Display','final');
[W,H,D] = nnmf(Data,K,'algorithm','als','w0',W0,'h0',H0,'options',OP);
%% Plot the Ws
figure
j = K;
imagesc(W);ylabel('Channels'); title('Basis Vectors');
xticks([1:j]); WV = num2cell(string(zeros(1,j)));
for i = 1:j
    y = sprintf('W%d', i); WV{i} = y;
end
xticklabels(WV);colormap(flipud(gray(256)));
%% Plot the Hs
figure
HV = num2cell(string(zeros(1,j)));
for i = 1:j
    x = sprintf('H%d', i); HV{i} = x;
end
mydata = H(:,(25000*l/d+1):25000*(l+2)/d);
t = linspace(0,length(mydata),length(mydata));
stackedPlot(t,HV,mydata);title('Activation Matrix')
xlabl = floor(length(mydata)*d/fs);X=[0:xlabl];string(X);
xticks([0:(fs/d):length(mydata)]);xticklabels({X});xlabel('time (s)');
%% Identify H of interest for each rank visually? How to pick spikes in this?
% Plot partial reconstructions
l = 0;
for j=1:K
    figure
    Hs = zeros(size(H)); Hs(j,:) = H(j,:);
    pV = W*Hs;
    mydata = [1:size(pV,1)]' + pV(:,round(25000*l/d+1):round((25000*l+25000)/d))*0.1;
    t = 1:length(mydata);
    plot(t,mydata);
    title(sprintf('Reconstruction with basis %d',j))
    xlabl1 = floor(length(mydata)*d/fs);X1 = [0:xlabl1];string(X1);
    xticks([0:fs/d:length(mydata)]); xticklabels({X1}); xlabel('time (s)')
    grid on
end
%% Compare reconstructions
d = 3;
dLFP = dec(LFP,d);
%%
l = 0;
mydata = [1:size(dLFP,1)]' + dLLA(:,round(25000*l/d+1):round((25000*l+25000)/d))*0.05;
t = 1:length(mydata);
figure
plot(t,mydata,'b'); title('dLLA')
xlabl1 = floor(length(mydata)*d/fs);X1 = [0:xlabl1];string(X1);
xticks([0:fs/d:length(mydata)]); xticklabels({X1}); xlabel('time (s)')
grid on

mydata = [1:size(LFP,1)]' + W*H(:,round(25000*l/d+1):round((25000*l+25000)/d))*0.05;
t = 1:length(mydata);
figure
plot(t,mydata,'b'); title('W*H')
xlabl1 = floor(length(mydata)*d/fs);X1 = [0:xlabl1];string(X1);
xticks([0:fs/d:length(mydata)]); xticklabels({X1}); xlabel('time (s)')
grid on
%%
1-sumsqr(W*H)/sumsqr(dLLA)
%%
l = 0;

figure
subplot('Position',[0.1 0.4 0.2 0.55])
j = K; imagesc(flipud(W));ylabel('Channels'); title('Basis Vectors');
xticks([1:j]); WV = num2cell(string(zeros(1,j)));
for i = 1:j
    y = sprintf('W%d', i); WV{i} = y;
end
xticklabels(WV);colormap(flipud(gray(256)));

a1 = subplot('Position',[0.4 0.4 0.55 0.55])
mydata = [1:size(LFP,1)]' + dLFP(:,(25000*l+1):(25000*l+25000))*3000;
t = 1:length(mydata); plot(t,mydata,'b'); title('dLLA')
xlabl1 = floor(length(mydata)*d/fs);X1 = [0:xlabl1];string(X1);
xticks([0:fs/d:length(mydata)]); xticklabels({X1}); xlabel('time (s)')
grid on; xlim([0 length(t)])

a2 = subplot('Position',[0.4 0.1 0.55 0.2])
mydata = [1:size(H,1)]' + H(:,(25000*l+1):(25000*l+25000))*70;
t = 1:length(mydata); plot(t,mydata,'b'); title('H')
xlabl1 = floor(length(mydata)*d/fs);X1 = [0:xlabl1];string(X1);
xticks([0:fs/d:length(mydata)]); xticklabels({X1}); xlabel('time (s)')
grid on; xlim([0 length(t)])

linkaxes([a1,a2], 'x');
%% Maybe use all activation vectors and threshold them
% Lets take the first basis and threshold it to detect spikes

% Plot distribution of data and PDFs and detect events
%[His_PDF,Events_Detected,Event_indices,Mask] = Detect_Events(Data,Hs,Optimal_Rank)
[Events_Detected,Event_indices,Mask] = Detect_Events(dLLA,H,1);

%% Plot the Ws
figure
j = K;
imagesc(W);ylabel('Channels'); title('Basis Vectors');
xticks([1:j]); WV = num2cell(string(zeros(1,j)));
for i = 1:j
    y = sprintf('W%d', i); WV{i} = y;
end
xticklabels(WV);colorbar
%% Threshold W
WW = W;
stdw = std(W,0,1);
medw = median(W,1);
for i = 1:size(W,2)
    AB = (W(:,i)-medw(i))/stdw(i);
    AB(AB < 1) = 0; % Threshold our activation vector
    WW(:,i) = W(:,i).*double(logical(AB));
end

figure
j = K;
imagesc(WW);ylabel('Channels'); title('Basis Vectors');
xticks([1:j]); WV = num2cell(string(zeros(1,j)));
for i = 1:j
    y = sprintf('W%d', i); WV{i} = y;
end
xticklabels(WV);colorbar
%% If no events, do not proceed
M = [1:10];
HoI = H(M,:);
HTH = zeros(size(H));
for j = 1:size(HoI,1)
    % Use Kernel density estimation for PDF
    [f,xi,bw] = ksdensity(HoI(j,:)','Kernel','normal');
    % Smooth the PDF
    fsmooth = f;
    for i = 1%:10
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
    ABC(ABC < 0.25*Threshold) = 0; % Threshold our activation vector
    HTH(M(j),:) = ABC;
    
%     His_PDF = figure
%     subplot(1,2,1)
%     h = histogram(HoI(j,:),50); title('Histogram of H');
%     counts = h.Values; xlabel('Amplitude of H'); ylabel('Count')
%     subplot(1,2,2)
%     plot(xi,f,'linewidth', 2); hold on; plot(xi,fsmooth,'linewidth', 2);
%     plot(xi,yp,'linewidth', 2); title(sprintf('PDF of H%d',M(j)));xlabel(sprintf('Amplitude of H%d',M(j)));
%     ylim([0 max(f)+100]); legend('Original Density','Smoothed Density', 'Line Fit')
end

% M = [10];
% HoI = H(M,:);
% for j = 1:size(HoI,1)
%     % Use Kernel density estimation for PDF
%     [f,xi,bw] = ksdensity(HoI(j,:)','Kernel','normal');
%     % Smooth the PDF
%     fsmooth = f;
% %     for i = 1:10
% %         fss = smoothdata(fsmooth,'movmean',3);
% %         fsmooth = fss;
% %     end
%     % Find line that crosses through max point and max rate of change
%     df = [0 diff(fsmooth)];[mdr, I] = min(df);[mfs, J] = max(fsmooth);
%     % J = 14; % If there are miltiple peaks, set J at the last peak
%     p = polyfit(xi(J:I),fsmooth(J:I),1); yp = p(1)*xi+p(2);
%     % Threshold is where line crosses abcissa
%     TH = find(yp < 0); Threshold = xi(TH(1));
%     
%     ABC = HoI(j,:);
%     ABC(ABC < 0.25*Threshold) = 0; % Threshold our activation vector
%     HTH(M(j),:) = ABC;
%     
% %     His_PDF = figure
% %     subplot(1,2,1)
% %     h = histogram(HoI(j,:),50); title('Histogram of H');
% %     counts = h.Values; xlabel('Amplitude of H'); ylabel('Count')
% %     subplot(1,2,2)
% %     plot(xi,f,'linewidth', 2); hold on; plot(xi,fsmooth,'linewidth', 2);
% %     plot(xi,yp,'linewidth', 2); title(sprintf('PDF of H%d',M(j)));xlabel(sprintf('Amplitude of H%d',M(j)));
% %     ylim([0 max(f)+100]); legend('Original Density','Smoothed Density', 'Line Fit')
% end
OM = double(logical(WW*HTH));
HTH = sum(HTH,1);
%
LLA = dLLA;
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
LLAMZTH = LLAMZ; LLAMZTH(LLAMZTH < 20) = 0; % Important

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
    BLLAMZTH = LLAMZTH; % Make thresholded matrix a logical
    for i = 1:min(size(LLAMZTH))
        BLLAMZTH(i,:) = LLAMZTH(i,:).*HTHNew;
    end
    Mask = logical(BLLAMZTH);
end
%%
imagesc(LLAMZTH);colorbar
%%
l = 0;
LabelEvents(dLFP(:,(25000*l/d+1):25000*(l+2)/d),Mask(:,(25000*l/d+1):25000*(l+2)/d),round(fs/d),2.1);
title('NMF Detected Spikes') % Only works with LLA for some reason

LabelEvents(dLLA(:,(25000*l/d+1):25000*(l+2)/d),Mask(:,(25000*l/d+1):25000*(l+2)/d),round(fs/d),0.0001);
title('NMF Detected Spikes') % Only works with LLA for some reason
%%
mask = Mask.*OM;
l = 0;
LabelEvents(dLFP(:,(25000*l/d+1):25000*(l+2)/d),mask(:,(25000*l/d+1):25000*(l+2)/d),round(fs/d),2.2);
title('NMF Detected Spikes') % Only works with LLA for some reason

LabelEvents(dLLA(:,(25000*l/d+1):25000*(l+2)/d),mask(:,(25000*l/d+1):25000*(l+2)/d),round(fs/d),0.0001);
title('NMF Detected Spikes') % Only works with LLA for some reason