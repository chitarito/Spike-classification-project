% Ankit already processed data
load EC202_B9andB11_Alldata.mat % load in preprocessed time series data
%% Compute the linelength transform of data
[XX, LLR, LFPsmooth] = LineLength(LFP,fs);
% Downsample by d
d = 3; dLLA = dec(LLR,d); dLFP = dec(LFP,d);
dLFPs = dec(LFPsmooth,d); dsig = dec(XX,d);


%% ----------------------- Procedure for choosing K -----------------------
%% Cross validation
K = 30; nReps = 8;
CrossVal(dLLA,K,nReps)
%%
savefig('CVK_K30_EC202_B9B11.fig')
%% Choose K with Dissimilarity score
[DK,Ws,Hs,Diss] = ChooseK_Diss(dLLA,K,nReps);
%% can change to plot median or mean
figure
plot(1:K,Diss,'ko'); title('Mean Diss'); hold on
h1 = plot(1:K,mean(Diss,1),'k-','linewidth',2);
legend(h1, {'mean Diss'}); xlabel('K'); ylabel('Diss')


%% Select K with lowest Diss score or lowest Cross Validation and run NMF
%takes a couple minutes to run on an hour of data sampled at (fs/d)Hz data
%[W,H,D] = NMFdecomp(data,K)
K = 4; [W,H,D] = NMFdecomp(dLLA,K);
%%% sumsqr error
%error = 1-sumsqr(W*H)/sumsqr(dLLA)

%% Threshold W
WW = W; stdw = std(W,0,1); medw = median(W,1);
for i = 1:size(W,2)
    AB = (W(:,i)-medw(i))/stdw(i); AB(AB < 0.5) = 0; % Threshold basis
    WW(:,i) = W(:,i).*double(logical(AB));
end

%% Threshold  activation vectors based on their interspike mean
% [mask,HsTH,OM] = thresH(dLLA,Th1,Th2,H,WW)
% [mask,HsTH,OM] = thresH(dLLA,1,1,H,W);
[mask,HsTH,OM] = thresH(dLLA,1,1,H,WW); % Don't thresh W, and lower H thresh
%% Now we can now combine events within 100ms or 200ms to detect candidate events
%[MASK,MASKf] = Cand_Events(LLA,dLLA,dsig,mask,fs,d,T1,T2,Tpad)
[MASK,CEvents,UpSamEvents] = Cand_Events(LLR,dLLA,XX,dsig,mask,fs,d,0.025,0.25,0.5);
%%
l = 0;
% LabelEvents(dsig(:,(25000*l/d+1):25000*(l+2)/d),MASK(:,(25000*l/d+1):25000*(l+2)/d),round(fs/d),0.07);
% title('Amp.*Diff') % Only works with LLA for some reason
LabelEvents(dLFPs(:,(25000*l/d+1):25000*(l+2)/d),MASK(:,(25000*l/d+1):25000*(l+2)/d),round(fs/d),2.2);
title('NMF Detected Spikes') % Only works with LLA for some reason

%% Now clip the events
[DwnsmClippedEvents,DwnsmClippedEventsNorm] = clip_events(dLFP,CEvents);
[ClippedEvents,ClippedEventsNorm] = clip_events(LFPsmooth,UpSamEvents);
%%
FS = fs/3;
save('clipped_spikesV3.mat','FS','fs','ClippedEvents','ClippedEventsNorm','DwnsmClippedEvents','DwnsmClippedEventsNorm')

%%
figure
for i = 1:100
    CC = ClippedEvents{67}{i};
    if sum(size(CC)) > 2
        subplot(10,10,i)
        t = 1:length(CC); plot(t,CC); xlim([0,length(t)])
    end
end
sgtitle('Example Spikes Waveforms')
figure
for i = 1:100
    CC =DwnsmClippedEvents{67}{i};
    if sum(size(CC)) > 2
        subplot(10,10,i)
        t = 1:length(CC); plot(t,CC); xlim([0,length(t)])
    end
end
sgtitle('Example Spikes Waveforms Dwnsm')
%%
figure
for i = 1:1000
    CC = ClippedEvents{67}{i};
    if sum(size(CC)) > 2
        t = 1:length(CC); plot(t,CC);title('Spike Waveforms')
    end
    hold on
end
hold off
%%
figure
for i = 1:100
    CC = ClippedEventsNorm{67}{i};
    if sum(size(CC)) > 2
        subplot(10,10,i)
        t = 1:length(CC); plot(t,CC)
        xlim([0,length(t)]); xlim([0,length(t)])
        yline(1);yline(-1);xline(299);ylim([-1.5 1.5])
    end
end