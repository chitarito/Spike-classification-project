function [W,H,D] = NMFdecomp(data,K)

%takes a couple minutes to run on an hour of data sampled at 198Hz data
Data = data;
% In this case K = 4 seems to be the best option
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
end