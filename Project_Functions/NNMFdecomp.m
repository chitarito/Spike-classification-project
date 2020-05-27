function BFs = NNMFdecomp(Data,fs,Rank) 
% Algorith to decompose the data using NNMF
    %---------- Find Basis matrices for ranks 1 through 'Rank' ------------
    BW = num2cell(zeros(1,Rank)); BH = num2cell(zeros(1,Rank)); BD = num2cell(zeros(1,Rank));
    for i = 1:Rank
        % First, we obtained RMS of the reconstructions for 30 sets of randomly...
        ...initialized BFs and five iterations with multiplicative updates.  
        opt = statset('MaxIter',5); %,'Display','final');
        [W0,H0] = nnmf(Data,i,'replicates',30,'options',opt,'algorithm','mult');
        % This will give us output the basis with the smallest RMS to use as...
        ...initial basis for the full decomposition
  
        % Second, we obtained a full decomposition using alternating least squares...
        ...after initialization with the best set of bases found in the first step.
        OP = statset('Display','final');
        [W,H,D] = nnmf(Data,i,'algorithm','als','w0',W0,'h0',H0,'options',OP);
        BW{i} = W; BH{i} = H; BD{i} = D;
    end
    BFs = struct('BWs',BW,'BHs',BH,'BDs',BD); % Create a structure to hold the Basis functions
end