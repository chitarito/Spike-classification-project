function [MASKf,Events,USEvents] = Cand_Events(LLR,dLLA,XX,dsig,mask,fs,d,T1,T2,Tpad)
% it's time to clip them for each channel
    % first we must go through each channel and merge events within 200ms of each other
    % if our sampling rate is FS = fs/3, 200ms = 0.2s*FS = 39.697 samples
    % thus if two events are within 200ms, they will be combined to a single event
    FS = fs/d;
    % to do this we will use the mask for labeling events
    MasK = mask;
    % for each channel we need to find the time points of events

    % For each channel, create an object to hold the events for each channel
    CHEvents = cell(1,size(MasK,1));

    for j = 1:size(MasK,1)
    % let's start with a single channel
    
    % If the zeros are less than 200ms apart merge ones
    % so I want to find the zeros, then I want to find the consecutive
    % zeros which will have a difference of 1. Then if that string of ones
    % is less than .2*FS samples, then I'll set those zeros equal to one.
        a = MasK(j,:);
        aa = find(~a);
        
    % Find time points of events
        k = [true;diff(aa(:))~=1]; % find the difference in IsZ where its not a difference of 1...
    ...indicating a new event. Start it with a true because the first part is also and event.
        s = cumsum(k); % counts the number of events
        x =  histc(s,1:s(end)); %counts the number of indices in each event
        idx = find(k); % Find the index of the beginning of each event in IsZ
        N = length(x);
        AAA = cell(1,N);
        for i = 1:N
            index = find(s == i); % find the indices for each event in s
            AAA{i} = aa(index); % find actual indices in data from IsZ...
        ...now these indices are the indices for each event. Put them in cell array
        end

    % boom, we've got the indices for the events of string of zeros
    % Now combine them together if under T1ms (i.e. 200ms)
        for p = 1:length(AAA)
            if length(AAA{p}) < round(T1*FS)
                a(AAA{p}) = 1;
            end
        end
        AC = a;

        IsZ = find(AC);
        if isempty(IsZ)
            CHEvents{j} = 0; 
            display(sprintf('Channel %d has no spikes',j))
        else
        % Find time points of events
            k = [true;diff(IsZ(:))~=1]; % find the difference in IsZ where its not a difference of 1...
        ...indicating a new event. Start it with a true because the first part is also and event.
            s = cumsum(k); % counts the number of events
            x =  histc(s,1:s(end)); %counts the number of indices in each event
            idx = find(k); % Find the index of the beginning of each event in IsZ
            NumberofEvents = length(x);
            IVEvents = cell(1,length(x));
            for ii = 1:length(x)
                index = find(s == ii); % find the indices for each event in s
                IVEvents{ii} = IsZ(index); % find actual indices in data from IsZ...
            ...now these indices are the indices for each event. Put them in cell array
            end
            IVEN = IVEvents;

            % now go through the events in the channel and if it is longer than
            % 250ms, remove event
            for kk = 1:length(IVEN)
                Yy = IVEN{kk};
            % If event is longer than 250ms, remove
                if length(Yy) > round(T2*FS)
                    IVEN{kk} = 0;
                end
            end
            
            Bb = zeros(1,length(IVEN));
            for kK = 1:length(IVEN)
                if sum(size(IVEN{kK})) > 2
                    Bb(kK) = 1;
                end
            end
            % find how many empty elements there are
            BB = find(~Bb);
            IVN = cell(1,(length(IVEN)-length(BB))); % subtract length of empty elements

            idn = zeros(size(IVN));
            for q = 1:length(IVEN)
                Cc = IVEN{q};
                if sum(size(Cc)) > 2
                    idn(q) = 1; % now this is the vector of new events for the channel
                end
            end
            IDN = find(idn); % find the non zero elements of idn
            for r = 1:length(IDN)
                IVN{r} = IVEN{IDN(r)};
            end
            CHEvents{j} = IVN;
        end
    end
%
    % Let's see what this mask looks like: Generate new mask
    mASk = zeros(size(MasK));
    for jj = 1:size(MasK,1)
        A = mASk(jj,:);
        % Now go through each event for each channel
        C = CHEvents{jj};
        if sum(size(C)) > 2
            for i = 1:length(CHEvents{jj})
                B = [C{i}];
                A(B) = 1;
            end
            mASk(jj,:) = A;
        else
            mASk(jj,:) = mASk(jj,:);
        end
    end
%
    % Now I have to find the peak of each of these events and extend them 
    ...then find a peak again and extend to generate clip
    % Then, for each channel, and for each event, find the maximum amplitude (or first?)
    ... then add 100ms before and 200ms after
    CHE = CHEvents;
    for ij = 1:size(dLLA,1) % for each channel
        C = CHEvents{ij};
        for jn = 1:length(C)
            if sum(size(C)) > 2
                % find the peak of the amp.*diff in the event time points
                Zt = dsig(ij,C{jn});
                [~,Index] = max(Zt);
                % Then create a new event containing cell and change the ...
                ...events time points by selecting the 'peak' time point...
                ... and then add 200ms before it and 200ms after
                L = (C{jn}(Index)-round(FS*0.2)); L(L < 0) = 0; 
                LN = L(L ~= 0); LLN = min(LN);
                
                F = (C{jn}(Index)+round(FS*0.2)); F(F > length(dLLA)) = 0; 
                FN = F(F ~= 0); FFN = max(FN);
                
                CHE{ij}{jn} = [LLN:FFN];
            end
        end
    end
 %   
    % Now create mask and see
    MASK = zeros(size(mASk));
    for je = 1:size(mASk,1)
        A = MASK(je,:);
        % Now go through each event for each channel
        C = CHE{je};
        if sum(size(C)) > 2
            for iii = 1:length(C)
                B = [C{iii}];
                A(B) = 1;
            end
            MASK(je,:) = A;
        else
            MASK(je,:) = MASK(je,:);
        end
    end
%
    % Now repeat procedure to make sure we are on the peak
    Events = CHE;
    for ie = 1:size(dLLA,1) % for each channel
        C = CHE{ie};
        for jq = 1:length(C)
            if sum(size(C)) > 2
                % find the peak of the amp.*diff in the event time points
                Zs = dsig(ie,C{jq});
                [~,J] = max(Zs);
                % Then create a new event containing cell and change the ...
                ...events time points by selecting the 'peak' time point...
                    ... and then add 500ms before it and 500ms after
                L = (C{jq}(J)-round(FS*Tpad)); L(L < 0) = 0; 
                LN = L(L ~= 0); LLN = min(LN);
                
                F = (C{jq}(J)+round(FS*Tpad));
                F(F > length(dLLA)) = 0; FN = F(F ~= 0); FFN = max(FN);
               
                Events{ie}{jq} = [LLN:FFN];
            end
        end
    end
    % Now create mask and see
    MASKf = zeros(size(MASK));
    for jr = 1:size(MASK,1)
        A = MASKf(jr,:);
        % Now go through each event for each channel
        C = Events{jr};
        if sum(size(C)) > 2
            for i = 1:length(C)
                B = [C{i}];
                A(B) = 1;
            end
            MASKf(jr,:) = A;
        else
            MASKf(jr,:) = MASKf(jr,:);
        end
    end
    
    % Now create another events cell containing the regular sampling indices
    UpSamEvents = CHE;
    for ir = 1:length(CHE)
        if sum(size(CHE{ir})) > 2
            for js = 1:length(CHE{ir})
                if sum(size(CHE{ir}{js})) > 2
                    C = CHE{ir}{js};
                    Zz = median([1:length(C)]);
                    % now take the min and max indices and multiply by d
                    Lq = (C(Zz)*d-round(fs*Tpad)); Lq(Lq < 0) = 0;
                    LNq = Lq(Lq ~= 0); LLNq = min(LNq);
                    
                    Fq = (C(Zz)*d+round(fs*Tpad)); Fq(Fq > length(LLR)) = 0;
                    FNq = Fq(Fq ~= 0); FFNq = max(FNq);
                    
                    UpSamEvents{ir}{js} = [LLNq:FFNq];
                end
            end
        end
    end
    
    USEvents = UpSamEvents;
    for il = 1:size(dLLA,1) % for each channel
        Cq = UpSamEvents{il};
        for jl = 1:length(Cq)
            if sum(size(Cq)) > 2
                Z = XX(il,Cq{jl});
                [~,Jg] = max(Z);
                % Then create a new event containing cell and change the ...
                ...events time points by selecting the 'peak' time point...
                    ... and then add 100ms before it and 200ms after
                Ln = (Cq{jl}(Jg)-round(fs*Tpad)); Ln(Ln < 0) = 0; 
                LNn = Ln(Ln ~= 0); LLNn = min(LNn);
                
                Fn = (Cq{jl}(Jg)+round(fs*Tpad));
                Fn(Fn > length(LLR)) = 0; FNn = Fn(Fn ~= 0); FFNn = max(FNn);
               
                USEvents{il}{jl} = [LLNn:FFNn];
            end
        end
    end
end