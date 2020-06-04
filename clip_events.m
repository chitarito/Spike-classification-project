function [ClippedEvents,ClippedEventsNormalized] = clip_events(dLFP,CEvents)
    ClippedE = CEvents;
    ClippedEvents = {};
    ClippedEventsNormalized = {};
    % Now repeat procedure
    for i = 1:size(dLFP,1) % for each channel
        C = CEvents{i}; % obtain the events of this channel
        for j = 1:length(C) % for each event in the channel
            if sum(size(C)) > 2 % if it doesn't have events, make zero
                CC = C{j};
                if sum(size(CC)) > 2 % if the event is not empty
                    ClippedE{i}{j} = dLFP(i,C{j});
                end
            end
        end
        
        CCC = ClippedE{i};
        ZX = zeros(length(CCC),2);
        for jj = 1:length(CCC)
            if sum(size(CCC)) > 2
                ZX(jj,:) = size(CCC{jj});
            end
        end
        zx = max(ZX); % This gives the length of an event to create matrix
        
        % Now go over again and if the events are the exact same only keep one
        b = zeros(length(CCC),max(zx)); % create a matrix for each event
        for ii = 1:length(CCC) % Now go through each event that is non empty
            if sum(size(CCC)) > 2
                if sum(size(CCC{ii})) > 2
                    b(ii,:) = CCC{ii}; % find empty events
                end
            end
        end
        % Now find the unique rows (unique events because we may have multiple of the same ones)
        A = b; [C,ia,ic] = unique(A,'rows') ;
        iwant = zeros(size(A)); iwant(ia,:) = C ;
        Non_copyE = iwant(any(iwant,2),:);
        
        % Now each row is an event so clip it
        ClippedEvents{i} = cell(1,size(Non_copyE,1));
        ClippedEventsNormalized{i} = cell(1,size(Non_copyE,1));
        for k = 1:size(Non_copyE,1) % for each event in the channel
            Z = median([1:length(Non_copyE(k,:))]);
            ClippedEvents{i}{k} = Non_copyE(k,:);
            ClippedEventsNormalized{i}{k} = Non_copyE(k,:)/abs(Non_copyE(k,Z));
        end
    end
end