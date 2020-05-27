function DataV = PlotData(Data,fs,t0,tf,Ws,Hs,gain1,gain2,d) % gain1 for data, gain2 for H, d is downsampling factor
% give time you want to see in minutes start is 0

    %----------------------Identify Bais Set of interest-----------------------
    % Plot the Data, H and W matrices
    DataV = figure;
    % Plot from time I to time J. Time is in minutes
    I = t0;
    J = tf;
    data = Data(:,(1+I*60*fs):(J*60*fs)); % Line Length only active electrodes
    gain1 = 3.5;
    gain2 = 1.3;
    % get the maximum range of the input data
    mx1 = max(data(:));mn1 = min(data(:));spacing1 = 0.1;
    [Ny1, Nx1] = size(data); y1 = 1:Ny1; x1 = 1:Nx1;
    rng1 = (1 + spacing1) * (mx1 - mn1);
    
    a1 = subplot('Position',[0.5 0.4 0.45 0.55])
    % create a series of linear plots of each row in data offset by the maximum data range
    plot(x1,data*gain1+repmat(rng1*(Ny1:-1:1).', 1, Nx1), 'k-','Linewidth',.5);
    axis tight; ylims1 = get(gca, 'YLim'); ylim([(ylims1(1)-10) ylims1(2)+10]);
    yticks1 = (1:Ny1) * rng1 + (mx1 + mn1) / 2; title('LFP');
    try
        % use new name for YTickLabel
        set(gca, 'YTick', yticks1, 'YTickLabel', fliplr(y1), 'YLim', ylims1);
        
    catch
        % use old name for YTickLabel
        set(gca, 'YTick', yticks1, 'YTickLabels', fliplr(y1), 'YLim', ylims1);
    end
    xlabl1 = floor(length(data)/fs);X1 = [(0+I*60):(xlabl1+I*60)];string(X1);
    xticks([0:fs:length(data)]); xticklabels({X1}); xlabel('time (s)');
    set(gca, 'YGrid', 'off', 'XGrid', 'on')
    %
    wx = 0.02; wy = 0.011;
    ws = {};
    
    % first select iteration
    for i = 1:size(Ws,2)
        ws{i} = squeeze(Ws(:,i,:));
    end
    % Now ws is a cell with each element being the k basis Ws of each iteration
    for i = 1:length(ws)
        vector = squeeze(Ws(:,i,:));
        if (sum(vector~=0) == 0)
            data = squeeze(Ws(:,i,:));
            display(sprintf('Basis %d is zero',i))
            
        else
            data = ws{i};
            gain = 2;
            % get the maximum range of the input data
            mx1 = max(data(:));mn1 = min(data(:));spacing1 = 0.1;
            [Ny1, Nx1] = size(data); y1 = 1:Ny1; x1 = 1:Nx1;
            rng1 = (1 + spacing1) * (mx1 - mn1);
            
            % create a series of linear plots of each row in data offset by the maximum data range
            subplot('Position',[(i*wy+wx*(i-1)), 0.4, wx, 0.55])
            plot(x1,data*gain+repmat(rng1*(Ny1:-1:1).', 1, Nx1),'b-','Linewidth',.5);
            
            axis tight; ylims1 = get(gca, 'YLim'); ylim([0 ylims1(2)+10]);
            yticks1 = (1:Ny1) * rng1 + (mx1 + mn1) / 2;
            try
                % use new name for YTickLabel
                set(gca, 'YTick', yticks1, 'YTickLabel', fliplr(y1), 'YLim', ylims1);
                
            catch
                % use old name for YTickLabel
                set(gca, 'YTick', yticks1, 'YTickLabels', fliplr(y1), 'YLim', ylims1);
            end
            title(sprintf('W%d', i));
            
            set(gca, 'YGrid', 'off', 'XGrid', 'on')
        end
    end
    
    % k = 5;
    % %
    % HV = num2cell(string(zeros(1,k)));
    % for i=1:k
    %     x = sprintf('H%d', i); HV{i} = x;
    % end
    
    % plot the Hs
    Hsr={};
    for i = 1:size(Hs,1)
        Hsr{i} = resample(Hs(i,:),d,1);
    end
    
    H = cell2mat(Hsr'); % Uses a finite impulse response filter...
    ...designed using the window method with a hamming window. The ...
        ... filter has an order of 30
    H = H(:,1:length(Data));
    H = H(any(H,2),:);
    Data = H(:,(1+I*60*fs):(J*60*fs)); % Line Length only active electrodes
    %gain2 = 1.5;
    
    % get the maximum range of the input data
    mx1 = max(Data(:));mn1 = min(Data(:));spacing1 = 0.1;
    [Ny1, Nx1] = size(Data); y1 = 1:Ny1; x1 = 1:Nx1;
    rng1 = (1 + spacing1) * (mx1 - mn1);
    
    a2 = subplot('Position',[0.5 0.05 0.45 0.29])
    plot(x1,Data*gain2+repmat(rng1*(Ny1:-1:1).', 1, Nx1), 'k-','Linewidth',.5);
    axis tight; ylims1 = get(gca, 'YLim'); ylim([(ylims1(1)-10) ylims1(2)+10]);
    yticks1 = (1:Ny1) * rng1 + (mx1 + mn1) / 2;
    try
        % use new name for YTickLabel
        set(gca, 'YTick', yticks1, 'YTickLabel', fliplr(y1), 'YLim', ylims1);
        
    catch
        % use old name for YTickLabel
        set(gca, 'YTick', yticks1, 'YTickLabels', fliplr(y1), 'YLim', ylims1);
    end
    title('H');
    xlabl1 = floor(length(Data)/fs); X1 = [(0+I*60):(xlabl1+I*60)];string(X1);
    xticks([0:fs:length(Data)]); xticklabels({X1}); xlabel('time (s)');
    set(gca, 'YGrid', 'off', 'XGrid', 'on')
    linkaxes([a1,a2], 'x');
end