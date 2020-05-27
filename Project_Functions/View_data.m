function DataVis = View_data(Data,fs,BFs,Desired_Rank,gain1)
    % View imagesc of data
    %----------------------Identify Bais Set of interest-----------------------
    % Plot the Data, H and W matrices
    data = Data; % Line Length only active electrodes
    k = Desired_Rank; % Plot Basis for this Rank
    DataVis = figure
    
    % get the maximum range of the input data
    mx1 = max(data(:));mn1 = min(data(:));spacing1 = 0.1;
    [Ny1, Nx1] = size(data); y1 = 1:Ny1; x1 = 1:Nx1;
    rng1 = (1 + spacing1) * (mx1 - mn1);

    a1 = subplot(2,2,1)
    % create a series of linear plots of each row in data offset by the maximum data range
    plot(x1,data*gain1+repmat(rng1*(Ny1:-1:1).', 1, Nx1), 'k-','Linewidth',.5);
    axis tight; ylims1 = get(gca, 'YLim'); ylim([(ylims1(1)-10) ylims1(2)+10]);
    yticks1 = (1:Ny1) * rng1 + (mx1 + mn1) / 2; title('ECoG')
    try
        % use new name for YTickLabel
        set(gca, 'YTick', yticks1, 'YTickLabel', fliplr(y1), 'YLim', ylims1);

    catch
        % use old name for YTickLabel
        set(gca, 'YTick', yticks1, 'YTickLabels', fliplr(y1), 'YLim', ylims1); 
    end
    xlabl1 = floor(length(data)/fs);X1 = [0:xlabl1];string(X1);
    xticks([0:fs:length(data)]); xticklabels({X1}); xlabel('time (s)');

    
    WV = num2cell(string(zeros(1,k)));
    for i=1:k
        y = sprintf('W%d', i); WV{i} = y;
    end
    a2 = subplot(2,2,2)
    imagesc(BFs(k).BWs); title('Basis Vectors'); ylabel('Channels');
    xticks([1:k]); xticklabels(WV); colormap(flipud(gray(256)));axis tight;

    HV = num2cell(string(zeros(1,k)));
    for i=1:k
        x = sprintf('H%d', i); HV{i} = x;
    end
    a3 = subplot(2,2,3)
    H = BFs(k).BHs(:,1:length(Data)); t = linspace(0,length(H),length(H)); stackedPlot(t,HV,H); title('H');
    xlabl = floor(length(H)/fs); X=[0:xlabl]; string(X);
    xticks([0:fs:length(H)]); xticklabels({X}); xlabel('time (s)');

    linkaxes([a1,a3], 'x');
end