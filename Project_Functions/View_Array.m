function DataVis = View_Array(Data,fs,gain1)
    % View imagesc of data
    %----------------------Identify Bais Set of interest-----------------------
    % Plot the Data, H and W matrices
    data = Data; % Line Length only active electrodes
    
    DataVis = figure
    
    % get the maximum range of the input data
    mx1 = max(data(:));mn1 = min(data(:));spacing1 = 0.1;
    [Ny1, Nx1] = size(data); y1 = 1:Ny1; x1 = 1:Nx1;
    rng1 = (1 + spacing1) * (mx1 - mn1);

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
    xticks([0:fs:length(data)]); xticklabels({X1}); xlabel('time (s)')
    grid on
end