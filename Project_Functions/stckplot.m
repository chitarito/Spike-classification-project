function PLOT = stckplot(Data,fs,time,gain) % Vector of data (i.e can give LL and EEG) 
% Make so it can take in interval?
    idx = time; data = Data;
    PLOT = figure
    if   1 <idx | length(Data) < idx
        print('Time outside of bounds')
    else
        % get the maximum range of the input data
        mx1 = max(data(:));mn1 = min(data(:));spacing1 = 0.1;

        [Ny1, Nx1] = size(data); y1 = 1:Ny1; x1 = 1:Nx1;
        rng1 = (1 + spacing1) * (mx1 - mn1);

        b1 = subplot(2,1,2)
        % create a series of linear plots of each row in data offset by the maximum data range
        plot(x1,data*gain+repmat(rng1*(Ny1:-1:1).', 1, Nx1), 'k-','Linewidth',.5);
        axis tight; ylims1 = get(gca, 'YLim'); ylim([(ylims1(1)-10) ylims1(2)+10]);
        yticks1 = (1:Ny1) * rng1 + (mx1 + mn1) / 2;
        try
            % use new name for YTickLabel
            set(gca, 'YTick', yticks1, 'YTickLabel', fliplr(y1), 'YLim', ylims1);

        catch
            % use old name for YTickLabel
            set(gca, 'YTick', yticks1, 'YTickLabels', fliplr(y1), 'YLim', ylims1); 
        end
        xlabl1 = floor(length(data)/fs);X1 = [0:xlabl1];string(X1);
        xticks([0:fs:length(data)]); xticklabels({X1}); xlabel('time (s)');
    end
end