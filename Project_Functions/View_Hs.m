% Plot Hs
function VisHs = View_Hs(Hs,fs,gain)
    if sum(Hs,'all') == 0;
        
    else
        Hs = Hs(any(Hs,2),:);
        Data = Hs; % Line Length only active electrodes
        gain2 = gain;
        
        % get the maximum range of the input data
        mx1 = max(Data(:));mn1 = min(Data(:));spacing1 = 0.1;
        [Ny1, Nx1] = size(Data); y1 = 1:Ny1; x1 = 1:Nx1;
        rng1 = (1 + spacing1) * (mx1 - mn1);
        
        figure
        plot(x1,Data*gain2+repmat(rng1*(Ny1:-1:1).', 1, Nx1), 'k-','Linewidth',1);
        axis tight; ylims1 = get(gca, 'YLim'); ylim([(ylims1(1)-10) ylims1(2)+10]);
        yticks1 = (1:Ny1) * rng1 + (mx1 + mn1) / 2; title('H')
        try
            % use new name for YTickLabel
            set(gca, 'YTick', yticks1, 'YTickLabel', fliplr(y1), 'YLim', ylims1);
            
        catch
            % use old name for YTickLabel
            set(gca, 'YTick', yticks1, 'YTickLabels', fliplr(y1), 'YLim', ylims1);
        end
        xlabl1 = floor(length(Data)/fs); X1 = [1:length(Hs)];string(X1);
        xticks([0:fs:length(Data)]); xticklabels({X1}); xlabel('time (s)');
end