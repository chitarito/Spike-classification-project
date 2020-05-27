function Visfigs = View_figs(Data,fs,On)
    Visfigs = figure;
    
    for i = 1:length(Data)
        data = Data{i};
        gain = 1;
        % get the maximum range of the input data
        mx1 = max(data(:));mn1 = min(data(:));spacing1 = 0.1;
        [Ny1, Nx1] = size(data); y1 = 1:Ny1; x1 = 1:Nx1;
        rng1 = (1 + spacing1) * (mx1 - mn1);
        
        % create a series of linear plots of each row in data offset by the maximum data range
        subplot(size(Data,1),1,i)
        plot(x1,data*gain+(repmat(rng1*(Ny1:-1:1).', 1, Nx1)),'b-','Linewidth',1.2);
%         hold on
%         dat = zeros(size(data)); dat([8 10 12:15],:) = data([8 10 12:15],:);
%         Mask = double(logical(dat)); Mask(Mask == 0) = NaN;
%         plot(x1,(data*gain+(repmat(rng1*(Ny1:-1:1).', 1, Nx1))).*Mask,'k-','Linewidth',2);
%         
        axis tight; ylims1 = get(gca, 'YLim'); ylim([0 ylims1(2)+10]);
        yticks1 = (1:Ny1) * rng1 + (mx1 + mn1) / 2;
        try
            % use new name for YTickLabel
            set(gca, 'YTick', yticks1, 'YTickLabel', fliplr(y1), 'YLim', ylims1);
            
        catch
            % use old name for YTickLabel
            set(gca, 'YTick', yticks1, 'YTickLabels', fliplr(y1), 'YLim', ylims1);
        end
%         xlabl1 = floor(length(data)*9/fs); X1 = [0:10:xlabl1];string(X1);
%         xticks([0:(10*fs):(length(data))]); xticklabels({X1});
        %set(gca,'xtick',[])
        set(gca,'xticklabel',[])
        xlabel('time (s)');
        %ylabel('Activation Vectors')
        ylabel('Spike Rate')
       
        for j = 1:length(On);
            xline(7934,'r','Linewidth',2);
            xline((7934+11910),'b','Linewidth',2);
        end
    end
end