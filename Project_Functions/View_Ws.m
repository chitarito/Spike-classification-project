function VisWs = View_Ws(Ws,fs)
    VisWs = figure;
    
    W = Ws;
    Zy = sum(Ws,[1 3]);
    Zyy = Zy(:,any(Zy,1));
    for i = 1:size(W,2)
        vector = squeeze(W(:,i,:));
        if (sum(sum(vector)) == 0)
            Data = squeeze(W(:,i,:));
            display(sprintf('Basis %d is zero',i))
        else
            data = squeeze(W(:,i,:));
            gain = 2;
            % get the maximum range of the input data
            mx1 = max(data(:));mn1 = min(data(:));spacing1 = 0.1;
            [Ny1, Nx1] = size(data); y1 = 1:Ny1; x1 = 1:Nx1;
            rng1 = (1 + spacing1) * (mx1 - mn1);
            
            % create a series of linear plots of each row in data offset by the maximum data range
            subplot(1,length(Zyy),i)
            plot(x1,data*gain+(repmat(rng1*(Ny1:-1:1).', 1, Nx1)),'b-','Linewidth',.5);
            
            axis tight; ylims1 = get(gca, 'YLim'); ylim([0 ylims1(2)+10]);
            yticks1 = (1:Ny1) * rng1 + (mx1 + mn1) / 2;
            try
                % use new name for YTickLabel
                set(gca, 'YTick', yticks1, 'YTickLabel', fliplr(y1), 'YLim', ylims1);
                
            catch
                % use old name for YTickLabel
                set(gca, 'YTick', yticks1, 'YTickLabels', fliplr(y1), 'YLim', ylims1);
            end
            xlabl1 = floor(1000*size(data,2)/fs); X1 = [0:50:xlabl1];string(X1);
            xticks([0:(size(data,2)/(xlabl1/50)):size(data,2)]); xticklabels({X1}); xlabel('time (ms)');
            title(sprintf('W %d',i))
        end
    end
end