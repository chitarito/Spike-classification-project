function [Ws,Hs] = plotBFs(BFs,Rank,fs)
% Now plot these BWs and BHs to visualize
    Ws = figure(1)
    for j = 1:length(BFs)
        WV = num2cell(string(zeros(1,j)));
        for i = 1:j
            y = sprintf('W%d', i); WV{i} = y;
        end
        ran = Rank; i = 1;
        while i < ran
            if rem(sqrt(ran),1)==0
                break
            else 
                i = i+1;
                ran = ran-1;
            end
        end
        subplot(1,length(BFs),j)
        imagesc(BFs(j).BWs);ylabel('Channels'); title('Basis Vectors');
        xticks([1:j]); xticklabels(WV);colormap(flipud(gray(256)));
    end 

    Hs = figure(2)
    for j = 1:length(BFs)
        HV = num2cell(string(zeros(1,j)));
        for i = 1:j
            x = sprintf('H%d', i); HV{i} = x;
        end
        ran = Rank; i = 1;
        while i < ran
            if rem(sqrt(ran),1)==0
                break
            else 
                i = i+1;
                ran = ran-1;
            end
        end
        subplot(1,length(BFs),j)
        t = linspace(0,length(BFs(j).BHs),length(BFs(j).BHs));
        stackedPlot(t,HV,BFs(j).BHs);title('Activation Matrix')
        xlabl = floor(length(BFs(j).BHs)/fs);X=[0:xlabl];string(X);
        xticks([0:fs:length(BFs(j).BHs)]);xticklabels({X});xlabel('time (s)');
    end
end