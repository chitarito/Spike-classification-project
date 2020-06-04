function plotW = plotW(W)
    j = size(W,2);
    plotW = figure
    imagesc(flipud(W));ylabel('Channels'); title('Basis Vectors');
    xticks([1:j]); WV = num2cell(string(zeros(1,j)));
    for i = 1:j
        y = sprintf('W%d', i); WV{i} = y;
    end
    xticklabels(WV);colormap(flipud(gray(256)));
    yticks([1:size(W,1)]);yticklabels(num2cell(fliplr(string(1:size(W,1)))))
end