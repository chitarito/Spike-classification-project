function plotH = plotH(H,fs,d)
    j = size(H,1);
    plotH = figure
    HV = num2cell(string(zeros(1,j)));
    for i = 1:j
        x = sprintf('H%d', i); HV{i} = x;
    end
    mydata = H;
    t = linspace(0,length(mydata),length(mydata));
    stackedPlot(t,HV,mydata);title('Activation Matrix')
    xlabl = floor(length(mydata)*d/fs);X=[0:xlabl];string(X);
    xticks([0:(fs/d):length(mydata)]);xticklabels({X});xlabel('time (s)');
end
