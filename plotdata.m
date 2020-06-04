function Plot = plotdata(Data,fs,gain)
    mydata = [1:size(Data,1)]' + Data*gain; t = 1:length(mydata);
    Plot = figure
    plot(t,mydata,'b'); title('dLLA')
    xlabl1 = floor(length(mydata)*fs);X1 = [0:xlabl1];string(X1);
    xticks([0:fs:length(mydata)]); xticklabels({X1}); xlabel('time (s)')
    grid on
end