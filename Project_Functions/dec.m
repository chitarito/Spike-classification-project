function deLL = dec(Data,d)
% Give data and decimation fator d.
    LLAll = Data;
    DeLLA={};
    for i = 1:size(LLAll,1)
        DeLLA{i} = decimate(LLAll(i,:),d,'fir');
    end
    deLL = cell2mat(DeLLA'); % Uses a finite impulse response filter...
...designed using the window method with a hamming window. The ...
    ... filter has an order of 30
end