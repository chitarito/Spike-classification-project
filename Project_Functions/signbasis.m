function SigBs = signbasis(BWCorr,I)
    zzz = BWCorr(:,:,I);
    
    % "normalize" the columns
    for i = 1:size(zzz,2)
        zzz(:,i) = zzz(:,i)./norm(zzz(:,i));
    end
    % Now find max of each column 1st
    ZZZ = max(zzz,[],1);
    for i = 1:size(zzz,2)
        A = zzz(:,i); A(A < ZZZ(i)) = 0; zzz(:,i) = A;
    end
    % Now find max of each row 2nd as it may have correlated across, but
    % only one truly correlates?
    ZZZ2 = max(zzz,[],2);
    for i = 1:size(zzz,2)
        A = zzz(i,:); A(A < ZZZ2(i)) = 0; zzz(i,:) = A;
    end
    SigBs = zzz;
    % Now the spiking basis are the rows with non zero elements and their
    % column is which basis set they belong to. Can assign a color to each
    % column.
    % label significant basis in set
    str = string(find(sum(SigBs,2) ~=0)'); 
    display(['Significant Basis:' str]);
end