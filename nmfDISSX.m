function diss = nmfDISSX(W1,W2)  
% USAGE: -----------------------------------------------------------------
% INPUTS
% W1 is the output of a NMF fit with K factors
% W2 is the output of a different NMF with the same K
% ------------------------------------------------------------------------
% OUTPUTS
% diss:     Diss is a measure of the dissimilarity of two factorizations.
% ------------------------------------------------------------------------
    [~,K] = size(W1);
    for i = 1:K
        for j = 1:K
        % We need to find the cross correlation between the columns of the
        % matrix W1 and the columns of the matrix W2
        A = W1(:,i); B = W2(:,j);
        C(i,j) = xcorr(A,B,0,'coeff');
        end
    end
    maxrow = max(C,[],1); maxcol = max(C,[],2); 
    maxrow(isnan(maxrow)) = 0; maxcol(isnan(maxcol)) = 0;
    diss = (1/(2*K))*(2*K - sum(maxrow) -sum(maxcol));
end