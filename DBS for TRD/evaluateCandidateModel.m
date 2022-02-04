function [modelNRMSE, DMax, TOuter, kSelected] = evaluateCandidateModel(Y, blockLengths, scores)

    iBlock = [0 cumsum(blockLengths)];
    blockAvgs = zeros(length(blockLengths), size(Y, 2));
    for b = 1:length(blockLengths)
        blockAvgs(b,:) = mean(Y((iBlock(b)+1):iBlock(b+1),:));
    end

    %%% Ridge Regression with inner level of cross-validation to determine optimal k
    K = logspace(-4, 5, 100);
    predScoresLOK = zeros(length(blockAvgs), length(K));
    for b = 1:length(blockAvgs)
        blockAvgsLOO = blockAvgs; blockAvgsLOO(b,:) = [];
        scoresLOO = scores; scoresLOO(b) = [];
        TK = ridge(scoresLOO, blockAvgsLOO, K, 0); 
        predScoresLOK(b, :) = blockAvgs(b,:) * TK(2:end, :) + TK(1,:);
    end
    NRMSEK = ModifiedFunction3_performanceMeasurement(scores, predScoresLOK);
    [~, iKMin] = min(NRMSEK);
    kSelected = K(iKMin);
    
    %Fit a model to all training IMS
    TOuter = ridge(scores, blockAvgs, kSelected, 0);
    predScoresOuter = blockAvgs * TOuter(2:end) + TOuter(1);
    
    %%% Evaluate candidate model prediction error and sensitivity on training data with inner level CV
    NRMSE = zeros(1,length(scores)); D = zeros(1,length(scores));
    dMat = [ones(1,size(blockAvgs, 2)); blockAvgs]; d = svd(dMat);
    f = sum((d.^2)./(d.^2 + kSelected));
    sig2 = sum((predScoresOuter - scores).^2) / (length(scores) - f);
    
    predScoresLO = zeros(size(scores));
    for j = 1:length(scores)
        scoresLOO = scores; scoresLOO(j) = []; %s0LOO = mean(scoresLOO);
        blockAvgsLOO = blockAvgs; blockAvgsLOO(j,:) = [];
        TInner = ridge(scoresLOO, blockAvgsLOO, kSelected, 0);
        predScoresLO(j) = blockAvgs(j,:) * TInner(2:end) + TInner(1); 
        
        %sensitivity metric
        predScoresInner = blockAvgs * TInner(2:end) + TInner(1);
        D(j) = sum((predScoresOuter - predScoresInner).^2)/(f * sig2);        
    end
    DMax = max(D);
    modelNRMSE = ModifiedFunction3_performanceMeasurement(scores, predScoresLO);
end