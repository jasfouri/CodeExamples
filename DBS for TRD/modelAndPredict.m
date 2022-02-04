function [testNRMSE, selectedModel, coeffMag, PAll, TAll] = modelAndPredict(scores, features, PCOptions, channels, regions, iSelectedRegions)
    disp("Selecting candidate models")
    allCandidateNRMSE = NaN(length(PCOptions), length(regions), length(scores)); allCandidateD = allCandidateNRMSE;
    selectedModel = zeros(length(scores), 2); %first col is number of PCs, second col is region
    
    %for CV, consider creating matrix where each row is scores for each stage of cross validation and operating on this matrix - may be faster
    for s = 1:length(scores) %outer level of CV (LOO = leave-one-out)
        scoresLOO = scores; scoresLOO(s) = [];
        featuresLOO = features; featuresLOO(s) = []; featuresLOOConcat = vertcat(featuresLOO{:});
        blockLengths = cellfun('size', featuresLOO, 1);
        iNewRegions = 1:length(regions); 
        if ~isempty(iSelectedRegions); iNewRegions(ismember(iNewRegions, iSelectedRegions(s,:))) = []; end %only executes past first stage of progressive region selection
        for r = iNewRegions %ofc, acc, dpfc, vpfc
            
            %for smaller regions 002:
%             if (r == 7 || r == 10); PCOptions = [1 5 10]; else; PCOptions = [1 5 10 15]; end
            
            iRegion = find(ismember(channels, regions{r}));
            if ~isempty(iSelectedRegions); iRegion = union(iRegion, find(ismember(channels, horzcat(regions{ismember(1:length(regions), iSelectedRegions(s,:))})))); end

            %%% PCA
            Z = reshape(featuresLOOConcat(:, iRegion, :), size(featuresLOOConcat, 1), length(iRegion)*size(featuresLOOConcat, 3)); %all spectral features from all channels combined into a 2D matrix
            Q = cov(Z); [V, D] = eig(Q); %first col of V is eigenvector corresponding to smallest eigenvalue
            eigenvecs = flip(V, 2)'; %rows of eigenvecs are PCs in descending order of eigenvalues

            for p = 1:length(PCOptions)
                P = eigenvecs(1:PCOptions(p), :);
                Y = (P * Z')'; 
                [allCandidateNRMSE(p,r,s), allCandidateD(p,r,s), ~, ~] = evaluateCandidateModel(Y, blockLengths, scoresLOO);
            end
        end
        
        %select model based on sensitivity and prediction error metrics
        candidateNRMSE = allCandidateNRMSE(:,:,s);
        if ~any(any(allCandidateD(:,:,s) < 1))
            [selectedModel(s,1), selectedModel(s,2)] = find(candidateNRMSE == min(candidateNRMSE, [], 'all') & candidateNRMSE ~= 0); %first element of selectedModel is number of PCs, second element is region
        end
        DMax = 0.5; 
        while selectedModel(s,1) == 0
            minNRMSE = min(candidateNRMSE(allCandidateD(:,:,s) < DMax));
            if ~isempty(minNRMSE)
                [selectedModel(s,1), selectedModel(s,2)] = find(candidateNRMSE == minNRMSE);
            end
            DMax = DMax + 0.1;
        end       
        
%         disp(num2str(s) + "/" + num2str(length(scores)))
    end

    %% Use candidate model (one per fold) to predict test CAT-DI for each fold of CV
    disp("Predicting test CAT-DI with selected candidate models")
    testPred = zeros(size(scores));
    kSelected = zeros(size(scores));
    coeffMag = {};
    PAll = {}; TAll = {};

    for s = 1:length(scores) 
        iNewRegions = 1:length(regions);
        if ~isempty(iSelectedRegions); iNewRegions(ismember(iNewRegions, iSelectedRegions(s,:))) = []; end %only executes past first stage of progressive region selection
        p = selectedModel(s,1); r = selectedModel(s,2);
        
        %for smaller regions 002:
%         if (r == 7 || r == 10); PCOptions = [1 5 10]; else; PCOptions = [1 5 10 15]; end
        
        %first recalculate parameters for the selected model:
        iRegion = find(ismember(channels, regions{r}));
        if ~isempty(iSelectedRegions); iRegion = union(iRegion, find(ismember(channels, horzcat(regions{ismember(1:length(regions), iSelectedRegions(s,:))})))); end
        scoresLOO = scores; scoresLOO(s) = [];
        featuresLOO = features; featuresLOO(s) = []; featuresLOOConcat = vertcat(featuresLOO{:});
        blockLengths = cellfun('size', featuresLOO, 1);
        
        ZLOO = reshape(featuresLOOConcat(:, iRegion, :), size(featuresLOOConcat, 1), length(iRegion)*size(featuresLOOConcat, 3));
        Q = cov(ZLOO); [V, D] = eig(Q);
        eigenvecs = flip(V, 2)';
        P = eigenvecs(1:PCOptions(p), :); PAll{s} = P;
        Y = (P * ZLOO')';
        
        [~, ~, T, kSelected(s)] = evaluateCandidateModel(Y, blockLengths, scoresLOO); TAll{s} = T;

        %now calculate test CAT-DI point
        featuresTest = features{s};
        ZTest = reshape(featuresTest(:, iRegion, :), size(featuresTest, 1), length(iRegion)*size(featuresTest, 3));
        YTest = (P * ZTest')';
        
        testBlockAvg = mean(YTest);
        testPred(s) = testBlockAvg * T(2:end) + T(1);
        
        coeffMag{s} = T(2:end)' * P;
    end
    
    testNRMSE = ModifiedFunction3_performanceMeasurement(scores, testPred);
    
    plot(scores, testPred, 'ob'); hold on; 
    minBound = min([scores; testPred]) - 2; maxBound = max([scores; testPred]) + 2;
    plot((minBound):(maxBound), (minBound):(maxBound));
    xlim([minBound maxBound]); ylim([minBound, maxBound]);
    xlabel("True Score"); ylabel("Predicted Score");
    if ~isempty(iSelectedRegions)
        title("NRMSE = " + num2str(testNRMSE) + " | 2 Regions");
    else 
        title("NRMSE = " + num2str(testNRMSE) + " | 1 Region"); 
    end
    hold off; 
    
end