% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Omid Sani, Yuxiao Yang, Maryam Shanechi
%   Shanechi Lab, University of Southern California, 2018
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NRMSE = ModifiedFunction3_performanceMeasurement(trueVals, predictedValsMat)
    
    % Computes the evaluation metric (NRMSE) for the decoding performance
    % measurement
    % Inputs:
    % - (1) trueVals: vector of true measured values (IMS scores)
    % - (2) predictedVals: matrix where each col is vector of predicted
    % vals corresponding to a ridge parameter k
    % Outputs:
    % - (1) NRMSE: normalized root mean-squared error of prediction


    % Compute RMSE for the main predictions
    error = predictedValsMat - repmat(trueVals, 1, size(predictedValsMat, 2));
    RMSE = sqrt(mean(error.^2, 1));
    
    % Find predictions of a mean-predictor (predicting each IMS as the mean 
    % of others)
    N = length(trueVals);
    meanPred = [];
    for i = 1:N
        meanPred(i, 1) = mean( trueVals(~ismember(1:N, i)) );
    end
    
    % Find NRMSE of the mean-predictor
    meanPredError = meanPred - trueVals;
    meanPredRMSE = sqrt(mean(meanPredError.^2, 1));
    
    % Compute NRMSE
    NRMSE = RMSE / meanPredRMSE;

end