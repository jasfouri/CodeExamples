%% Load Channel Data and Identify Regions
patient = "001";
projectdir = cd;
electrodes = readtable(fullfile(projectdir, "Electrodes_" + patient + ".csv"));

get_channels_area = @(r) electrodes.Electrode((electrodes.area == r) & (electrodes.Label ~= "empty") & (electrodes.Electrode ~= 156))';
ofc_channels = get_channels_area("OFC");
acc_channels = get_channels_area("ACC");
dpfc_channels = get_channels_area("DPFC");
vpfc_channels = get_channels_area("VPFC");

channels = sort([ofc_channels, acc_channels, vpfc_channels, dpfc_channels]);
regions = {ofc_channels, vpfc_channels, dpfc_channels, acc_channels};

%% Load CAT-DI Data
catdi = readtable(fullfile(cd, "CATDI_scores_" + patient + ".csv"));

if patient == "001"
    scores = str2double(extractBetween(catdi.Result, 1, 4)); %incomplete scores are NaN
    omitBlocks = [25, 26, 32, 34, 38, 39];
elseif patient == "002"
    scores = catdi.Result;
    scores(13) = []; %no data for this block
    omitBlocks = [2, 14, 15, 19];
end
scores(omitBlocks) = NaN;

%% Load neural data already processed
datadir = fullfile(cd, "LoadedData_" + patient);
if patient == "001"
    load(fullfile(datadir, "allFeatures_fs=500Dec_Twindow=5_numBands=6.mat"));
end

if patient == "002"
    load(fullfile(datadir, "allFeatures_fs=500Dec_Twindow=5_numBands=6.mat"));
end

features(isnan(scores)) = [];
scores(isnan(scores)) = [];

%% Select candidate single-region models for each fold of CV
PCOptions = [1 5 10 15]; %can increase this higher because we have more CATDI testpoints
[trueTestNRMSE1, selectedModel1, coeffMag1, PAll, TAll] = modelAndPredict(scores, features, PCOptions, channels, regions, []);
%% Assess the resulting NRMSE with random and permuted tests
% minScore = min(scores); maxScore = max(scores);
% % numRandTests = 1000; randTestNRMSE = zeros(length(numRandTests), 1);
% % for i = 1:numRandTests
% %     randScores = round(rand(size(scores)) * (maxScore - minScore) + minScore); %need to be integers
% %     randTestNRMSE(i) = modelAndPredict(randScores, features, PCOptions, channels, regions, []);
% %     disp(i)
% % end
load(fullfile(datadir,'1000randTestNRMSE_20210802')) %saved after execution of the above loop
figure; histogram(randTestNRMSE, 'Normalization', 'pdf'); hold on

randTestDist = fitdist(randTestNRMSE', 'Gamma');
rangeRandNRMSE = floor(min(randTestNRMSE)):1e-6:ceil(max(randTestNRMSE));
randTestPDF = pdf(randTestDist, rangeRandNRMSE);
[~, pADGamma] = adtest(randTestNRMSE, 'Distribution', randTestDist);
gammaFit = (pADGamma < 0.05);
if ~gammaFit
    randTestDist = fitdist(mink(randTestNRMSE, 250), 'GeneralizedPareto');
    randTestPDF = pdf(randTestDist, rangeRandNRMSE);
    [~, pADPareto] = adtest(randTestNRMSE, 'Distribution', randTestDist);
    paretoFit = (pADPareto < 0.05);
end
plot(rangeRandNRMSE, randTestPDF); xline(trueTestNRMSE1, 'LineWidth', 2, 'Color', 'r');
xlabel('Test Prediction Error (NRMSE)'); ylabel('Probability Density');
% xlim([min(rangeRandNRMSE) max(rangeRandNRMSE)]);

randTestCDF = cdf(randTestDist, rangeRandNRMSE);
pRandTest = randTestCDF(rangeRandNRMSE == round(trueTestNRMSE1, 4));
title("p = " + num2str(pRandTest))

[~, iCrit] = min(abs(randTestCDF - 0.05)); sigNRMSE = rangeRandNRMSE(iCrit);
disp("p Value: " + num2str(pRandTest));
disp("True test NRMSE: " + num2str(trueTestNRMSE1) + "; Required NRMSE for significance: " + num2str(sigNRMSE));

%% Select candidate two-region models
figure; [trueTestNRMSE2, selectedModel2] = modelAndPredict(scores, features, PCOptions, channels, regions, selectedModel1(:,2));

%% Select candidate three-, four-, or five-region models
% [trueTestNRMSE3, selectedModel3] = modelAndPredict(scores, features, PCOptions, channels, regions, [selectedModel1(:,2), selectedModel2(:,2)]);
% [trueTestNRMSE4, selectedModel4] = modelAndPredict(scores, features, PCOptions, channels, regions, [selectedModel1(:,2), selectedModel2(:,2), selectedModel3(:,2)]);
% [trueTestNRMSE5, selectedModel5] = modelAndPredict(scores, features, PCOptions, channels, regions, [selectedModel1(:,2), selectedModel2(:,2), selectedModel3(:,2), selectedModel4(:,2)]);

%% Determine which channels and frequencies are most significant predictors
num_bands = 6;
figure; selected_r = mode(selectedModel1(:,2)); selected_num = mode(selectedModel1(:,1));
coeffMagArray = vertcat(coeffMag1{selectedModel1(:,1) == selected_num});
imagesc(abs(coeffMagArray)); hold on
for i = 1:(length(regions{selected_r})-1)
    plot([(i*num_bands+ 0.5), (i*num_bands + 0.5)], [0, (length(scores)+0.5)], 'k-', 'LineWidth', 2.5);
end

ylabel('Fold of cross-validation'); xlabel('Feature (Freq bands of channels of selected region)');
h = colorbar;
ylabel(h, '|Regression coefficient|');