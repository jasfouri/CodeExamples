amfFrequency = 55e3; 

voltage = "11";
%voltage options: "03", "05", "075", "1", "3", "5", "6", "7", "8", "9", "11", "13"

%%% Select one:
% sample = "FeOCluster1.25mL"; conc = "10.5";
sample = "FeOCluster1.25mL"; conc = "22";
% sample = "Cobalt"; conc = "10.5";
% sample = "Cobalt"; conc = "19";

samplechar = extractBetween(sample, 1, 1);
amfPeriod = 1/amfFrequency;
freq = num2str(amfFrequency/1e3) + "kHz";

allWaterData = []; allNpData = [];

PathName = cd;
FileName = freq + "/17uLWater/W" + voltage;
waterData1 = readtable(fullfile(PathName, FileName));
allWaterData(:,1,1) = waterData1.Source - min(waterData1.Source);
allWaterData(:,2,1) = waterData1.CH1;
allWaterData(:,3,1) = waterData1.CH2;

FileName = freq + "/17uLWater/WW" + voltage;
waterData2 = readtable(fullfile(PathName, FileName));
allWaterData(:,1,2) = waterData2.Source - min(waterData2.Source);
allWaterData(:,2,2) = waterData2.CH1;
allWaterData(:,3,2) = waterData2.CH2;

PathName = cd;
FileName = freq + "/17uL" + sample + conc + "/" + samplechar + voltage;
npData1 = readtable(fullfile(PathName, FileName));
allNpData(:,1,1) = npData1.Source - min(npData1.Source);
allNpData(:,2,1) = npData1.CH1;
allNpData(:,3,1) = npData1.CH2;

FileName = freq + "/17uL" + sample + conc + "/" + samplechar + samplechar + voltage;
npData2 = readtable(fullfile(PathName, FileName));
allNpData(:,1,2) = npData2.Source - min(npData2.Source);
allNpData(:,2,2) = npData2.CH1;
allNpData(:,3,2) = npData2.CH2;

allWaterData(:,1,3) = npData1.Source - min(npData1.Source);
allWaterData(:,2,3) = npData1.CH1;
allWaterData(:,3,3) = 0;

p = 1;
for w = 1:size(allWaterData, 3)
    for n = 1:size(allNpData, 3)
    waterTimeData = allWaterData(:,1,w);
    waterBData = allWaterData(:,2,w);
    waterMData = allWaterData(:,3,w);

    npTimeData = allNpData(:,1,n);
    npBData = allNpData(:,2,n);
    npMData = allNpData(:,3,n);

    samplingPeriod = (waterTimeData(end) - waterTimeData(1)) / length(waterTimeData);
    samplingFrequency = 1/samplingPeriod;

    %find at what timepoints zeros occur (i.e., timepoint of each half-period)
    
    %align waterBData and npBData
    q = 1; r = 1;
    iZerosWater = [];
    for k = 1:(length(waterTimeData)-1)
        if (waterBData(k) == 0) || (waterBData(k)*waterBData(k+1)<0)
            if ceil(waterTimeData(k)/(amfPeriod/2)) > r %assign zero values to new row for every half-period
                r = r + 1;
                q = 1;
            end
            iZerosWater(r, q) = k;
            q=q+1;
        end
    end

q = 1; r = 1;
    iZerosNP = [];
    for k = 1:(length(npTimeData)-1)
        if (npBData(k) == 0) || (npBData(k)*npBData(k+1)<0)
            if ceil(npTimeData(k)/(amfPeriod/2)) > r %assign zero values to new row for every half-period
                r = r + 1;
                q = 1;
            end
            iZerosNP(r, q) = k;
            q=q+1;
        end
    end
    
    iZerosWater(iZerosWater == 0) = NaN; iZeroWater = round(median(iZerosWater, 2, 'omitnan')); %indices of every half-period
    iZerosNP(iZerosNP == 0) = NaN; iZeroNP = round(median(iZerosNP, 2, 'omitnan'));
    iZeroWater = iZeroWater(~isnan(iZeroWater));
    firstZeroWater = find(waterBData(iZeroWater) == 0 & waterBData(iZeroWater + 1) > 0, 1);
    iZeroNP = iZeroNP(~isnan(iZeroNP));
    firstZeroNP = find(npBData(iZeroNP) == 0 & npBData(iZeroNP + 1) > 0, 1);
    
    waterBData = waterBData(iZeroWater(firstZeroWater):end);
    waterMData = waterMData(iZeroWater(firstZeroWater):end);
    npBData = npBData(iZeroNP(firstZeroNP):end); 
    npMData = npMData(iZeroNP(firstZeroNP):end);
    waterTimeData = waterTimeData(iZeroWater(firstZeroWater):end); waterTimeData = waterTimeData - min(waterTimeData);
    npTimeData = npTimeData(iZeroNP(firstZeroNP):end); npTimeData = npTimeData - min(npTimeData);
    
    if length(waterTimeData) > length(npTimeData)
        waterBData = waterBData(1:length(npTimeData)); waterMData = waterMData(1:length(npTimeData));
        waterTimeData = waterTimeData(1:length(npTimeData));  
        iZero = iZeroNP - iZeroNP(firstZeroNP) + 1; iZero(iZero < 0) = [];
    else
        npBData = npBData(1:length(waterTimeData)); npMData = npMData(1:length(waterTimeData));
        npTimeData = npTimeData(1:length(waterTimeData));
        iZero = iZeroWater - iZeroWater(firstZeroWater) + 1; iZero(iZero < 0) = [];
    end
    iPeriod = iZero(1:2:length(iZero));
    numPeriods = length(iPeriod) - 1;

    %center all data at 0V by subtracting average value
    RC_waterBData = waterBData - trapz(waterBData)/length(waterBData);
    RC_waterMData = waterMData - trapz(waterMData)/length(waterMData);
    RC_npBData = npBData - trapz(npBData)/length(npBData);
    RC_npMData = npMData - trapz(npMData)/length(npMData);

    netM = RC_npMData - RC_waterMData;
    time = linspace(0, samplingPeriod, length(waterTimeData));
    
    figure(3); subplot(size(allWaterData, 3), size(allNpData, 3), p); plot(time(iPeriod(1):iPeriod(10)), netM(iPeriod(1):iPeriod(10))); title('netM');
    
    M = cumtrapz(time, netM);
    B = cumtrapz(time, RC_npBData);

    RC_M = M - trapz(M)/length(M);
    RC_B = B - trapz(B)/length(B);

    %Re-adjust the signal
    FieldProbeArea = 6.698500E-06; % m^2 180305 Wire centered area
    %171024 calibrated FieldProbeArea = 0.000006668501646142370 % m^2
    Bfinal = RC_B/(20*FieldProbeArea);
    %Signal reduced by 20 times (from amplification) and divided by (area) * (number of turns = 1) => Tesla
    Mfinal = RC_M/(20)*246877.578285656;

    if mod(length(iZero), 2) == 1; iZero(end) = []; end
    iZeroR = reshape(iZero, 2, length(iZero)/2);

    lengthPeriods = diff(iPeriod); trueLengthPeriod = mode(lengthPeriods);
    BfinalR = []; MfinalR = []; q = 1;
    for i = 1:(length(iPeriod)-1)
        if lengthPeriods(i) == trueLengthPeriod
            BfinalR(q,:) = Bfinal(iPeriod(i):iPeriod(i+1));
            MfinalR(q,:) = Mfinal(iPeriod(i):iPeriod(i+1));
            q = q + 1;
        end
    end

    %plot loops without centering
    figure(1); subplot(size(allWaterData, 3), size(allNpData, 3), p); plot(BfinalR', MfinalR');
    colororder(jet(size(BfinalR,1)));
    cb = colorbar; caxis([1 size(BfinalR, 1)]); cb.Label.String = 'Period';
    title("W" + w + ", N" + n);

    %center each loop on 0 for each period
    periodAvgM = trapz(MfinalR, 2)/size(MfinalR, 2); periodAvgB = trapz(BfinalR, 2)/size(BfinalR, 2);
    MfinalRC = MfinalR - repmat(periodAvgM, 1, size(MfinalR, 2));
    BfinalRC = BfinalR - repmat(periodAvgB, 1, size(BfinalR, 2));

    %plot average loop using centered
    BfinalMean = mean(BfinalRC); MfinalMean = mean(MfinalRC);
    figure(2); subplot(size(allWaterData, 3), size(allNpData, 3), p); plot(BfinalMean, MfinalMean, 'LineWidth', 1);
    xlabel('Magnetic Field strength (T)'); ylabel('M (emu)')
    title("W" + w + ", N" + n);
    p = p + 1;
    end
end

function [] = plot_period(periodNum, iPeriod, numPeriods, Bfinal, Mfinal)
    if periodNum > numPeriods; periodNum = numPeriods; end
    periodInd = iPeriod(periodNum):iPeriod(periodNum+1);
    plot(Bfinal(periodInd), Mfinal(periodInd));
    xlabel('Magnetic Field strength (T)'); ylabel('M (emu)');
    title("Period # " + num2str(periodNum) + "/" + num2str(numPeriods));
end
