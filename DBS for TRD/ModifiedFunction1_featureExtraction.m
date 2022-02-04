% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Omid Sani, Yuxiao Yang, Maryam Shanechi
%   Shanechi Lab, University of Southern California, 2018
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function features = ModifiedFunction1_featureExtraction(rawData, dataFs, T)

    % Extracts power features from neural signals
    % Inputs:
    % - (1) rawData (time x channel): ECoG data matrix 
    % - (2) dataFs (scalar number): Sampling rate of ECoG, in Hz (e.g. 256)

    % Outputs:
    % - (1) features (time x channel x band): Signal log power extracted
    %           from non-overlapping 10s windows, in each of the 5 frequency
    %           bands specified below (see "bands")
    
%     bands = [1 4; 4 8; 8 12; 12 30; 35 50; 70 124]; % 6 frequency bands - for ds to 250Hz
    bands = [1 4; 4 8; 8 12; 12 30; 35 50; 70 150]; % 6 frequency bands - for ds to 500Hz
    
%     bands = [1 8; 8 12; 12 30; 35 50; 70 124]; %5 freq bands (combined delta and theta) - for ds to 250Hz
%     bands = [1 8; 8 12; 12 30; 35 50; 70 150]; %5 freq bands (combined delta and theta) - for ds to 500Hz

%     bands = [1 8; 8 12; 12 30; 30 55; 65 100]; % Start-end frequency of bands - Shanechi
      
    timeStep = T; % Time step of feature extraction [in seconds]
    
    % Compute the number of T s windows that fit in the data
    stepSamples = round(dataFs * timeStep); 
    totalSteps = floor(size(rawData, 1)/stepSamples);
    totalDataSamples = totalSteps*stepSamples;

    features = [];
    for bi = 1:size(bands, 1) % Iterate over feature frequency bands
        % %%%%%%%%%%%%%%%%
        % Design filter
        % %%%%%%%%%%%%%%%%
        band = bands(bi, :);
        
        N   = 8;  % Order
        Fc1 = band(1);  % First Cutoff Frequency
        Fc2 = band(2);  % Second Cutoff Frequency

        % Construct an FDESIGN object and call its BUTTER method.
        h  = fdesign.bandpass('N,F3dB1,F3dB2', N, Fc1, Fc2, dataFs);
        Hd = design(h, 'butter');
        [b,a] = sos2tf(Hd.sosMatrix, Hd.ScaleValues);
        % %%%%%%%%%%%%%%%%
        % End of filter design
        % %%%%%%%%%%%%%%%%
        
        % %%%%%%%%%%%%%%%%
        % Apply filter
%         spectrum1 = abs(fft(rawData(:,1))); figure; plot(spectrum1(1:floor(length(rawData(:,1))/2))); hold on; 
%         [mag, ~] = freqz(b,a, length(rawData(:,1)), dataFs); plot(abs(mag)*10^(4.5));
%         filteredRawData = FiltFiltM(b, a, rawData);
        filteredRawData = filtfilt(b, a, rawData);

        % Separate non-overlapping T s windows in the data and put them in
        % columns of a matrix (time X window number X channel)
        filteredRawData = reshape(filteredRawData(1:totalDataSamples, :), [stepSamples, totalSteps, size(filteredRawData, 2)]);

        % Compute root mean-squared energy in each window
        % and then take db (i.e. 20*log10(x))
        dBPowerInBand = mag2db(rms(filteredRawData, 1) + eps);
        
        % Reshape data into a window number x channels
        dBPowerInBand = permute(dBPowerInBand, [2 3 1]);

        features = cat(3, features, dBPowerInBand);
    end
end