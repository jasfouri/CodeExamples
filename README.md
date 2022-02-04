# CodeExamples/DBS for TRD
Contains code from some of my work in Dr. Sameer Sheth's Lab. The aim is to decode mood state from neural activity (stereo-EEG signals) in treatment-resistant depression patients, eventually for closed-loop neurostimulation. The code performs dimensionality reduction on high-dimensional neural data then regresses mood state onto neural activity. It selects one model after fitting multiple candidate models which vary by region of neural activity and number of principal components retained. It uses leave-one-out cross-validation to evaluate the model due to the sparsity of mood state measurements, and a nested leave-one-out cross-validation stage to select the candidate model for each outer fold.

The script to run is "MoodDecoding.m". (It relies on my functions "modelAndPredict.m" and "evaluateCandidateModels.m".) "MoodDecoding.m" will output four figures. 
Figure 1 plots the prediction of mood state against true mood state using neural activity from one selected brain region (each test prediction datapoint corresponds to one fold of cross-validation): 

![001_1Region_ResultScatter](https://user-images.githubusercontent.com/98965657/152456252-d8cbf7ca-9f22-45e5-8493-b74465c7510a.jpg)

Figure 2 plots the results of a permutation test run with this model to evaluate the significance of its prediction accuracy:

![001_1Region_SignificanceTest](https://user-images.githubusercontent.com/98965657/152456335-5158a898-24e5-4485-80fa-20f9a880b375.jpg)

Figure 3 is a similar plot to Figure 1, except depicting the results of a model which uses two selected brain regions, performing worse:

![001_2Regions_ResultScatter](https://user-images.githubusercontent.com/98965657/152456519-e3993f19-fbcf-40a5-b65b-2f03d28a4078.jpg)

Figure 4 plots the magnitude of the regression coefficients for the power features back in the original feature space to interpret which features are most predictive of mood:

![001_1Region_CoeffSignif](https://user-images.githubusercontent.com/98965657/152456592-9024fe28-1896-4500-a597-1f7270f3a3ce.jpg)

The script can be modified to run on another dataset by changing the patient number on line 2 (options are "001" or "002").

# CodeExamples/AC Magnetometer
Contains code and data I collected from some of my work in Dr. Jacob Robinson's Lab. The code processes magnetic nanoparticle magnetization data (in the presence of an alternating magnetic field) to characterize the hysteresis loops of the nanoparticles. The file "HysteresisLoops.m" is the code and the folder "CodeExamples/AC Magnetometer/55kHz" contains the data. 

To collect a clean hysteresis loop from nanoparticles, I had to take two recordings of a nanoparticle sample and three recordings of a water sample (a control for further denoising the magnetization signal), then process all six combinations from these trials before selecting the cleanest hysteresis loop. "HysteresisLoops.m" yields three figures as output: 
Figure 1 plots the hysteresis loops colored by period of the magnetic field for each of the six combinations: 

![HysteresisLoopPeriodsEx](https://user-images.githubusercontent.com/98965657/152451283-ae753ae0-8d63-4d2d-bb54-a3a9c66256d5.jpg)


Figure 2 plots the average hysteresis loop from each combination after centering. Clearly, some loops are cleaner than others at the edges near saturation: 

![AvgHysteresisLoopsEx](https://user-images.githubusercontent.com/98965657/152451197-e2467b0b-f904-4f7d-b3bd-c300f8080543.jpg)


Figure 3 plots the average net magnetization signal (water minus nanoparticle) from a few periods of the recording: 

![AvgNetMEx](https://user-images.githubusercontent.com/98965657/152451208-f6fb9117-90ec-4e74-b3fb-6370e85bb2ee.jpg)

The script can be modified to run on data collected at a different magnetic field strength (change "voltage" on line 3), or data collected from other types of nanoparticles by selecting another one of the options under line 6.
