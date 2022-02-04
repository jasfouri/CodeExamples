# CodeExamples
"CodeExamples/AC Magnetometer" contains code and data I collected from some of my work in Jacob Robinson's Lab. The code processes magnetic nanoparticle magnetization data (in the presence of an alternating magnetic field) to characterize the hysteresis loops of the nanoparticles. The file "HysteresisLoops.m" is the code and the folder "CodeExamples/AC Magnetometer/55kHz" contains the data. 

It is difficult to collect a clean hysteresis loop from the nanoparticles. I took two recordings of a nanoparticle sample and three recordings of a water sample (a control for further denoising the magnetization signal), then processed all six combinations from these trials to select the cleanest hysteresis loop. "HysteresisLoops.m" yields three figures as output: 
Figure 1 plots the hysteresis loops colored by period of the magnetic field for each of the six combinations: ![HysteresisLoopPeriodsEx](https://user-images.githubusercontent.com/98965657/152451283-ae753ae0-8d63-4d2d-bb54-a3a9c66256d5.jpg)


Figure 2 plots the average hysteresis loop from each combination after centering. Clearly, some loops are cleaner than others at the edges near saturation: ![AvgHysteresisLoopsEx](https://user-images.githubusercontent.com/98965657/152451197-e2467b0b-f904-4f7d-b3bd-c300f8080543.jpg)


Figure 3 plots the average net magnetization signal (water minus nanoparticle) from a few periods of the recording: ![AvgNetMEx](https://user-images.githubusercontent.com/98965657/152451208-f6fb9117-90ec-4e74-b3fb-6370e85bb2ee.jpg)
