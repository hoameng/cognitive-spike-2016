# cognitive-spike-2016
These scripts were used to detect and collate spikes and conduct the main statistical analysis for 

•	Ung H, Cazares C, Nanivadekar A, Kini L, Wagenaar J, Becker D, Krieger A, Lucas T, Litt B, Davis K. Interictal epileptiform activity outside the seizure onset zone impacts cognition. Brain (2017) 

Disclaimer: These scripts aren't tailored for general use, but feel free to modify to your liking. Uses the spike detector from Janca et al. 2014. This spike detector has a lot of false positives, and the main reason we were able to use the algorithm was because the experiment incorporated an internal control, and the signal in our experiment was found despite the noise. A better algorithm was developed based on the vetted spikes from this manuscript, though would require vetted example spikes (it is briefly detailed in my thesis https://repository.upenn.edu/edissertations/2618/)

# Files:
- run_detections_spikes_par.m	- Iterates thorough all IEEG datasets, detects spikes, runs spatial integration, and uploads to IEEG portal using parallel computing. Requires spike_ja_wrapper.m, spike_detector_hilbert_v16_byISARG.m, and utility scripts from portal-matlab-tools
https://github.com/ieeg-portal/portal-matlab-tools/tree/master/Utilities

- collate_spikes.m	- loads spikes above as well as experimental labels, collates to table for statistical analysis  
- run_stat_analysis.R - runs statistical analyses based on the spike tables from collate_spikes

Please feel free to contact me at hoameng.ung@gmail.com with any questions or comments.





