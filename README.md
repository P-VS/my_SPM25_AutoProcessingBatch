# my_SPM25_AutoProcessingBatch
UZB-fmr-processing-scripts

This repository contains the scripts used at UZ Brussel to standardise, optimise and automate the processing of single echo and multi-echo fMRI data. Additionally, a batch script is provided to do the segmentation in CAT12 for VBM. The python and Matlab scripts can be used independently from each other. The scripts come without any warrant for the results. You may use them at your own responsibility. I only share them to help and inspire others with their own processing scripts.

The folders contains following tools:

Python -convert_dcm2niix: convert dicom to nii data format using dcm2niix in Python and organise the data for the processing (BIDS based)

Matlab -my_spmbatch/AutoDCM2nii: initiating the automatic ccoonversion of dicom data into nifti data and organising them for further processing with my_spmbatch

-my_spmbatch/AutoSPMpreprocessing_fmri: initiating the automatic preprocessing of fMRI data.

-my_spmbatch/AutoSPM1stlevel_fmri: initiating the automatic 1st level analysis of fMRI data.

-my_spmbatch/AutoSPMpreprocessing_vbm: initiating the automatic preprocessing of VBM data (T1 anatomical high resolution scans) in CAT12.

-my_spmbatch/AutoSPMpreprocessing_asl: initiating of the automatic preprocessing of 3D PCASL data from a GE scanner (WIP)

-my_spmbatch/AutoSPMpreprocessing_aslbold: initiating the automatic preprocessing of ASLBOLD (experimental sequence on GE) (WIP)

Each analysis can be started by setting all the parameters in the beginning of the AutoSPM… script and clic ‘Run’.

fMRI_Class The slides of my fMRI classes are provided. During these classes, a demo dataset is analysed manually in the same way as is done by my scripts.
Required software installed:

Matlab
SPM25 (https://www.fil.ion.ucl.ac.uk/spm/)
MRIcroGL (https://www.mccauslandcenter.sc.edu/mricrogl/) for dicom to nifti conversion
CAT12 (https://neuro-jena.github.io/cat/) for VBM
GIFT (https://trendscenter.org/software/gift/) to use ICA-AROMA denoising

IMPORTANT: !! Look at your data before starting any preprocessing. It makes no sense to lose time in trying to preprocess bad data !!

Automatic set the origin in the anterior comisura and align with the MNI template to improve coregistration and normalisation steps
Remove dummy scans (prefix to the file name for the combination of steps 1 and 2 = e)
Realignment to the first echo (prefix to the file name = r)
TOPUP like EPI geometric distortion correction based on a phase encoding gradient polarity reversed scan (pepolar) (prefix to the file name = u)
For ASLBOLD: filetering between BOLD (f<0.1Hz) (prefix to the file name = f, endfix is set to _bold) and ALS (f>0.1Hz) part (endfix is set to _asl)
denoising (prefix to the file name = d)
Bandpass filtering and detrending
Extension of motion regressors to 24 by adding the temporal derivative and the squared regressors)
aCompCor (noise componenten determined on no gray or white matter voxels)
ICA-AROMA for single echo fMRI / ME-ICA-AROMA foor multi-echo fMRI. A component is labeled as noise if at least 1 of following criteria is met:
Fraction of the component in non brain areas > 2 * fraction of the component in brain areas (grey and white matter areas)
The highest correlation between a component’s time course and the noise regressors (24 motion+aCompCor) (see Van Schuerbeek te al. The optimized combination of aCompCor and ICA-AROMA to reduce motion and physiologic noise in task fMRI data. Biomed. Physiologic. Eng. Express 2022, 8(5), doi:10.1088/2057-1976/ac63f0)
Non T2*-related signal as determined with ME-ICA: Rho > 1.25 * Kappa (see Kundu et al. Differentiating BOLD and non-BOLD signals in fMRI time series using multi-echo EPI. NeuroImage 2012, 60(3): 1759-1770, code based ons tedana: https://github.com/ME-ICA/tedana)
High frequency content > 50%
Noise regression of the 24 motion regressors, aCompCor regressors / Soft removal of the ICA noise components
For ASLBOLD: the S0 components as determined with ME-ICA (non ASL if Kappa > 1.25 * Rho) in step 6.3 are added to the filtered ASL component from step 5 to form the functional asl signal
If multi-echo fMRI: echo combination (prefix to the file name = c). Choices are
Simple averaging
TE weighted
T2* weighted (as in tedana: https://github.com/ME-ICA/tedana)
Dynamic T2* mapping (see Heunis et al. 2021. The effects of multi-echo fMRI combination and rapid T2*-mapping on offline and real-time BOLD sensitivity. NeuroImage 238, 118244)
Slice time correction (also possible with HyperBand/Simultaneous multislice) (prefix to the file name = a)
Normalisation to the MNI template (prefix to the file name = w)
Smoothing (prefix to the file name = s)
Denoising (if not done in step 6) (prefix to the file name = d)
IMPORTANT: !! Look at your the result of the preprocessing before starting any 1st level (individual subject) or groups analysis. Doing statistical tests on wrongly preprocessed data can affect your results!
