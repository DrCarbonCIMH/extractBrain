# extractBrain
Matlab function/app to extract rodent brain in MRI data based on PCNN algorithm.
Make sure you have SPM12 installed.

# Usage
[ Pout, finalBrainSize ] = ms_do_brainExtraction( P, BrSize, flag )
- P: MRI file
- BrSize: range in which the brain size is to be expected (e.g. [350,650]; multiply by 1000 if you "extend" voxelsize)
- flag: can have fields 'p', 'r', and 'preprocStyle' (see below for details

# Explanation
PCNN is an iterative procedure; it activates neurons/pixels based on
certain image properties (starting with the highest value and then influencing neighbours). 
the algorithm performs the PCNN stepwise and collects activated pixels; these pixels serve as mask;
idea: find the biggest cluster of activated pixels and define that as brain mask
you will see that the brain mask follows a decreasing slope till it
"breaks" and the brain mask combines with 'outer clusters' 
check the iterations with 'ms_app_checkBrainMasks'

** flag **
you can set some parameters using the 'flag' variable (see "set some
'standard' options" within the main file). Of most interest is probably..
- the p value which influences the smoothing of the mask; p=5 seems to be a good value -> but sometimes 3 is better (more details);
- r influences the PCNN itself; best to keep it 3 (or 2)
- preprocStyle is usually 'full' -> 'minimal' can be better if your image is already a nice one (e.g. bias-corrected, volume coil,..)

The best parameter set is probably more like an art.. but try to make the 'plateau' as wide as possible (see docs)

![alt text](https://github.com/DrCarbonCIMH/extractBrain/blob/main/doc/brainExtraction_result.png)
upper left: anatomical input image; upper right: preprocessed image (using 'preprocStyle' full); lower left: 'extracted' brain; lower right: estimated brain mask

In case you are not happy with the result you can use 'ms_app_checkBrainMasks' to define a new brain mask by choosing another iteration step. 'OptG' is the step which the algorithm chose initially. Load the created 'BrainMask*.mat' file (use SPM's recursive search function to load several .mat files at once).

![alt text](https://github.com/DrCarbonCIMH/extractBrain/blob/main/doc/checkBrainMasks_app.png)

!! Watch Out !! the given docs are not uptodate!
Use ms_app_checkBrainMasks.mlapp for checking the correct iterations

This code was inspired by the excellent paper "Robust Automatic Rodent Brain Extraction Using 3-D Pulse-Coupled Neural Networks (PCNN)" by
Nigel Chou, Jiarong Wu, Jordan Bai Bingren, Anqi Qiu, and Kai-Hsiang Chuang

When using this code, please cite [CIT]
