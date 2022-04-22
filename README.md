# extractBrain
Matlab function/app to extract rodent brain in MRI data based on PCNN algorithm.
Make sure you have SPM12 installed.

# Usage
[ Pout, finalBrainSize ] = ms_do_brainExtraction( P, BrSize, flag )

# Explanation
PCNN is an iterative procedure; it activates neurons/pixels based on
certain image properties (starting with the highest value and then influencing neighbours). 
the algorithm performs the PCNN stepwise and collects actived pixels; these pixels serve as mask;
idea: find the biggest cluster of activated pixels and define that as brain mask
you will see that the brain mask follows a decreasing slope till it
"breaks" and the brain mask combines with 'outer clusters' 
check the iterations with 'ms_gui_checkBrainMasks'

** flag **
you can set some parameters using the 'flag' variable (see "set some
'standard' options" within the main file). Of most interest is probably..
- the p value which influences the smoothing of the mask; p=5 seems to be a good value -> but sometimes 3 is better (more details);
- r influences the PCNN itself; best to keep it 3 (or 2)
- preprocStyle is usually 'full' -> 'minimal' can be better if your image is already a nice one (e.g. bias-corrected, volume coil,..)

The best parameter set is probably more like an art.. but try to make the 'plateau' as wide as possible (see docs)

![alt text](https://github.com/DrCarbonCIMH/extractBrain/blob/main/doc/brainExtraction_result.png)
![alt text](https://github.com/DrCarbonCIMH/extractBrain/blob/main/doc/checkBrainMasks_app.png)

!! Watch Out !! the given docs are not uptodate!
Use ms_app_checkBrainMasks.mlapp for checking the correct iterations

This code was inspired by the excellent paper "Robust Automatic Rodent Brain Extraction Using 3-D Pulse-Coupled Neural Networks (PCNN)" by
Nigel Chou, Jiarong Wu, Jordan Bai Bingren, Anqi Qiu, and Kai-Hsiang Chuang

When using this code, please cite [CIT]