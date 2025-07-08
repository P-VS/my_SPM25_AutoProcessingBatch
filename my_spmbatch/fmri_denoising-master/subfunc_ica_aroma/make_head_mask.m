function head_mask = make_head_mask(anat)

%MAKE_HEAD_MASK(ANAT) creates a head mask of the anatomical image 
%Inputs:
%   -anat: directory of the anatomical file
%Returns: 
%   -head_mask

Vanat = spm_vol(anat);
anatdat = spm_read_vols(Vanat);

% SPM tends to put NaNs in the data outside the brain
anatdat(~isfinite(anatdat)) = 0;
sorted_anatdat = sort(anatdat(:));% sort ascending

%calculate lower and upper cutoff indexes
lower_cutoff = floor(0.2 * length(sorted_anatdat));
upper_cutoff = min(floor(0.85 * length(sorted_anatdat)), length(sorted_anatdat));
%Calculate the differences between adjacent elements of the sorted input array
%between the given indices. 
%Each element delta[i] is the difference between sorted_anatdat(i+1) and
%sorted_anatdat(i) for i ranging from lower_cutoff to upper_cutoff 

delta = sorted_anatdat(lower_cutoff + 1 : upper_cutoff) - sorted_anatdat(lower_cutoff:upper_cutoff -1);

[~, ia] = max(delta); %index of the maximum element of delta
% ia: position where the derivative of the sorted_anatdat array changes most rapidly. 

%threshold: average of the two adjacent elements around the ia-th index of sorted_anatdat, 
%which is the midpoint of the steepest ascent in the array.

threshold = 0.5 *(sorted_anatdat(ia + lower_cutoff) + sorted_anatdat(ia + lower_cutoff + 1));

hmask = anatdat >= threshold;

% Opening 
SE = strel('square', 3);
hmask_open = imopen(hmask, SE);
% Closing 
SE = strel('square', 2);
head_mask = imclose(hmask_open, SE);