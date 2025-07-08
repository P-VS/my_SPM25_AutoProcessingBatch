function ica_aroma_removal(param_file, sub, ses, noise_ICs, denoiseDir)
%ICA_AROMA_REMOVAL Removes the noise ICs 

% param_file: path to parameter initialization file created in do_ica
% noise_ICs: path of file with indices of the ICs classified as noise
% denoiseDir: directory where the denoised files are stored

mkdir(denoiseDir);

icatb_removeArtifact(param_file, denoiseDir, sub, ses, noise_ICs);

end


