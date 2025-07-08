function [nonBOLDICdata,nonASLICdata] = ica_aroma_classification(ppparams, params, funcmask, ica_dir, t_r)

%This function classifies the ica components extracted by fmri_do_ica into
% noise/no-noise (both physiological and motion noise) based on:
%- The maximum robust correlation with the noise regressors (confounds):
%         - aCompCor confounds (physiological noise components) identified by aCompCor (the temporal signals of 5PCAcomponents determined
%           in the CSF using aCompCor, as physiological noise regressors.)
%         - 24 the motion regressors computed from the realignment
%           parameters
%- The amount of high frequency (> 1 Hz)
%- Their location at the brain edge

aroma_dir = fullfile(ica_dir,'ICA-AROMA');
if ~exist(aroma_dir,"dir"), mkdir(aroma_dir); end

GM = ppparams.fc1im;
WM = ppparams.fc2im;
CSF = ppparams.fc3im;

gmdat = spm_read_vols(spm_vol(GM));
wmdat = spm_read_vols(spm_vol(WM));
csfdat = spm_read_vols(spm_vol(CSF));

braindat = spm_read_vols(spm_vol(funcmask));
braindat((gmdat + wmdat) < 0.2) = 0;
braindat(braindat > 0.0) = 1;
    
Vbr = spm_vol(funcmask);
Vbr.fname = fullfile(aroma_dir, 'braindata.nii'); 
Vbr.descrip = 'braindata';
Vbr = rmfield(Vbr, 'pinfo'); %remove pixel info so that there is no scaling factor applied so the values
spm_write_vol(Vbr, braindat);

nobrain = spm_read_vols(spm_vol(funcmask)); 
nobrain(braindat > 0.0) = 0; 
nobrain(nobrain > 0.0) = 1;

Vcsf = spm_vol(funcmask);
Vcsf.fname = fullfile(aroma_dir, 'nobraindata.nii'); 
Vcsf.descrip = 'nobraindata';
Vcsf = rmfield(Vcsf, 'pinfo'); %remove pixel info so that there is no scaling factor applied so the values
spm_write_vol(Vcsf, nobrain);

icaparams_file = fullfile(ica_dir,'ica_aroma_ica_parameter_info.mat');

load(icaparams_file);

%%%%%%%%%% Get the required variables from sesInfo structure %%%%%%%%%%
% Number of subjects
numOfSub = sesInfo.numOfSub;
numOfSess = sesInfo.numOfSess;

% Number of components
numComp = sesInfo.numComp;

dataType = sesInfo.dataType;

mask_ind = sesInfo.mask_ind;

% First scan
structFile = deblank(sesInfo.inputFiles(1).name(1, :));

%%%%%%%% End for getting the required vars from sesInfo %%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% Get component data %%%%%%%%%%%%%%%%%

% Get the ICA Output files
icaOutputFiles = sesInfo.icaOutputFiles;

[subjectICAFiles, meanICAFiles, tmapICAFiles, meanALL_ICAFile] = icatb_parseOutputFiles('icaOutputFiles', icaOutputFiles, 'numOfSub', ...
        numOfSub, 'numOfSess', numOfSess, 'flagTimePoints', sesInfo.flagTimePoints);

% component files
if ~exist(meanALL_ICAFile.name)
    compFiles = subjectICAFiles(1).ses(1).name;
else
    compFiles = meanALL_ICAFile.name;
end

compFiles = icatb_fullFile('directory', ica_dir, 'files', compFiles);

compData = spm_read_vols(spm_vol(compFiles));
dim = size(compData);
HInfo.V = spm_vol(structFile);
HInfo.DIM = dim(1:3);
HInfo.VOX = double(HInfo.V(1).private.hdr.pixdim(2:4)); HInfo.VOX = abs(HInfo.VOX);

compData = permute(compData, [4 1 2 3]);
compData = reshape(compData, size(compData, 1), prod(dim(1:3)));
[compData] = icatb_applyDispParameters(compData, 1, 3, 1, dim(1:3), HInfo);
compData = reshape(compData, [size(compData,1), dim(1), dim(2), dim(3)]);
compData = permute(compData, [2 3 4 1]);

if ppparams.save_intermediate_results
    VC = spm_vol(compFiles);
    for ic=1:numComp
        VC(ic) = spm_vol(funcmask);
        VC(ic).fname = fullfile(aroma_dir, 'thres_compdata.nii'); 
        VC(ic).n = [ic 1];
        spm_write_vol(VC(ic), compData(:,:,:,ic));
    end
end

% load time course
icaTimecourse = icatb_loadICATimeCourse(compFiles, 'real', [], [1:numComp]);

% Structural volume
dim = HInfo.DIM;
tdim = size(icaTimecourse(:,1));

% Reshape to 2d
compData = reshape(compData, [prod(dim(1:3)),numComp]);

%%%%%%%% End for getting the component data %%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%% Decision tree %%%%%%%%%%%%%%%%%

BOLDComp = ones([numComp,1]);
ASLComp = ones([numComp,1]);

% 1. """Fraction component outside GM or WM""
fprintf('Brain/No brain fractions \n') 

nobrainFract = zeros(numComp,1);
brainFract = zeros(numComp,1);

for i = 1:numComp       
    
    Compdat = compData(:,i);

    totComp = sum(Compdat(Compdat>0),"all");
    brainComp = sum(Compdat(braindat>0),"all");
    nbrainComp = sum(Compdat(nobrain>0),"all");

    nobrainFract(i) = nbrainComp/totComp;
    brainFract(i) = brainComp/totComp;

    clear Compdat
end

% Componets which falls mainly outside the brain are considered as noise
tmp = find(nobrainFract > 3*brainFract);
if ~isempty(tmp)
    BOLDComp(tmp) = 0;
    ASLComp(tmp) = 0;
end

%  2. """Maximum robust correlation with confounds"""  
fprintf('Correlation with noise regresors \n') 

%     Read confounds
    if isfield(ppparams,'acc_file'), conffile = ppparams.acc_file; 
    elseif isfield(ppparams,'der_file'), conffile = ppparams.der_file; 
    elseif isfield(ppparams,'rp_file'), conffile = ppparams.rp_file; end

    conf_model = load(conffile); 

%     Determine the maximum correlation between confounds and IC time-series
    [nmixrows, nmixcols] = size(icaTimecourse);
    [nconfrows, nconfcols] = size(conf_model);

    corr_mat = zeros(nmixcols,nconfcols);
    for icomp=1:nmixcols
        iicaTC = icaTimecourse(:,icomp)-mean(icaTimecourse(:,icomp),'all');
        iicaTC = iicaTC / std(iicaTC);

        for inoise=1:nconfcols
            iConf = conf_model(:,inoise)-mean(conf_model(:,inoise),'all');
            iConf = iConf / std(iConf);

            ijcorr = iicaTC .* iConf;
            corr_mat(icomp,nconfcols) = sum(ijcorr,'all')/nmixrows;
        end
    end
    
%     Resample the mix and conf_model matrices to have the same number of columns
    
    %corr_mat = corr(nmix, nconf_model);
    max_correls = max(corr_mat, [], 2);
    
    max_correls = double(reshape(max_correls,[nmixcols,int64(size(max_correls, 1)/nmixcols)]));
    maxRPcorr = max_correls(:,1);

tmp = find(abs(maxRPcorr) > 0.75);
if ~isempty(tmp)
    BOLDComp(tmp) = 0;
    ASLComp(tmp) = 0;
end

% 3.    if ME-fMRI do MEICA
if numel(ppparams.echoes)>2
    fprintf('MEICA \n') 

    [kappas,rhos,kappa_elbow,rhos_elbow] = my_spmbatch_meica(compData,icaTimecourse,ppparams,params);

    tmp = find(and(3.0*kappas<rhos,kappas<kappa_elbow));
    if ~isempty(tmp), BOLDComp(tmp) = 0; end

    tmp = find(and(3.0*rhos<kappas,rhos<rhos_elbow));
    if ~isempty(tmp), ASLComp(tmp) = 0; end
else
    kappas=zeros([numComp,1]);
    rhos=zeros([numComp,1]);
    kappa_elbow=zeros([numComp,1]);
    rhos_elbow=zeros([numComp,1]);
end

% 4. """High frequency content"""
fprintf('High frequency content \n') 

FT = abs(fft(icaTimecourse, [], 1));% FT of each IC along the firt dimension (time) 
FT = FT(1:(length(FT)/2) +1,:); % keep postivie frequencies (Hermitian symmetric) +1 to get Nyquist frequency bin as well 
       
%   Determine sample frequency
    Fs = 1/t_r;
    
%   Determine Nyquist-frequency
    Ny = Fs/2;
    
%   Determine which frequencies are associated with every row in the melodic_FTmix file  (assuming the rows range from 0Hz to Nyquist)
    f = Ny * (1 : size(FT, 1)) / size(FT, 1);
    
%   Only include frequencies higher than 0.01Hz
    fincl = find(f > 0.01); %get indices
    FT = FT(fincl, :);
    f = f(fincl);
    
%     Set frequency range to [0-1]
    f_norm = (f - 0.01) / (Ny - 0.01);
    
%     For every IC; get the cumulative sum as a fraction of the total sum
    fcumsum_fract = cumsum(FT,1) ./ sum(FT,1);
    
%     Determine the index of the frequency with the fractional cumulative sum closest to 0.5
    [~, idx_cutoff] = min(abs(fcumsum_fract - 0.5));
    
%     Now get the fractions associated with those indices, these are the final feature scores
    HFC = f_norm(idx_cutoff)';
    
    thr_HFC = f_norm(min(find((f-0.12)>0)));

    tmp = (HFC > thr_HFC);
    if ~isempty(tmp), BOLDComp(tmp) = 0; end

%   Classify the ICs

    nonBOLDICs = find(BOLDComp<1);
    nonBOLDICdata = icaTimecourse(:,nonBOLDICs);

    nonASLICs = find(ASLComp<1);
    nonASLICdata = icaTimecourse(:,nonASLICs);
    
%   Save results

    if ppparams.save_intermediate_results
        if numel(ppparams.echoes)>1
            if params.func.isaslbold
                varnames = {'NoBrain_fraction','Brain_fraction','max_correlations','kappa','rho','high_frequency_content','BOLD_comp','ASL_comp'};
                T = table(nobrainFract,brainFract,maxRPcorr,kappas,rhos,HFC,BOLDComp,ASLComp,'VariableNames',varnames);
            else
                varnames = {'NoBrain_fraction','Brain_fraction','max_correlations','kappa','rho','high_frequency_content','BOLD_comp'};
                T = table(nobrainFract,brainFract,maxRPcorr,kappas,rhos,HFC,BOLDComp,'VariableNames',varnames);
            end
        else
            if params.func.isaslbold
                varnames = {'NoBrain_fraction','Brain_fraction','max_correlations','high_frequency_content','BOLD_comp','ASL_comp'};
                T = table(nobrainFract,brainFract,maxRPcorr,HFC,BOLDComp,ASLComp,'VariableNames',varnames);
            else
                varnames = {'NoBrain_fraction','Brain_fraction','max_correlations','high_frequency_content','BOLD_comp'};
                T = table(nobrainFract,brainFract,maxRPcorr,HFC,BOLDComp,'VariableNames',varnames);
            end
        end
    
        aroma_file = fullfile(aroma_dir,'AROMA_desission.csv');
    
        writetable(T,aroma_file,'WriteRowNames',false);
               
        funcfname = split([ppparams.func(ppparams.echoes(1)).prefix ppparams.func(ppparams.echoes(1)).funcfile],'.nii');
    
        if ~params.use_parallel
            fg = spm_figure('FindWin','Graphics');
        
            for icomp=1:numComp
                spm_figure('Clear','Graphics');
                figure(fg);
        
                plot([0:t_r:t_r*(tdim(1)-1)],icaTimecourse(:,icomp))
                
                saveas(fg,fullfile(aroma_dir,['comp-' num2str(icomp,'%03d') '_time.png']));
        
                spm_figure('Clear','Graphics');
                figure(fg);
                
                plot(f,FT(:,icomp)./max(FT(:,icomp),[],'all'),f,fcumsum_fract(:,icomp))
                xline(thr_HFC*Ny)
                
                saveas(fg,fullfile(aroma_dir,['comp-' num2str(icomp,'%03d') '_frequency.png']));
            end
        end
        
        % Path to noise ICs (might want to change path)
        noise_ICs_dir = fullfile(aroma_dir, strcat('ICA-AROMA_ICs_noise_',funcfname{1},'.txt'));
        
        if ~isempty(nonBOLDICs) 
            if length(size(nonBOLDICs)) > 0
                dlmwrite(noise_ICs_dir, nonBOLDICs(:), 'precision', "%i"); %write matrix to text file with values separated by a ',' as a delimiter
            else
                dlmwrite(noise_ICs_dir, int64(nonBOLDICs), 'precision', "%i");
            end
        end
    end
end