function [ppparams,delfiles,keepfiles] = my_spmbatch_asl_M0correction(ppparams,params,delfiles,keepfiles)

if contains(params.asl.T1correctionM0,'tisssue_maps')
    preproc.channel.vols = {ppparams.rm0scan};
    preproc.channel.biasreg = 0.001;
    preproc.channel.biasfwhm = 60;
    preproc.channel.write = [0 0];
    preproc.tissue(1).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,1')};
    preproc.tissue(1).ngaus = 1;
    preproc.tissue(1).native = [1 0];
    preproc.tissue(1).warped = [0 0];
    preproc.tissue(2).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,2')};
    preproc.tissue(2).ngaus = 1;
    preproc.tissue(2).native = [1 0];
    preproc.tissue(2).warped = [0 0];
    preproc.tissue(3).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,3')};
    preproc.tissue(3).ngaus = 2;
    preproc.tissue(3).native = [1 0];
    preproc.tissue(3).warped = [0 0];
    preproc.tissue(4).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,4')};
    preproc.tissue(4).ngaus = 3;
    preproc.tissue(4).native = [0 0];
    preproc.tissue(4).warped = [0 0];
    preproc.tissue(5).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,5')};
    preproc.tissue(5).ngaus = 4;
    preproc.tissue(5).native = [0 0];
    preproc.tissue(5).warped = [0 0];
    preproc.tissue(6).tpm = {fullfile(spm('Dir'),'tpm','TPM.nii,6')};
    preproc.tissue(6).ngaus = 2;
    preproc.tissue(6).native = [0 0];
    preproc.tissue(6).warped = [0 0];
    preproc.warp.mrf = 1;
    preproc.warp.cleanup = 1;
    preproc.warp.reg = [0 0.001 0.5 0.05 0.2];
    preproc.warp.affreg = 'mni';
    preproc.warp.fwhm = 0;
    preproc.warp.samp = 3;
    preproc.warp.write = [0 1];
    preproc.warp.vox = NaN;
    preproc.warp.bb = [NaN NaN NaN;NaN NaN NaN];
    
    spm_preproc_run(preproc);
    
    gmmap = spm_file(ppparams.rm0scan, 'prefix','c1');
    wmmap = spm_file(ppparams.rm0scan, 'prefix','c2');
    csfmap = spm_file(ppparams.rm0scan, 'prefix','c3');

    [ppth,fname,~] = fileparts(ppparams.rm0scan);
     
    delfiles{numel(delfiles)+1} = {fullfile(ppth,[fname '._seg8.mat'])};
    delfiles{numel(delfiles)+1} = {fullfile(ppth,['y_' fname '.nii'])};  
    delfiles{numel(delfiles)+1} = {gmmap};
    delfiles{numel(delfiles)+1} = {wmmap};
    delfiles{numel(delfiles)+1} = {csfmap};
end

% The T1 values used, are the averaged T1 values reported in the review of 
% Bojorquez et al. 2017. What are normal relaxation times of tissues at 3 T? Magnetic Resonance Imaging 35:69-80
% (https://mri-q.com/uploads/3/4/5/7/34572113/normal_relaxation_times_at_3t.pdf)

T1gm = 1.459;
T1wm = 0.974;
T1csf = 4.190;

switch params.asl.T1correctionM0
    case 'tisssue_maps'
        
        GM = spm_vol(gmmap);
        WM = spm_vol(wmmap);
        CSF = spm_vol(csfmap);
        
        GMdat = spm_read_vols(GM);
        WMdat = spm_read_vols(WM);
        CSFdat = spm_read_vols(CSF);

        T1dat = (T1gm * GMdat + T1wm * WMdat + T1csf * CSFdat);
        
    case 'T1_map'
        % WIP: not yet implemented
    case 'average_GM'
        M0I = spm_vol(ppparams.rm0scan);
        m0dat = spm_read_vols(M0I);
        
        m0mask = my_spmbatch_mask(m0dat);

        m0mask(m0mask>0) = 1;

        T1dat = m0mask*T1gm;

    case 'average_WM'
        M0I = spm_vol(ppparams.rm0scan);
        m0dat = spm_read_vols(M0I);
        
        m0mask = my_spmbatch_mask(m0dat);

        T1dat = m0mask*T1wm;
end

m0dat = spm_read_vols(spm_vol(ppparams.rm0scan));

corr_T1 = 1 ./ (1-exp(-2.025./T1dat)); % M0 sscan @GE is a saturation rescovery SE with ST=2.0
m0dat = m0dat .* corr_T1;

VM0 = spm_vol(ppparams.rm0scan);

VM0.fname = spm_file(ppparams.rm0scan, 'prefix','t');
VM0.descrip = 'T1 corrected M0 scan';
VM0 = rmfield(VM0,'pinfo');
VM0 = spm_create_vol(VM0);
VM0 = spm_write_vol(VM0,m0dat);

delfiles{numel(delfiles)+1} = {spm_file(ppparams.rm0scan, 'prefix','t')};

ppparams.tm0scan = VM0.fname;