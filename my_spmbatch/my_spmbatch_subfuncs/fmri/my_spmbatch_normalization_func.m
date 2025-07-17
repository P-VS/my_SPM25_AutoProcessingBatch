function [ppparams,delfiles,keepfiles] = my_spmbatch_normalization_func(ne,ppparams,params,delfiles,keepfiles)

Vfunc = spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(ne).prefix ppparams.func(ne).funcfile]));

for i=1:numel(Vfunc)
    wfuncfiles{i,1} = [Vfunc(i).fname ',' num2str(i)];
end

%% Normalization of the functional scan
if ne==ppparams.echoes(1)
    ppparams.deffile = fullfile(ppparams.subfuncdir,['y_' ppparams.reffunc]);
    if ~exist(ppparams.deffile,"file")
        reffile = fullfile(ppparams.subfuncdir,ppparams.reffunc);
    
        funcnormest.subj.vol = {reffile}; % {fullfile(ppparams.subfuncdir,[ppparams.func(ne).prefix ppparams.func(ne).funcfile ',1'])};
        funcnormest.eoptions.biasreg = 0.0001;
        funcnormest.eoptions.biasfwhm = 60;
        funcnormest.eoptions.tpm = {fullfile(spm('Dir'),'tpm','TPM.nii')};
        funcnormest.eoptions.affreg = 'mni';
        funcnormest.eoptions.reg = [0 0 0.1 0.01 0.04];
        funcnormest.eoptions.fwhm = 0;
        funcnormest.eoptions.samp = 3;
    
        spm_run_norm(funcnormest);
    end

    delfiles{numel(delfiles)+1} = {ppparams.deffile};
end

%% Normalise func

%Write the spatially normalised data

funcnormw.woptions = spm_get_defaults('normalise.write');
funcnormw.woptions.vox = params.func.normvox;

dt = Vfunc(1).dt;
if ~(dt(1)==spm_type('uint16'))
    funcnormw.woptions.interp = 1;
end

funcnormw.subj.def = {ppparams.deffile};
funcnormw.subj.resample = wfuncfiles(:,1);

spm_run_norm(funcnormw);

keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subfuncdir,['w' ppparams.func(ne).prefix ppparams.func(ne).funcfile])};

ppparams.func(ne).prefix = ['w' ppparams.func(ne).prefix];

clear Vfunc Vtemp