function [ppparams,delfiles,keepfiles] = my_spmbatch_aslbold_normalization(ppparams,params,delfiles,keepfiles)

Vcbf = spm_vol(fullfile(ppparams.subperfdir,[ppparams.perf(1).cbfprefix ppparams.perf(1).cbffile]));

for i=1:numel(Vcbf)
    waslfiles{i,1} = [Vcbf(i).fname ',' num2str(i)];
end

%% Normalization of the M0 scan

m0normest.subj.vol = {fullfile(ppparams.subperfdir,[ppparams.perf(1).m0scanprefix ppparams.perf(1).m0scanfile ',1'])};
m0normest.eoptions.biasreg = 0.0001;
m0normest.eoptions.biasfwhm = 60;
m0normest.eoptions.tpm = {fullfile(spm('Dir'),'tpm','TPM.nii')};
m0normest.eoptions.affreg = 'mni';
m0normest.eoptions.reg = [0 0 0.1 0.01 0.04];
m0normest.eoptions.fwhm = 0;
m0normest.eoptions.samp = 3;

spm_run_norm(m0normest);

ppparams.deffile = fullfile(ppparams.subperfdir,['y_' ppparams.perf(1).m0scanprefix ppparams.perf(1).m0scanfile ',1']);
delfiles{numel(delfiles)+1} = {ppparams.deffile};

%% Normalise CBF data

%Write the spatially normalised data

cbfnormw.woptions = spm_get_defaults('normalise.write');
cbfnormw.woptions.vox = params.func.normvox;

dt = Vcbf(1).dt;
if dt(1)==spm_type('uint16')
    cbfnormw.woptions.interp = 4;
else
    cbfnormw.woptions.interp = 1;
end

cbfnormw.subj.def = {ppparams.deffile};
cbfnormw.subj.resample = waslfiles(:,1);

spm_run_norm(cbfnormw);

keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subperfdir,['w' ppparams.perf(1).cbfprefix ppparams.perf(1).cbffile])};

ppparams.perf(1).wcbffile = ['w' ppparams.perf(1).cbfprefix ppparams.perf(1).cbffile];

clear Vcbf
