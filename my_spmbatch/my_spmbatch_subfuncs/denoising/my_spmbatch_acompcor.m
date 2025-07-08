function [ppparams,keepfiles] = my_spmbatch_acompcor(funcdat,ppparams,params,keepfiles)

if isfield(ppparams,'der_file')
    confounds = load(ppparams.der_file);
    [~,cfname,~] = fileparts(ppparams.der_file);
elseif isfield(ppparams,'rp_file')
    confounds = load(ppparams.rp_file);
    [~,cfname,~] = fileparts(ppparams.rp_file);
else
    confounds = [];
    tcfname = [ppparams.func(ppparams.echoes(1)).prefix ppparams.func(ppparams.echoes(1)).funcfile];
    cfname = split(tcfname,'.nii');
    cfname = cfname{1};
end

GM = spm_vol(ppparams.fc1im);
WM = spm_vol(ppparams.fc2im);
CSF = spm_vol(ppparams.fc3im);

gmdat = spm_read_vols(GM);
wmdat = spm_read_vols(WM);
csfdat = spm_read_vols(CSF);

braindat = gmdat+wmdat;
braindat(braindat<0.2)=0;
braindat(braindat>0.0)=1;

csfdat(braindat>0.0)=0;
csfdat(csfdat<0.2)=0;
csfdat(csfdat>0.0)=1;

acc_confounds = fmri_acompcor(funcdat(:,:,:,:),{csfdat},params.denoise.Ncomponents,'confounds',confounds,'filter',[],'PolOrder',1);

if ~isempty(confounds)
    new_confounds = cat(2,confounds,acc_confounds);
else
    new_confounds = acc_confounds;
end

accf = ['acc_' cfname];
accname = split(accf,'_echo-');
ppparams.acc_file = fullfile(ppparams.subfuncdir,[accname{1} '.txt']);

writematrix(new_confounds,ppparams.acc_file,'Delimiter','tab');

keepfiles{numel(keepfiles)+1} = {ppparams.acc_file};

clear csfdat braindat gmdat wmdat