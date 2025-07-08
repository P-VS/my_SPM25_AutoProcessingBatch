function [ppparams,delfiles,keepfiles] = my_spmbatch_asl_cbfquantification(ppparams,delfiles,keepfiles)

try
    jsondat = fileread(ppparams.deltamjson);
    jsondat = jsondecode(jsondat);
    
    LabelingDuration = jsondat.LabelingDuration;
    PLD = jsondat.PostLabelDelay;
    NumberOfAverages = jsondat.NumberOfAverages;
catch
    LabelingDuration = 1.5;
    PLD = 1.525;
    NumberOfAverages = 1;
end

%% Step 1: make head mask

M0I = spm_vol(ppparams.tm0scan);
m0dat = spm_read_vols(M0I);

m0mask = my_spmbatch_mask(m0dat);

% from Alsop 2015 MRM
BloodT1 = 1.650; 
LabelingEfficiency = 0.85; %labeling efficiency for PCASL
SupressionEfficiency = 0.75; %Effect of background suppression on labeled spins (0.75 for GE 3D PCASL)
Lambda = 0.9; %blood-brain partition coefficient

%% Step 2: Calculating the CBF images

% CBF = (6000 * Lambda * DM * exp(PLD/BloodT1)) ./ (2 * LabelingEfficiency * SupressionEfficiency * SCF * NumberOfAverages * BloofT1 * M0 * (1-exp(-LabelingDuration/BloodT1)))

DM = spm_vol(ppparams.edeltam);
deltamdat = double(spm_read_vols(DM));

asldim = size(deltamdat);
if numel(asldim)==4, tdim=asldim(4); else tdim=1; end
if tdim==1, SCF = 32; else SCF=1; end %@GE the deltam images are upscalled by a factor 32 in PCASL

M0 = spm_vol(ppparams.tm0scan);
m0dat = double(spm_read_vols(M0));

tcbfdat = zeros(asldim(1),asldim(2),asldim(3),tdim);

CBF=DM;
rmfield(CBF,'pinfo');
for it=1:tdim
    dmdat = reshape(deltamdat(:,:,:,it),[asldim(1),asldim(2),asldim(3)]);

    cbfdat = (6000*Lambda*dmdat*exp(PLD/BloodT1)) ./ (2*LabelingEfficiency*SupressionEfficiency*SCF*NumberOfAverages*BloodT1*m0dat*(1-exp(-LabelingDuration/BloodT1)));
    tcbfdat(:,:,:,it) = cbfdat .* (m0mask>0);
    
    [dmpth,dmname,~] = fileparts(ppparams.edeltam);
    dmparts = split(dmname,'_deltam');
    
    cbffile = fullfile(dmpth,['q' dmparts{1} '_cbf.nii']);

    CBF(it).fname = cbffile;
    CBF(it).descrip = 'CBF';
    CBF(it).n = [it 1];
end

CBF = myspm_write_vol_4d(CBF,tcbfdat);

ppparams.qcbfmap = cbffile;
delfiles{numel(delfiles)+1} = {cbffile};