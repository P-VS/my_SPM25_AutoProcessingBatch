function [ppparams,delfiles,keepfiles] = my_spmbatch_fasl_cbfmapping(ppparams,params,delfiles,keepfiles)

%used CBF model ass in Alsop et al. Recommended Implementation of Arterial Spin-Labeled Perfusion MRI for Clinical Applications: 
% %A Consensus of the ISMRM Perfusion Study Group and the European
% Consortium for ASL in Dementia. Magnetic Resonance in Medicine 73:102–116 (2015)
%(https://onlinelibrary.wiley.com/doi/epdf/10.1002/mrm.25197?src=getftr)

%% Coorecting PLD for slice time differences

Vdeltam=spm_vol(fullfile(ppparams.subperfdir,[ppparams.perf(1).deltamprefix ppparams.perf(1).deltamfile]));

nfname = split(ppparams.perf(1).deltamfile,'_deltam');

tdim = numel(Vdeltam);
voldim = Vdeltam(1).dim;

jsondat = fileread(ppparams.func(1).jsonfile);
jsondat = jsondecode(jsondat);
tr = jsondat.RepetitionTime;
te = 1000*jsondat.EchoTime;
if isfield(jsondat,'LabelingDuration'), LD = jsondat.LabelingDuration; else LD = params.asl.LabelingDuration; end
if isfield(jsondat,'PostLabelDelay'), PLD = jsondat.PostLabelDelay; else PLD = params.asl.PostLabelDelay; end

if isfield(jsondat,'SliceTiming'), SliceTimes = jsondat.SliceTiming; else SliceTimes = []; end
if ~(numel(SliceTimes)==voldim(3)), SliceTimes = []; end

isfield(jsondat,'MultibandAccellerationFactor')
if isempty(SliceTimes)
    if isfield(jsondat,'MultibandAccelerationFactor') || isfield(jsondat,'MultibandAccellerationFactor')
        if isfield(jsondat,'MultibandAccelerationFactor'), hbf = jsondat.MultibandAccelerationFactor; 
        else hbf = jsondat.MultibandAccellerationFactor; end
        nslex = ceil(voldim(3)/hbf);
        isl = zeros([1,nslex]);
        isl(1:2:nslex)=[0:1:(nslex-1)/2];
        isl(2:2:nslex)=[ceil(nslex/2):1:nslex-1];
        isl=repmat(isl,[1,hbf]);
        isl = isl(1:voldim(3));
    else 
        isl = [1:2:voldim(3) 2:2:voldim(3)];
        isl = isl-1;
        nslex = voldim(3);
    end

    TA = tr-LD-PLD;
    SliceTimes = isl*TA/nslex;
else
    if max(SliceTimes)>(tr-LD-PLD)
        TA = tr-LD-PLD;
    
        SliceTimes = SliceTimes * TA/tr;
    end
end

vol_PLD = zeros(voldim(1),voldim(2),voldim(3));
for is=1:voldim(3)
    vol_PLD(:,:,is) = PLD+SliceTimes(is);
end

reshape(vol_PLD,[voldim(1)*voldim(2),voldim(3)]);

%% Correct M0 for T1 effects
% The T1 values used, are the averaged T1 values reported in the review of 
% Bojorquez et al. 2017. What are normal relaxation times of tissues at 3 T? Magnetic Resonance Imaging 35:69-80
% (https://mri-q.com/uploads/3/4/5/7/34572113/normal_relaxation_times_at_3t.pdf)

T1a = 1.650; %longitudinal relaxation time of arterial blood
lambda = 0.9; %blood-brain partition coefficient for gray matter
alpha = 0.85; %laeling efficiency

Vasl=spm_vol(fullfile(ppparams.subperfdir,[ppparams.perf(1).aslprefix ppparams.perf(1).aslfile]));
fasldata = spm_read_vols(Vasl);

conidx = 2:2:numel(Vasl);

m0vol = mean(fasldata(:,:,:,conidx),4);

mask = my_spmbatch_mask(m0vol);

GM = spm_vol(fullfile(ppparams.subperfdir,ppparams.perf(1).c1m0scanfile));
WM = spm_vol(fullfile(ppparams.subperfdir,ppparams.perf(1).c2m0scanfile));
CSF = spm_vol(fullfile(ppparams.subperfdir,ppparams.perf(1).c3m0scanfile));

gmim = spm_read_vols(GM);
wmim = spm_read_vols(WM);
csfim = spm_read_vols(CSF);

T1gm = 1.459;
T1wm = 0.974;
T1csf = 4.190;

T1dat = (T1gm * gmim + T1wm * wmim + T1csf * csfim);

corr_T1 = zeros(voldim);
corr_T1(T1dat>0) = 1 ./ (1-exp(-tr./T1dat(T1dat>0)));
m0vol = m0vol .* corr_T1;

clear fasldata gmim wmim csfim T1dat corr_T1 GM WM CSF Vm0 Vasl

%% CBF calculations series
cm0vol = 2*alpha*m0vol*T1a.*(exp(-vol_PLD/T1a)-exp(-(LD+vol_PLD)/T1a));
cm0vol = reshape(cm0vol,[voldim(1)*voldim(2)*voldim(3),1]);

Vout = Vdeltam;
rmfield(Vout,'pinfo');

nvols = params.loadmaxvols;
for ti=1:nvols:tdim
    if ti+nvols>tdim, nvols=tdim-ti+1; end

    fprintf(['CBF vols: ' num2str(ti) '-' num2str(ti+nvols-1) '\n'])
    
    cbfdata = spm_read_vols(Vdeltam(ti:ti+nvols-1));
    cbfdata = reshape(cbfdata,[voldim(1)*voldim(2)*voldim(3),nvols]);

    cbfdata = lambda*6000*cbfdata;
    cbfdata(cm0vol>0,:) = cbfdata(cm0vol>0,:)./repmat(cm0vol(cm0vol>0),[1,nvols]);
    cbfdata = reshape(cbfdata,[voldim(1),voldim(2),voldim(3),nvols]);

    cbfdata = cbfdata.*repmat(mask,[1,1,1,nvols]);
    cbfdata(cbfdata<-40) = 0;
    cbfdata(cbfdata>150) = 0;

    for iv=1:nvols
        Vout(iv).fname = fullfile(ppparams.subperfdir,[ppparams.perf(1).deltamprefix nfname{1} '_cbf.nii']);
        Vout(iv).descrip = 'my_spmbatch - cbf';
        Vout(iv).dt = [spm_type('float32'),spm_platform('bigend')];
        Vout(iv).n = [ti+iv-1 1];
        Vout(iv) = spm_write_vol(Vout(iv),cbfdata(:,:,:,iv));
    end

    clear cbfdata
end

clear  Vout Vdeltam

%% CBF calculations mean deltaM
Vmeandm=spm_vol(fullfile(ppparams.subperfdir,ppparams.perf(1).meandeltam));
meand_nfname = split(ppparams.perf(1).meandeltam,'_deltam');

cbfdata = spm_read_vols(Vmeandm);
cbfdata = reshape(cbfdata,[voldim(1)*voldim(2)*voldim(3),1]);

cbfdata = lambda*6000*cbfdata;
cbfdata(cm0vol>0) = cbfdata(cm0vol>0)./cm0vol(cm0vol>0);
cbfdata = reshape(cbfdata,[voldim(1),voldim(2),voldim(3)]);

cbfdata = cbfdata.*mask;
cbfdata(cbfdata<-40) = 0;
cbfdata(cbfdata>150) = 0;

Vout = Vmeandm;
Vout.fname = fullfile(ppparams.subperfdir,[meand_nfname{1} '_cbf.nii']);
Vout.descrip = 'my_spmbatch - cbf';
Vout.dt = [spm_type('float32'),spm_platform('bigend')];
Vout.n = [1 1];
Vout = spm_write_vol(Vout,cbfdata);

clear cbfdata Vout Vmean0 m0vol cm0vol vol_PLD mask

ppparams.perf(1).meancbf = [meand_nfname{1} '_cbf.nii'];
ppparams.perf(1).cbfprefix = ppparams.perf(1).deltamprefix;
ppparams.perf(1).cbffile = [nfname{1} '_cbf.nii'];

delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,[ppparams.perf(1).cbfprefix ppparams.perf(1).cbffile])};    
delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,ppparams.perf(1).meancbf)};   