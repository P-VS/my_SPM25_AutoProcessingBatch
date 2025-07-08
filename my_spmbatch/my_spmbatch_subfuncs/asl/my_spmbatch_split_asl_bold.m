function [ppparams,delfiles,keepfiles] = my_spmbatch_split_asl_bold(params,ppparams,ie,delfiles,keepfiles)

fprintf('Split ASL/BOLD \n')

Vfunc = spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(ie).prefix ppparams.func(ie).funcfile]));
funcdat = spm_read_vols(Vfunc);

mask = my_spmbatch_mask(funcdat);

jsondat = fileread(ppparams.func(ppparams.echoes(1)).jsonfile);
jsondat = jsondecode(jsondat);

tr = jsondat.RepetitionTime;
Ny = 1/(2*tr);

s = size(funcdat);
funcdat = reshape(funcdat(:,:,:,:),[prod(s(1:end-1)),s(end)]);

tmp=find(mask>0);
mfuncdat = funcdat(tmp,:);

mbolddat = my_spm_batch_filt_funcdata(mfuncdat, 0.008, Ny-0.008, tr);

bolddat = zeros(size(funcdat));
bolddat(tmp,:) = mbolddat;
bolddat = reshape(bolddat(:,:),s);

fname = split(ppparams.func(ie).funcfile,'_aslbold.nii');

Vbold = Vfunc;
for iv=1:numel(Vbold)
    Vbold(iv).fname = fullfile(ppparams.subfuncdir,['f' ppparams.func(ie).prefix fname{1} '_aslbold.nii']);
    Vbold(iv).descrip = 'my_spmbatch - split bold';
    Vbold(iv).pinfo = [1,0,0];
    Vbold(iv).n = [iv 1];
end

Vbold = myspm_write_vol_4d(Vbold,bolddat);

clear bolddat Vbold

if ~exist(fullfile(ppparams.subpath,'perf'),'dir'), mkdir(fullfile(ppparams.subpath,'perf')); end
ppparams.subperfdir = fullfile(ppparams.subpath,'perf');

tmpdat = my_spm_batch_filt_funcdata(mfuncdat, 0.008, Inf, tr);
masldat = tmpdat-mbolddat;

asldat = zeros(size(funcdat));
asldat(tmp,:) = masldat;
asldat = reshape(asldat(:,:),s);

clear masldat mfuncdat tmpdat mbolddat

Vasl = Vfunc;
for iv=1:numel(Vasl)
    Vasl(iv).fname = fullfile(ppparams.subperfdir,['f' ppparams.func(ie).prefix fname{1} '_label.nii']);
    Vasl(iv).descrip = 'my_spmbatch - label asl';
    Vasl(iv).pinfo = [1,0,0];
    Vasl(iv).n = [iv 1];
end

Vasl = myspm_write_vol_4d(Vasl,asldat);

delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir,['f' ppparams.func(ie).prefix fname{1} '_aslbold.nii'])};

ppparams.func(ie).funcfile = [fname{1} '_aslbold.nii'];
ppparams.func(ie).prefix = ['f' ppparams.func(ie).prefix];

clear funcdat asldat Vfunc Vasl