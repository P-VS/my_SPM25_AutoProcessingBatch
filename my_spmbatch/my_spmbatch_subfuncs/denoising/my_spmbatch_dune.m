function [ppparams,keepfiles,delfiles] = my_spmbatch_dune(ppparams,params,keepfiles,delfiles)

if params.func.isaslbold && contains(ppparams.func(1).funcfile,'_aslbold.nii')
    fname = split(ppparams.func(1).funcfile,'_aslbold.nii');
    ppparams.func(1).funcfile = [fname{1} '_bold.nii'];
    ppparams.perf(1).perffile = [fname{1} '_asl.nii'];
    ppparams.perf(1).prefix = ppparams.func(1).prefix;
end

if contains(ppparams.func(1).funcfile,'_echo-')
    nfname = split(ppparams.func(1).funcfile,'_echo-');
    ppparams.func(1).funcfile = [nfname{1} '_bold.nii'];
end

if contains(ppparams.func(1).funcfile,'_bold')
    nfname = split(ppparams.func(1).funcfile,'_bold');
    ppparams.func(1).funcfile = [nfname{1} '_dune-' params.denoise.DUNE_part '_bold.nii'];
end

ppparams.func = ppparams.func(1);
ppparams.func(1).prefix = ['cd' ppparams.func(1).prefix];

if params.func.isaslbold
    nfname = split(ppparams.perf(1).perffile,'_echo-');
    ppparams.perf(1).perffile = [nfname{1} '_dune-asl_asl.nii'];

    ppparams.perf(1).prefix = ['cd' ppparams.perf(1).prefix];

    if ~exist(fullfile(ppparams.subpath,'perf'),'dir'), mkdir(fullfile(ppparams.subpath,'perf')); end
    ppparams.subperfdir = fullfile(ppparams.subpath,'perf');

    oldfile = fullfile(ppparams.subfuncdir,[ppparams.perf(1).prefix ppparams.perf(1).perffile]);
    newfile = fullfile(ppparams.subperfdir,[ppparams.perf(1).prefix ppparams.perf(1).perffile]);

    if exist(oldfile, "file"), movefile(oldfile,newfile); end
end

ppparams.echoes = 1;
ppparams.meepi = false;

delfiles{numel(delfiles)+1} = {fullfile(ppparams.subfuncdir,[ppparams.func(1).prefix ppparams.func(1).funcfile])};