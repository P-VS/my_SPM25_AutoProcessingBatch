function [delfiles,keepfiles] = my_spmbatch_cat12vbm(ppparams,params,delfiles,keepfiles)

fname = split(ppparams.subanat,'.nii');

delfiles{numel(delfiles)+1} = {fullfile(ppparams.subanatdir,['cat_' fname{1} '.xml'])};
delfiles{numel(delfiles)+1} = {fullfile(ppparams.subanatdir,['catlog_' fname{1} '.txt'])};

catestwrite.data = {fullfile(ppparams.subanatdir,ppparams.subanat)};
catestwrite.data_wmh = {''};
catestwrite.nproc = 4;
catestwrite.useprior = '';
catestwrite.opts.tpm = {fullfile(spm('Dir'),'tpm','TPM.nii')};
catestwrite.opts.affreg = 'mni';
catestwrite.opts.biasacc = 0.5;
catestwrite.extopts.restypes.optimal = [1 0.3];

catestwrite.extopts.setCOM = 1;

catestwrite.extopts.APP = 1070;
catestwrite.extopts.affmod = 0;
catestwrite.extopts.spm_kamap = 0;
catestwrite.extopts.LASstr = 0.5;
catestwrite.extopts.LASmyostr = 0;
catestwrite.extopts.gcutstr = 2;
catestwrite.extopts.WMHC = 2;
catestwrite.extopts.registration.shooting.shootingtpm = {fullfile(spm('Dir'),'toolbox','cat12','templates_MNI152NLin2009cAsym','Template_0_GS.nii')};
catestwrite.extopts.registration.shooting.regstr = 0.5;
catestwrite.extopts.vox = params.vbm.normvox;
catestwrite.extopts.bb = 12;
catestwrite.extopts.SRP = 22;
catestwrite.extopts.ignoreErrors = 1;

catestwrite.output.BIDS.BIDSno = 0;

if params.vbm.do_surface
    catestwrite.output.surface = 1;
    catestwrite.output.surf_measures = 1;

    keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subanatdir,['lh.central.' fname{1} '.gii'])};
    keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subanatdir,['rh.central.' fname{1} '.gii'])};
    keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subanatdir,['lh.pial.' fname{1} '.gii'])};
    keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subanatdir,['rh.pial.' fname{1} '.gii'])};
    keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subanatdir,['lh.sphere.' fname{1} '.gii'])};
    keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subanatdir,['rh.sphere.' fname{1} '.gii'])};
    keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subanatdir,['lh.sphere.reg.' fname{1} '.gii'])};
    keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subanatdir,['rh.sphere.reg.' fname{1} '.gii'])};
    keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subanatdir,['lh.white.' fname{1} '.gii'])};
    keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subanatdir,['rh.white.' fname{1} '.gii'])};
    keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subanatdir,['lh.pbt.' fname{1}])};
    keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subanatdir,['rh.pbt.' fname{1}])};
    keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subanatdir,['lh.thickness.' fname{1}])};
    keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subanatdir,['rh.thickness.' fname{1}])};
else
    catestwrite.output.surface = 0;
    catestwrite.output.surf_measures = 0;
end

if params.vbm.do_roi_atlas
    catestwrite.output.ROImenu.atlases.neuromorphometrics = 1;
    catestwrite.output.ROImenu.atlases.lpba40 = 1;
    catestwrite.output.ROImenu.atlases.cobra = 1;
    catestwrite.output.ROImenu.atlases.hammers = 0;
    catestwrite.output.ROImenu.atlases.thalamus = 1;
    catestwrite.output.ROImenu.atlases.thalamic_nuclei = 1;
    catestwrite.output.ROImenu.atlases.suit = 1;
    catestwrite.output.ROImenu.atlases.ibsr = 0;
    catestwrite.output.ROImenu.atlases.ownatlas = {''};

    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subanatdir,['catROI_' fname{1} '.xml'])};
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subanatdir,['catROIs_' fname{1} '.xml'])};

    keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subanatdir,['catROI_' fname{1} '.mat'])};
    keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subanatdir,['catROIs_' fname{1} '.mat'])};
else
    catestwrite.output.ROImenu.noROI = struct([]);
end

catestwrite.output.GM.mod = 1;
catestwrite.output.GM.dartel = 0;
catestwrite.output.WM.mod = 1;
catestwrite.output.WM.dartel = 0;
catestwrite.output.CSF.mod = 1;
catestwrite.output.CSF.dartel = 0;
catestwrite.output.CSF.warped = 0;

keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subanatdir,['wm' fname{1} '.nii'])};
keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subanatdir,['mwp1' fname{1} '.nii'])};
keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subanatdir,['mwp2' fname{1} '.nii'])};
keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subanatdir,['mwp3' fname{1} '.nii'])};

if params.vbm.do_normalization
    catestwrite.output.GM.native = 0;
    catestwrite.output.GM.mod = 1;
    catestwrite.output.WM.native = 0;
    catestwrite.output.WM.mod = 1;
    catestwrite.output.CSF.native = 0;
    catestwrite.output.CSF.mod = 1;
    catestwrite.output.labelnative = 0;
else
    catestwrite.output.GM.native = 1;
    catestwrite.output.GM.mod = 0;
    catestwrite.output.WM.native = 1;
    catestwrite.output.WM.mod = 0;
    catestwrite.output.CSF.native = 1;
    catestwrite.output.CSF.mod = 0;
    catestwrite.output.labelnative = 1;

    keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subanatdir,['p0' fname{1} '.nii'])};
    keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subanatdir,['p1' fname{1} '.nii'])};
    keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subanatdir,['p2' fname{1} '.nii'])};
    keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subanatdir,['p3' fname{1} '.nii'])};
end

catestwrite.output.ct.native = 0;
catestwrite.output.ct.warped = 0;
catestwrite.output.ct.dartel = 0;
catestwrite.output.pp.native = 0;
catestwrite.output.pp.warped = 0;
catestwrite.output.pp.dartel = 0;
catestwrite.output.WMH.native = 0;
catestwrite.output.WMH.warped = 0;
catestwrite.output.WMH.mod = 0;
catestwrite.output.WMH.dartel = 0;
catestwrite.output.SL.native = 0;
catestwrite.output.SL.warped = 0;
catestwrite.output.SL.mod = 0;
catestwrite.output.SL.dartel = 0;
catestwrite.output.TPMC.native = 0;
catestwrite.output.TPMC.warped = 0;
catestwrite.output.TPMC.mod = 0;
catestwrite.output.TPMC.dartel = 0;
catestwrite.output.atlas.native = 0;
catestwrite.output.label.native = 1;
catestwrite.output.label.warped = 0;
catestwrite.output.label.dartel = 0;
catestwrite.output.bias.warped = 1;
catestwrite.output.las.native = 0;
catestwrite.output.las.warped = 0;
catestwrite.output.las.dartel = 0;
catestwrite.output.jacobianwarped = 0;
catestwrite.output.warps = [0 0];
catestwrite.output.rmat = 0;

cat_run(catestwrite);

delfiles{numel(delfiles)+1} = {fullfile(ppparams.subanatdir,'mri')};
delfiles{numel(delfiles)+1} = {fullfile(ppparams.subanatdir,['catlog_' fname{1} '.txt'])};

keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subanatdir,['cat_' fname{1} '.mat'])};
keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subanatdir,['catreport_' fname{1} '.pdf'])};
keepfiles{numel(keepfiles)+1} = {fullfile(ppparams.subanatdir,['catreportj_' fname{1} '.jpg'])};