function [ppparams,delfiles,keepfiles] = my_spmbatch_fieldmap(ne,ppparams,delfiles,keepfiles)

mbstep = 1;

%%Load fieldmap data per echo and coregister to func
 
e1jsondat = fileread(ppparams.fmap.json1);
e1jsondat = jsondecode(e1jsondat);
te1 = e1jsondat.EchoTime*1000;

Ve1amp = spm_vol(ppparams.fmap.amim1);
Ve1ph  = spm_vol(ppparams.fmap.phim1);

if ppparams.reorient
    nfname = split(Ve1amp.fname,'.nii');

    transfile = [nfname(1) '_reorient.mat'];
    if isfile(transfile)
        load(transfile,'M')
        transM = M;
    else
        transM = my_spmbatch_vol_set_com(Ve1amp.fname);
        transM(1:3,4) = -transM(1:3,4);
    end

    MM = Ve1amp.private.mat0;

    Ve1amp = my_reset_orientation(Ve1amp,transM * MM);

    ve1amdat = spm_read_vols(Ve1amp);
    
    Ve1amp.fname = spm_file(Ve1amp.fname, 'prefix','e');;
    Ve1amp.descrip = 'my_spmbatch - reoriented';
    Ve1amp.n = [1 1];
    Ve1amp = spm_create_vol(Ve1amp);
    Ve1amp = spm_write_vol(Ve1amp,ve1amdat);

    clear ve1amdat
    
    Ve1ph = my_reset_orientation(Ve1ph,transM * MM);

    ve1phdat = spm_read_vols(Ve1ph);
    
    Ve1ph.fname = spm_file(Ve1ph.fname, 'prefix','e');;
    Ve1ph.descrip = 'my_spmbatch - reoriented';
    Ve1ph.n = [1 1];
    Ve1ph = spm_create_vol(Ve1ph);
    Ve1ph = spm_write_vol(Ve1ph,ve1phdat);    

    clear ve1phdat
end

fme1step = mbstep;

reffunc = fullfile(ppparams.subfuncdir,['e' ppparams.func(ne).funcfile ',1']);

matlabbatch{mbstep}.spm.spatial.coreg.estwrite.ref(1) = {reffunc};
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.source(1) = {Ve1amp.fname};
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.other(1) = {Ve1ph.fname};
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

mbstep = mbstep+1;

amp1file = spm_file(Ve1amp.fname, 'prefix','r');
ph1file = spm_file(Ve1ph.fname, 'prefix','r');

delfiles{numel(delfiles)+1} = {amp1file};
delfiles{numel(delfiles)+1} = {ph1file};

e2jsondat = fileread(ppparams.fmap.json2);
e2jsondat = jsondecode(e2jsondat);
te2 = e2jsondat.EchoTime*1000;

Ve2amp = spm_vol(ppparams.fmap.amim2);
Ve2ph  = spm_vol(ppparams.fmap.phim2);

if ppparams.reorient
    nfname = split(Ve2amp.fname,'.nii');

    transfile = [nfname(1) '_reorient.mat'];
    if isfile(transfile)
        load(transfile,'M')
        transM = M;
    else
        transM = my_spmbatch_vol_set_com(Ve2amp.fname);
        transM(1:3,4) = -transM(1:3,4);
    end

    MM = Ve2amp.private.mat0;

    Ve2amp = my_reset_orientation(Ve2amp,transM * MM);

    ve2amdat = spm_read_vols(Ve2amp);
    
    Ve2amp.fname = spm_file(Ve2amp.fname, 'prefix','e');;
    Ve2amp.descrip = 'my_spmbatch - reoriented';
    Ve2amp.n = [1 1];
    Ve2amp = spm_create_vol(Ve2amp);
    Ve2amp = spm_write_vol(Ve2amp,ve2amdat);

    clear ve2amdat
    
    Ve2ph = my_reset_orientation(Ve2ph,transM * MM);

    ve2phdat = spm_read_vols(Ve2ph);
    
    Ve2ph.fname = spm_file(Ve2ph.fname, 'prefix','e');
    Ve2ph.descrip = 'my_spmbatch - reoriented';
    Ve2ph.n = [1 1];
    Ve2ph = spm_create_vol(Ve2ph);
    Ve2ph = spm_write_vol(Ve2ph,ve2phdat); 

    clear ve2phdat
end

fme2step = mbstep;

matlabbatch{mbstep}.spm.spatial.coreg.estwrite.ref(1) = {reffunc};
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.source(1) = {Ve2amp.fname};
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.other(1) = {Ve2ph.fname};
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{mbstep}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

mbstep = mbstep+1;

amp2file = spm_file(Ve2amp.fname, 'prefix','r');
ph2file = spm_file(Ve2ph.fname, 'prefix','r');

delfiles{numel(delfiles)+1} = {amp2file};
delfiles{numel(delfiles)+1} = {ph2file};

%% Fieldmap

fmstep = mbstep;
  
jsondat = fileread(fullfile(ppparams.subfuncdir,ppparams.func(1).jsonfile));
jsondat = jsondecode(jsondat);
pedir = jsondat.PhaseEncodingDirection;

if contains(pedir,'-')
    blipdim = -1;
else
    blipdim = 1;
end

try
    trt = jsondat.TotalReadoutTime*1000;
catch
    trt = jsondat.AcquisitionMatrixPE*jsondat.EffectiveEchoSpacing*1000;
end

matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.shortphase(1) = {amp1file};
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.shortmag(1) = {ph1file};
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.longphase(1) = {amp2file};
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.data.phasemag.longmag(1) = {ph2file};
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = [te1 te2];
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 1;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = blipdim;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = trt;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 0;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm = 0;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = 'Mark3D';
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm = 10;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad = 0;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws = 1;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = {fullfile(spm('Dir'),'toolbox','FieldMap','T1.nii')};
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm = 5;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode = 2;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate = 4;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh = 0.5;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg = 0.02;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.session.epi(1) = {reffunc};
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'session';
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 0;
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.anat = '';
matlabbatch{mbstep}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 0;

mbstep = mbstep+1;

ppparams.func(ne).vdm_file = spm_file(amp1file, 'prefix','vdm5_sc');

delfiles{numel(delfiles)+1} = {spm_file(ph1file, 'prefix','m')};
delfiles{numel(delfiles)+1} = {spm_file(ph1file, 'prefix','bmask')};
delfiles{numel(delfiles)+1} = {spm_file(amp1file, 'prefix','sc')};
delfiles{numel(delfiles)+1} = {spm_file(amp2file, 'prefix','sc')};
delfiles{numel(delfiles)+1} = {spm_file(amp1file, 'prefix','fpm_sc')};
delfiles{numel(delfiles)+1} = {spm_file(amp1file, 'prefix','vdm5_sc')};

%%Run matlabbatch
if exist("matlabbatch",'var')
    spm_jobman('run', matlabbatch);
end