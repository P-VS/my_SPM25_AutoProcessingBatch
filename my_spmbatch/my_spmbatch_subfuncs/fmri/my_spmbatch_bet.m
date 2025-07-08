function [outfile,delfiles] = my_spmbatch_bet(infolder,infile,ppparams,params,delfiles,keepfiles)

% ---------------------------------------------------------------------
% Brain extraction using MRTOOL
% ---------------------------------------------------------------------
fprintf('Start brain extraction \n')

outfile = ['b' infile];

deltmp = false;

Vin = spm_vol(fullfile(infolder,infile));
if numel(Vin)>0
    indat = spm_read_vols(Vin(1));

    Vout = Vin(1);
    Vout.fname = fullfile(infolder,['tmp' infile]);
    Vout = spm_write_vol(Vout,indat);

    infile = ['tmp' infile];
    
    deltmp = true;
end

MRTool_brain.res_dir = {infolder};
MRTool_brain.t1w = {fullfile(infolder,infile)};

mrt_brain_extraction(MRTool_brain);

if deltmp, delete(fullfile(infolder,infile)); end

snm = split(infile ,'.nii');
tmpfile = [snm{1} '_masked.nii'];

movefile(fullfile(infolder,tmpfile),fullfile(infolder,outfile));

delfiles{numel(delfiles)+1} = {fullfile(infolder,outfile)};

fprintf('Done brain extraction \n')