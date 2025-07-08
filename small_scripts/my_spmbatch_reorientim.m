function my_spmbatch_reorientim(image,dir,contrast)

[pth nm ext] = fileparts(image);
transfile = fullfile(dir,[nm '_reorient.mat']);
if isfile(transfile)
    load(transfile,'M')
    transM = M;
else
    transM = eye(4);
end

Vanat = spm_vol(image);
MM = Vanat.private.mat0;

Vanat = my_reset_orientation(Vanat,transM * MM);

anatdat = spm_read_vols(Vanat);

Vanat.fname = spm_file(image, 'prefix','r');
Vanat.descrip = 'reoriented';
Vanat = spm_create_vol(Vanat);
Vanat = spm_write_vol(Vanat,anatdat);

image = spm_file(image, 'prefix','r');

auto_acpc_reorient(image,contrast);