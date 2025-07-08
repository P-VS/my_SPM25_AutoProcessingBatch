function [delfiles,keepfiles] = my_spmbatch_rfasl(ppparams,params,delfiles,keepfiles)

ppparams.subdmjson = ppparams.funcjsonfile;

jsondat = fileread(ppparams.funcjsonfile);
jsondat = jsondecode(jsondat);
tr = jsondat.RepetitionTime;
te1 = 1000.0*jsondat.EchoTime;

%% T2* corrected echo combbination
nechoes = numel(ppparams.fasl.rasl);

for ie=1:nechoes
    ppparams.subfuncstring = [ppparams.substring '_task-' ppparams.task '_bold_e'];
    jsonfile = fullfile(ppparams.subfmridir,[ppparams.subfuncstring num2str(ie) '.json']);
    
    jsondat = fileread(jsonfile);
    jsondat = jsondecode(jsondat);
    te(ie) = 1000.0*jsondat.EchoTime;

    if ie==1
        voldim = size(ppparams.fasl.rasl{ie}.dat);
        teasldat = zeros(voldim(1),voldim(2),voldim(3),voldim(4),nechoes);
        m0idat = zeros(voldim(1),voldim(2),voldim(3),nechoes);
    end

    teasldat(:,:,:,:,ie) = ppparams.fasl.rasl{ie}.dat;

    M0I = spm_vol(ppparams.fasl.rasl{ie}.subm0scan);
    m0idat(:,:,:,ie) = spm_read_vols(M0I);
end

mask = my_spmbatch_mask(teasldat(:,:,:,:,1));
mask_ind = find(mask>0);

%based on https://github.com/jsheunis/fMRwhy/tree/master
% Create "design matrix" X
X = horzcat(ones(nechoes,1), -te(:));

t2star = zeros(voldim(1)*voldim(2)*voldim(3),1);

Y=[];
for ne=1:nechoes
    tempteasldat = reshape(m0idat(:,:,:,ne),[voldim(1)*voldim(2)*voldim(3),1]);
    Y=[Y;reshape(tempteasldat(mask_ind,1),[1,numel(mask_ind)])];
end
Y = max(Y, 1e-11);

% Estimate "beta matrix" by solving set of linear equations
beta_hat = pinv(X) * log(Y);
 % Calculate S0 and T2star from beta estimation
T2star_fit = beta_hat(2, :); %is R2*

T2star_thresh_min = 1/500; % arbitrarily chosen, same as tedana
I_T2star_min = (T2star_fit < T2star_thresh_min); % vector of voxels where T2star value is negative
T2star_fit(I_T2star_min) = 0; % if values inside mask are zero or negative, set them to threshold_min value
T2star_thresh_max = 1/10; % arbitrarily chosen, same as tedana
I_T2star_max = (T2star_fit > T2star_thresh_max); % vector of voxels where T2star value is negative
T2star_fit(I_T2star_max) = 0; % if values inside mask are zero or negative, set them to threshold_min value

t2star(mask_ind) = T2star_fit;
    
asldat = zeros(voldim(1),voldim(2),voldim(3),voldim(4));

weights = zeros(voldim(1)*voldim(2)*voldim(3),nechoes);

for ne=1:nechoes
    weights(:,ne) = repmat((te(ne)),voldim(1)*voldim(2)*voldim(3),1) .* t2star(:,1); %-te(1)
    weights(:,ne) = exp(weights(:,ne));
end

weights = reshape(weights,[voldim(1),voldim(2),voldim(3),nechoes]);

sum_weights = sum(weights,4);
weights_mask = find(sum_weights>0);

m0dat = zeros(voldim(1),voldim(2),voldim(3));

for ne=1:nechoes
    mwidat = m0idat(:,:,:,ne);
    mwidat = mwidat .* weights(:,:,:,ne);
    mwidat(weights_mask) = mwidat(weights_mask) ./ nechoes;%sum_weights(weights_mask);

    m0dat = m0dat+mwidat; 
end

[fpath,fname,~] = fileparts(ppparams.fasl.rasl{1}.subm0scan);
nfname = split(fname,'bold_e');

M0 = M0I;
M0.fname = fullfile(fpath,['c' nfname{1} 'bold.nii']);
M0.descrip = 'my_spmbatch - m0scan';
M0.n = [1 1];
M0 = spm_create_vol(M0);
M0 = spm_write_vol(M0,m0dat);

ppparams.subm0scan = M0.fname;

delfiles{numel(delfiles)+1} = {ppparams.subm0scan};

for ti=1:voldim(4)
    for ne=1:nechoes
        asltidat = teasldat(:,:,:,ti,ne);
        asltidat = asltidat .* weights(:,:,:,ne);
        asltidat(weights_mask) = asltidat(weights_mask) ./ sum_weights(weights_mask);

        asldat(:,:,:,ti) = asldat(:,:,:,ti)+asltidat; 
    end 
end

Vfunc = spm_vol(ppparams.fasl.rasl{1}.file);

[fpath,fname,~] = fileparts(Vfunc(1).fname);
nfname = split(fname,'bold_e');

Vout = Vfunc;

for j=1:numel(Vout)
    Vout(j).fname = fullfile(fpath,['c' nfname{1} 'bold.nii']);
    Vout(j).descrip = 'my_spmbatch - combine echoes';
    Vout(j).pinfo = [1,0,0];
    Vout(j).dt = [spm_type('float32'),spm_platform('bigend')];
    Vout(j).n = [j 1];
end

Vout = myspm_write_vol_4d(Vout,asldat);

ppparams.aslfile = Vout.fname;

delfiles{numel(delfiles)+1} = {ppparams.aslfile};

%% Subtract labeled and control data into deltam
fprintf('Start subtraction control-label \n')

vol1 = asldat(:,:,:,1);%ppparams.fasl.rasldat(:,:,:,1);
vol2 = asldat(:,:,:,2);%ppparams.fasl.rasldat(:,:,:,2);

tmask = find(vol1>0);
dvol = mean(vol1(tmask)-vol2(tmask),'all');

if dvol>0, lorder=1; else lorder=0; end

niihdr = niftiinfo(ppparams.fasl.rasl{1}.file);
fov = niihdr.ImageSize(1) * niihdr.PixelDimensions(1);

ppparams.subdeltamdat = aslsub(asldat,'fstart',1,'fend','auto','order',lorder,'tr',tr,'fov',fov,'sur',1,'rel',0,'scaleoutput',1);
ppparams.subdeltamdat(ppparams.subdeltamdat<0) = 0;

mean_sub = mean(ppparams.subdeltamdat,4);

Vasl = spm_vol(ppparams.aslfile);

mVasl = Vasl(1);
mVasl.fname = spm_file(ppparams.aslfile, 'prefix','sub');
mVasl.descrip = 'asl - subtracted';
mVasl.n = [1 1];
mVasl = rmfield(mVasl,'pinfo');
mVasl = spm_create_vol(mVasl);
mVasl = spm_write_vol(mVasl,mean_sub);

ppparams.subdeltam = mVasl.fname;
delfiles{numel(delfiles)+1} = {mVasl.fname};

%% Corregistration M0map to deltam
fprintf('Start corregistration \n')

estwrite.ref(1) = {ppparams.subdeltam};
estwrite.source(1) = {ppparams.subm0scan};
estwrite.other = {ppparams.subm0scan};
estwrite.eoptions.cost_fun = 'nmi';
estwrite.eoptions.sep = [4 2];
estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
estwrite.eoptions.fwhm = [7 7];
estwrite.roptions.interp = 4;
estwrite.roptions.wrap = [0 0 0];
estwrite.roptions.mask = 0;
estwrite.roptions.prefix = 'r';

out_coreg = spm_run_coreg(estwrite);

ppparams.subm0scan = spm_file(ppparams.subm0scan, 'prefix','r');
delfiles{numel(delfiles)+1} = {ppparams.subm0scan};

%% T1 correction of M0
fprintf('Start M0 preprocessing \n')
[ppparams,delfiles,keepfiles] = my_spmbatch_asl_M0correction(ppparams,params,delfiles,keepfiles);

%% Quantification of CBF
fprintf('Start CBF quantification \n')
[ppparams,delfiles,keepfiles] = my_spmbatch_asl_cbfquantification(ppparams,delfiles,keepfiles);