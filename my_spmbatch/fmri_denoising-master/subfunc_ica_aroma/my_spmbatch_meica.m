function [kappas,rhos,kappa_elbow,rhos_elbow] = my_spmbatch_meica(z_maps,icaTimecourse,ppparams,params)

% This function is based on MEICA implementation in tedana (https://github.com/ME-ICA/tedana)

for ie=1:numel(ppparams.echoes)
    jsondat = jsondecode(fileread(ppparams.func(ppparams.echoes(ie)).jsonfile));
    tes(ie) = jsondat.("EchoTime")*1000;

    Vtemp = spm_vol(fullfile(ppparams.subfuncdir,[ppparams.func(ppparams.echoes(ie)).prefix ppparams.func(ppparams.echoes(ie)).funcfile]));
    if ie==1
        voldim = Vtemp(1).dim;
        nvox = prod(voldim(1:3));
        [ntime,ncomp] = size(icaTimecourse);
        mask = zeros([voldim(1),voldim(2),voldim(3)]);
        mebeta = zeros([nvox,numel(ppparams.echoes),ncomp]);
        mu = zeros([nvox,numel(ppparams.echoes)]);
    end

    mefuncdat = spm_read_vols(Vtemp);

    iemask = my_spmbatch_mask(mefuncdat);
    mask = mask + iemask;

    mefuncdat = reshape(mefuncdat,[nvox,ntime]);

    mu(:,ie) = mean(mefuncdat,2);

    mebeta(:,ie,:) = (pinv(icatb_remove_mean(icaTimecourse)) * icatb_remove_mean(mefuncdat'))';

    clear mefuncdat Vtemp iemask
end

tes = reshape(tes,[numel(tes),1]);

% Model 1 & 2
x1 = mu';
x2 = repmat(tes,[1,nvox]) .* mu';

f_t2_maps = zeros([nvox, ncomp]);
f_s0_maps = zeros([nvox, ncomp]);

for i_comp=1:ncomp
    comp_betas = mebeta(:, :, i_comp)';
    alpha = sum(abs(comp_betas).^2,1);

    % Only analyze good echoes at each voxel
    for j_echo=3:numel(tes)
        mask_idx = find(mask>=j_echo);
        alpha = sum(abs(comp_betas(1:j_echo,:)).^2,1);

        % S0 Model
        % (S,) model coefficient map
        coeffs_s0 = sum(comp_betas(1:j_echo,:) .* x1(1:j_echo,:),1) ./ sum(x1(1:j_echo,:).^2,1);
        pred_s0 = x1(1:j_echo,:) .* repmat(coeffs_s0, [j_echo,1]);
        sse_s0 = (comp_betas(1:j_echo,:) - pred_s0).^2;
        sse_s0 = sum(sse_s0,1);  % (S,) prediction error map
        f_s0 = (alpha - sse_s0) * (j_echo - 1) ./ (sse_s0);
        f_s0(f_s0 > 500) = 500;
        f_s0_maps(mask_idx, i_comp) = f_s0(mask_idx);
        f_s0_maps(isnan(f_s0_maps(:,i_comp)),i_comp)=0;

        clear coeffs_s0 pred_s0 sse_s0 f_s0

        % T2 Model
        coeffs_t2 = sum(comp_betas(1:j_echo,:) .* x2(1:j_echo,:),1) ./ sum(x2(1:j_echo,:).^2,1);
        pred_t2 = x2(1:j_echo,:) .* repmat(coeffs_t2, [j_echo,1]);
        sse_t2 = (comp_betas(1:j_echo,:) - pred_t2).^2;
        sse_t2 = sum(sse_t2,1);
        f_t2 = (alpha - sse_t2) * (j_echo - 1) ./ (sse_t2);
        f_t2(f_t2 > 500) = 500;
        f_t2_maps(mask_idx, i_comp) = f_t2(mask_idx);
        f_t2_maps(isnan(f_t2_maps(:,i_comp)),i_comp)=0;

        clear coeffs_t2 pred_t2 sse_t2 f_t2
    end

    clear comp_betas alpha
end

clear mebeta mu x1 x2

mask_idx = find(mask>=1);

z_maps(isnan(z_maps))=0;
weight_maps = z_maps.^2.0;
kappas = zeros([ncomp,1]); rhos = zeros([ncomp,1]);
for i_comp=1:ncomp
    kappas(i_comp) = mean(f_t2_maps(mask_idx, i_comp).*weight_maps(mask_idx, i_comp));
    rhos(i_comp) = mean(f_s0_maps(mask_idx, i_comp).*weight_maps(mask_idx, i_comp));
end

%kappa and rho elbow calculations
kappa_elbow = getelbow(kappas,ppparams.echoes);
rhos_elbow = getelbow(rhos,ppparams.echoes);

function elbow = getelbow(arr,nechoes)

F01 = spm_invFcdf((1-0.01),1,numel(nechoes)-1);

arr = reshape(arr,[1,numel(arr)]);
nonsig_arr = arr(arr<F01);
nonsig_arr = sort(nonsig_arr,'descend');
n_nns_arr = numel(nonsig_arr);
coords = [0:n_nns_arr-1; nonsig_arr];

p = coords-reshape(coords(:,1),[2,1]);
b=p(:,end);
b_hat=reshape(b./sqrt(sum(b.^2,'all')),[2,1]);
proj_p_b = p - b_hat'*p.*repmat(b_hat,[1,n_nns_arr]);
d = sum(abs(proj_p_b),1);
[~,k_min_ind] = max(d);

elbow = nonsig_arr(k_min_ind);