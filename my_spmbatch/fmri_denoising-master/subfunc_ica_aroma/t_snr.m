clear all; close all; clc;

mask_file = 'D:\fmri_data\sub-01\ses-001\masks\matlab\func_mask_sub1.nii';
func_file = 'D:\fmri_data\sub-01\ses-001\preproc_func\wuaesub-01_task-SE-EmoFaces_bold.nii';

denoised_func_mel = 'D:\fmri_data\sub-01\ses-001\preproc_func_melodic_denoised\wuaesub-01_task-SE-EmoFaces_bold_regfilt.nii';
denoised_func_gift = 'C:\Users\Sara\Downloads\preproc_func_ica_aroma_denoised\R_wuaesub-01_task-SE-EmoFaces_bold.nii';

% Load mask
Vmask= spm_vol(mask_file);
mask_data = mat2gray(spm_read_vols(Vmask));
mask_data_1d = reshape(mask_data,[],1);

% Preprocessed func before denoising:
Vfunc= spm_vol(func_file);
func_data = spm_read_vols(Vfunc);
% Apply mask
func_data = func_data .* mask_data;
% func_data_2d = reshape(func_data, size(func_data, 1)*size(func_data, 2)*size(func_data, 3), size(func_data, 4));
% Compute the mean over all time points
mean_func = mean(func_data, 4);

% After denoising func GIFT:
Vfunc_den = spm_vol(denoised_func_gift);
denoised_gift = spm_read_vols(Vfunc_den);
% Apply mask
denoised_gift = denoised_gift .* mask_data;
% Compute the mean over all time point
mean_denoised_gift = mean(denoised_gift, 4);

% After denoising func MELODIC:
Vfunc_den = spm_vol(denoised_func_mel);
denoised_mel = spm_read_vols(Vfunc_den);
% Apply mask
denoised_mel = denoised_mel .* mask_data;
% Compute the mean over all time point
mean_denoised_mel = mean(denoised_mel, 4);

% %% Compute tSNR maps

% Pre denoised
std_pre_den_3d = std(func_data, 0, 4);
t_snr_pre_den = mean_func ./ std_pre_den_3d;

% Post denoised - GIFT
std_post_den_gift_3d = std(denoised_gift, 0, 4);
t_snr_post_den_gift = mean_denoised_gift ./ std_post_den_gift_3d;

% Post denoised - MELODIC
std_post_den_mel_3d = std(denoised_mel, 0, 4);
t_snr_post_den_mel = mean_denoised_mel ./ std_post_den_mel_3d;

%% Paired t-test

% Reshape 3D tSNR maps to 2D
t_snr_pre_den_1d = reshape(t_snr_pre_den(~isnan(t_snr_pre_den)), [], 1);
t_snr_post_den_melodic_1d = reshape(t_snr_post_den_mel(~isnan(t_snr_post_den_mel)), [], 1);
t_snr_post_den_gift_1d = reshape(t_snr_post_den_gift(~isnan(t_snr_post_den_gift)), [], 1);

% t_snr_pre_den_1d = reshape(t_snr_pre_den, [], 1);
% t_snr_post_den_melodic_1d = reshape(t_snr_post_den_mel, [], 1);
% t_snr_post_den_gift_1d = reshape(t_snr_post_den_gift, [], 1);

[~, p, ci, stats] = ttest(t_snr_pre_den_1d, t_snr_post_den_gift_1d);
disp("Paired t-test Pre vs Post(GIFT):");
disp("p-value: " + num2str(p));
disp("t-statistic: " + num2str(stats.tstat));
disp(" ");

[~, p, ci, stats] = ttest(t_snr_pre_den_1d, t_snr_post_den_melodic_1d);
disp("Paired t-test Pre vs Post(MELODIC):");
disp("p-value: " + num2str(p));
disp("t-statistic: " + num2str(stats.tstat));
disp(" ");

[~, p, ci, stats] = ttest(t_snr_post_den_melodic_1d, t_snr_post_den_gift_1d);
disp("Paired t-test Post(GIFT) vs POST (MELODIC):");
disp("p-value: " + num2str(p));
disp("t-statistic: " + num2str(stats.tstat));
disp(" ");

[~, p, ci, stats] = ttest(t_snr_post_den_melodic_1d - t_snr_post_den_gift_1d);
disp("Paired t-test Post(GIFT) vs POST (MELODIC):");
disp("p-value: " + num2str(p));
disp("t-statistic: " + num2str(stats.tstat));
disp(" ");

%% Mean and median tsnr

% Exclude NaN from tsnr (background pixels)

% Pre denoising
disp("Pre denoising");
disp("Mean tsnr: " + mean(t_snr_pre_den_1d));
disp("Median tsnr: " + median(t_snr_pre_den_1d));
disp("Min tsnr: " + min(t_snr_pre_den_1d));
disp("Max tsnr: " + max(t_snr_pre_den_1d));
disp("Std: " + std(t_snr_pre_den_1d));
disp(" ");

% Post denoising - GIFT
disp("Post denoising - GIFT");
disp("Mean tsnr: " + mean(t_snr_post_den_gift_1d));
disp("Median tsnr: " + median(t_snr_post_den_gift_1d));
disp("Min tsnr: " + min(t_snr_post_den_gift_1d));
disp("Max tsnr: " + max(t_snr_post_den_gift_1d));
disp("Std: " + std(t_snr_post_den_gift_1d));
disp(" ");

% Post denoising - GIFT
disp("Post denoising - MELODIC");
disp("Mean tsnr: " + mean(t_snr_post_den_melodic_1d));
disp("Median tsnr: " + median(t_snr_post_den_melodic_1d));
disp("Min tsnr: " + min(t_snr_post_den_melodic_1d));
disp("Max tsnr: " + max(t_snr_post_den_melodic_1d));
disp("Std: " + std(t_snr_post_den_melodic_1d));
disp(" ");

% Mean difference Post denoising - GIFT
disp("Difference GIFT - MELODIC");
diff = t_snr_post_den_melodic_1d - t_snr_post_den_gift_1d;
% diff_3d = reshape(diff, size(mask_data,1), size(mask_data,2), size(mask_data,3));
disp("Mean difference tsnr: " + mean(diff));
disp("Median tsnr: " + median(diff));
disp("Min tsnr: " + min(diff));
disp("Max tsnr: " + max(diff));
disp("Std: " + std(diff));

%% Plot tSNR

% Pre denoising
figure;
imagesc(t_snr_pre_den(:,:,80)); %slice 80
colorbar;
colormap(hot); %hot % autum %jet
title('tSNR Map (Pre denoising)');

% Post denoising - GIFT
figure;
imagesc(t_snr_post_den_gift(:,:,80)); %slice 80
colorbar;
colormap(hot); %hot % autum %jet
title('tSNR Map (Post denoising - GIFT)');

% Post denoising - MELODIC
figure;
imagesc(t_snr_post_den_mel(:,:,80)); %slice 80
colorbar;
colormap(hot); %hot % autum %jet
title('tSNR Map (Post denoising - FSL-Melodic)');

% Difference
figure;
imagesc(diff_3d(:,:,80)); %slice 80
colorbar;
colormap(hot); %hot % autum %jet
title('tSNR Map (Difference Melodic - GIFT)');

%% Histogram

%Difference MELODIC - GIFT
figure;
histogram(diff,100, 'FaceColor', [0.3, 0.75, 0.75]);
xlabel('t-snr difference');
ylabel('Number of voxels');
title('t-snr Difference Denoised MELODIC - GIFT');
