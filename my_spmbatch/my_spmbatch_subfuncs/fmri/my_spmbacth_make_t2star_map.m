function t2star = my_spmbacth_make_t2star_map(tefuncdat,te)
%based on https://github.com/jsheunis/fMRwhy/tree/master

voldim = size(tefuncdat);
if numel(voldim)<5, tefuncdat = reshape(tefuncdat,[voldim(1),voldim(2),voldim(3),1,voldim(4)]); end
voldim = size(tefuncdat);

nechoes = numel(te);

mask = zeros([voldim(1),voldim(2),voldim(3)]);

for ie=1:nechoes
    iemask = my_spmbatch_mask(tefuncdat(:,:,:,:,ie));
    mask = mask + iemask;

    clear iemask
end

mask_ind = find(mask>nechoes-0.5);

t2star = zeros(voldim(1),voldim(2),voldim(3),voldim(4));

for ti=1:voldim(4)
    % Create "design matrix" X
    X = horzcat(ones(nechoes,1), -te(:));
    
    tempt2star = zeros(voldim(1)*voldim(2)*voldim(3),1);
    
    Y=[];
    for ne=1:nechoes
        temptefuncdat = reshape(tefuncdat(:,:,:,ti,ne),[voldim(1)*voldim(2)*voldim(3),1]);
        Y=[Y;reshape(temptefuncdat(mask_ind,1),[1,numel(mask_ind)])];
    end
    Y = max(Y, 1e-11);
    
    % Estimate "beta matrix" by solving set of linear equations
    beta_hat = pinv(X) * log(Y);
     % Calculate S0 and T2star from beta estimation
    T2star_fit = beta_hat(2, :); %is R2*
    
    tempt2star(mask_ind) = T2star_fit;
    tempt2star(mask_ind) = T2star_fit;
    zeromask = (tempt2star>0);
    
    tempt2star(zeromask) = 1 ./ tempt2star(zeromask);
    
    t2star_perc = prctile(tempt2star,99.5,'all');
    T2star_thresh_max = 10 * t2star_perc; %same as tedana
    tempt2star(tempt2star>T2star_thresh_max) = t2star_perc;
    t2star(:,:,:,ti) = reshape(tempt2star,[voldim(1),voldim(2),voldim(3)]);
    
    clear temptefuncdat tempt2star T2star_fit Y beta_hat
end

if numel(size(tefuncdat))<5, t2star = reshape(t2star,[voldim(1),voldim(2),voldim(3)]); end

clear mask mask_ind