function [tefasldat,Vout,ppparams,delfiles] = my_spmbatch_aslremovebold(tefasldata,ppparams,params,delfiles,is_asl)
%% Loading the ASL time series

nechoes = numel(params.func.echoes);

GM = spm_vol(fullfile(ppparams.subperfdir,ppparams.asl(1).c1m0scanfile));
WM = spm_vol(fullfile(ppparams.subperfdir,ppparams.asl(1).c2m0scanfile));
CSF = spm_vol(fullfile(ppparams.subperfdir,ppparams.asl(1).c3m0scanfile));

gmim = spm_read_vols(GM);
wmim = spm_read_vols(WM);
csfim = spm_read_vols(CSF);

mask = my_spmbatch_mask(tefasldata{params.func.echoes(1)}.data);
mask = gmim+wmim;
mask(gmim+wmim<0.1) = 0;
mask(mask>0) = 1;

mask_ind = find(mask>0);

csfim(gmim+wmim>0.1) = 0;

voldim = size(tefasldata{params.func.echoes(1)}.data);
Vasl = tefasldata{params.func.echoes(1)}.Vasl;

conidx = 2:2:numel(Vasl);
labidx = 1:2:numel(Vasl);

csfdata = sum(reshape(tefasldata{params.func.echoes(1)}.data .* repmat(csfim,1,1,1,numel(Vasl)),[voldim(1)*voldim(2)*voldim(3),numel(Vasl)]),1)/numel(find(csfim>0));

mean_csfcon = mean(csfdata(conidx));
mean_csflab = mean(csfdata(labidx));

if mean_csflab>mean_csfcon
    conidx = 1:2:numel(Vasl);
    labidx = 2:2:numel(Vasl);
end

switch params.asl.removebold
    case 't2star'
        for i=1:nechoes
            jsondat = fileread(ppparams.asl(params.func.echoes(i)).asljson);
            jsondat = jsondecode(jsondat);
        
            te(i) = 1000.0*jsondat.EchoTime;
        
            Vasl = tefasldata{params.func.echoes(i)}.Vasl;
        
            if i==1
                voldim = Vasl.dim;
                tefasldat = zeros(voldim(1),voldim(2),voldim(3),numel(Vasl),nechoes);
                utefasldat = tefasldat;
            end
        
            utefasldat(:,:,:,:,i) = tefasldata{params.func.echoes(i)}.data;
            otefasldat = utefasldat;
        
            for iv=1:numel(Vasl)
                otefasldat(:,:,:,iv,i) = my_spmbatch_smooth(utefasldat(:,:,:,iv,i),Vasl(iv),{},[3 3 3],0);
            end
        
            %if is_asl    
            %    if contains(params.asl.tagorder,'labeled')
            %        cidx=[-3 -2 -1 0 1 2];
            %        conidx = 2:2:numel(Vasl);
            %        labidx = 1:2:numel(Vasl);
            %    else
            %        cidx=[-2 -1 0 1 2 3];
            %        conidx = 1:2:numel(Vasl);
            %        labidx = 2:2:numel(Vasl);
            %    end
        
            %    for iv=1:numel(Vasl) % sinc interpolation to replace the labeled scans by a control scan
            %        if isempty(find(conidx==iv)), timeshift = 2.5; else timeshift = 3; end
            %        idx=ceil(iv/2)+cidx;
            %        idx(find(idx<1))=1;
            %        idx(find(idx>numel(conidx)))=numel(conidx);
            %        nimg = otefasldat(:,:,:,conidx(idx),i);
            %        nimg=reshape(nimg,size(nimg,1)*size(nimg,2)*size(nimg,3),size(nimg,4));
            %        clear tmpimg;
            %        [pn,tn]=size(nimg);
            %        tmpimg=sinc_interpVec(nimg(mask_ind,:),timeshift);
            %        Vconimg=zeros(size(nimg,1),1);
            %        Vconimg(mask_ind)=tmpimg;
            %        Vconimg=reshape(Vconimg,voldim(1),voldim(2),voldim(3));
            %        clear tmpimg pn tn;
            
            %        otefasldat(:,:,:,iv,i) = Vconimg;
            %    end
            %end
        end
        
        tefasldat = zeros(voldim(1),voldim(2),voldim(3),numel(Vasl));
        t2stardat = zeros(voldim(1),voldim(2),voldim(3),numel(Vasl));
        
        %based on https://github.com/jsheunis/fMRwhy/tree/master
        for ti=1:numel(Vasl)
        
            tifasldat = reshape(otefasldat(:,:,:,ti,:),[voldim(1),voldim(2),voldim(3),nechoes]);
        
            % Create "design matrix" X
            X = horzcat(ones(nechoes,1), -te(:));
        
            t2star = zeros(voldim(1)*voldim(2)*voldim(3),1);
        
            Y=[];
            for ne=1:nechoes
                temptefasldat = reshape(tifasldat(:,:,:,ne),[voldim(1)*voldim(2)*voldim(3),1]);
                Y=[Y;reshape(temptefasldat(mask_ind,1),[1,numel(mask_ind)])];
            end
            Y = max(Y, 1e-11);
        
            % Estimate "beta matrix" by solving set of linear equations
            beta_hat = pinv(X) * log(Y);
            % Calculate S0 and T2star from beta estimation
            lnS0_fit = beta_hat(1, :); %is lnS0
            T2star_fit = beta_hat(2, :); %is R2*
        
            T2star_thresh_min = 1/1500; % arbitrarily chosen, same as tedana
            I_T2star_min = (T2star_fit < T2star_thresh_min); % vector of voxels where T2star value is negative
            T2star_fit(I_T2star_min) = 0; % if values inside mask are zero or negative, set them to threshold_min value
        
            t2star(mask_ind) = T2star_fit;
        
            weights = zeros(voldim(1)*voldim(2)*voldim(3),nechoes);
        
            for ne=1:nechoes
                weights(:,ne) = repmat(te(ne),voldim(1)*voldim(2)*voldim(3),1) .* t2star(:,1);
                weights(mask_ind,ne) = exp(weights(mask_ind,ne));
            end
        
            weights = reshape(weights,[voldim(1),voldim(2),voldim(3),nechoes])/3;
            
            sum_weights = sum(weights,4);
            weights_mask = find(sum_weights>0);
        
            for ne=1:nechoes
                fasltidat = utefasldat(:,:,:,ti,ne);
                fasltidat = fasltidat .* weights(:,:,:,ne);
                fasltidat(fasltidat<0) = 0;
        
                tefasldat(:,:,:,ti) = tefasldat(:,:,:,ti)+fasltidat; 
            end 
        
            t2stardat(:,:,:,ti) = reshape(t2star,[voldim(1),voldim(2),voldim(3)]);
        end
        
        if params.save_intermediate_results
            Vasl = tefasldata{params.func.echoes(1)}.Vasl;
            
            if is_asl
                nfname = split(ppparams.asl(1).aslfile,'_echo-');
                endfix = '_asl';
            else
                nfname = split(ppparams.asl(1).m0scanfile,'_echo-');
                endfix = '_m0scan';
            end
            
            Vout2 = Vasl;
            
            for j=1:numel(Vout2)
                Vout2(j).fname = fullfile(ppparams.subperfdir,['t2' tefasldata{params.func.echoes(1)}.prefix nfname{1} endfix '.nii']);
                Vout2(j).descrip = 'my_spmbatch - T2star';
                Vout2(j).pinfo = [1,0,0];
                Vout2(j).dt = [spm_type('float32'),spm_platform('bigend')];
                Vout2(j).n = [j 1];
            end
            
            Vout2 = myspm_write_vol_4d(Vout2,t2stardat);
    
            delfiles{numel(delfiles)+1} = {Vout2.fname};
        end

        dt = [spm_type('float32'),spm_platform('bigend')];
    case 'filter'
        jsondat = fileread(ppparams.asl(params.func.echoes(1)).asljson);
        jsondat = jsondecode(jsondat);
    
        tr = jsondat.RepetitionTime;
        
        for i=1:nechoes
            ttefasldat = tefasldata{params.func.echoes(i)}.data;

            hpf = min([0.09,1/(2*tr)]);

            if is_asl
                s = size(ttefasldat);
                ttefasldat = reshape(ttefasldat(:,:,:,:),[prod(s(1:end-1)),s(end)]);
            
                [ttefasldat,~] = fmri_cleaning(ttefasldat(:,:),0,[tr hpf Inf],[],[],'restoremean','on');
                
                if i==1, tefasldat = reshape(ttefasldat(:,:),s); else tefasldat = tefasldat + reshape(ttefasldat(:,:),s); end
            else 
                if i==1, tefasldat = reshape(ttefasldat(:,:),s); else tefasldat = tefasldat + reshape(ttefasldat(:,:),s); end
            end     

            tefasldat = tefasldat/nechoes;

            clear 'ttefasldat'
        end

        dt = tefasldata{params.func.echoes(1)}.Vasl(1).dt;
end

Vout = tefasldata{params.func.echoes(1)}.Vasl;

if is_asl
    nfname = split(ppparams.asl(params.func.echoes(1)).aslfile,'_echo-');
    endfix = '_asl';

    ppparams.asl(1).caslfile = ['c' tefasldata{params.func.echoes(1)}.prefix nfname{1} '_asl.nii'];
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,ppparams.asl(1).caslfile)};

    ppparams.asl(1).aslfile = [nfname{1} '_asl.nii'];
    ppparams.asl(1).aslprefix = ['c' tefasldata{params.func.echoes(1)}.prefix];
    prefix = ppparams.asl(1).aslprefix;
else
    nfname = split(ppparams.asl(params.func.echoes(1)).m0scanfile,'_echo-');
    endfix = '_m0scan';

    ppparams.asl(1).cm0scanfile = ['c' tefasldata{params.func.echoes(1)}.prefix nfname{1} '_m0scan.nii'];
    delfiles{numel(delfiles)+1} = {fullfile(ppparams.subperfdir,ppparams.asl(1).cm0scanfile)};

    ppparams.asl(1).m0scanfile = [nfname{1} '_m0scan.nii'];
    ppparams.asl(1).m0prefix = ['c' tefasldata{params.func.echoes(1)}.prefix];
    prefix = ppparams.asl(1).m0prefix;
end

if params.save_intermediate_results
    for j=1:numel(Vout)
        Vout(j).fname = fullfile(ppparams.subperfdir,[prefix nfname{1} endfix '.nii']);
        Vout(j).descrip = 'my_spmbatch - BOLD removed';
        Vout(j).pinfo = [1,0,0];
        Vout(j).dt = dt;
        Vout(j).n = [j 1];
    end
    
    Vout = myspm_write_vol_4d(Vout,tefasldat);
end