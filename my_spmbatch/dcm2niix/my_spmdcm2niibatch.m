function my_spmdcm2niibatch(sub,params,useparfor)

%Code based on dcm2nii from
%(https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)
%(https://www.mathworks.com/matlabcentral/fileexchange/42997-xiangruili-dicm2nii)

substring = ['sub-' num2str(sub,['%0' num2str(params.sub_digits) 'd'])];

for si=1:numel(params.mridata)
    try
        sesstring = ['ses-' num2str(params.mridata(si).session,'%03d')];
    
        subpath = fullfile(params.datpath,substring,sesstring);
        
        infolder = fullfile(params.datpath,substring,params.mridata(si).folder);
        outfolder = fullfile(subpath,params.mridata(si).acqtype);
    
        if ~isfolder(outfolder)
            mkdir(outfolder);
        end
    
        %% nii file name
        outfname = [substring '_' sesstring];
        if contains(params.mridata(si).seqtype,'fmri') || contains(params.mridata(si).seqtype,'aslbold'), outfname=[outfname '_task-' params.mridata(si).task]; end
        if params.mridata(si).add_run, outfname=[outfname '_run-' num2str(params.mridata(1).run)]; end
    
        dirlist = dir(fullfile(infolder,['**' filesep '*.*']));
    
        tmp = find(strlength({dirlist.name})>4); %Remove '.' and '..'
        if ~isempty(tmp), dirlist = dirlist(tmp); end
           
        tmp = find(~contains({dirlist.name},'._')); %Remove the hiden files from Mac from the list
        if ~isempty(tmp), dirlist = dirlist(tmp); end
    
        dirlist = [dir(fullfile(infolder,['**' filesep '*MRDC*.*'])),...
                    dir(fullfile(infolder,['**' filesep '*.dcm'])),...
                    dir(fullfile(infolder,['**' filesep '*.DCM']))];

        if isempty(dirlist)
            tmp = find(contains({dirlist.name},'.zip'));    
            if ~isempty(tmp)
                for iz=1:numel(tmp)
                    fname = split(dirlist(tmp(iz)).name,'.zip');
                    if ~isfolder(fullfile(dirlist(tmp(iz)).folder,fname{1}))
                        unzip(fullfile(dirlist(tmp(iz)).folder,dirlist(tmp(iz)).name),infolder)
                    end
                end
            end
    
            dirlist = [dir(fullfile(infolder,['**' filesep '*MRDC*.*'])),...
                        dir(fullfile(infolder,['**' filesep '*.dcm'])),...
                        dir(fullfile(infolder,['**' filesep '*.DCM']))];

            if ~isempty(tmp)
                fprintf('No dicom files found')
                return
            end
        end
    
        my_spmdicm2niix(infolder, outfolder, '.nii', outfname, useparfor);
    
        namefilters(1).name = outfname;
        namefilters(1).required = true;
    
        niidirlist = my_spmbatch_dirfilelist(outfolder,'.nii',namefilters,false);
    
        for ifile=1:numel({niidirlist.name})
            fname = split(niidirlist(ifile).name,'.nii');
    
            if contains(fname{1},'_e'), params.mridata(si).add_echo = true; end
    
            niifile = fullfile(outfolder,[fname{1} '.nii']);
            jsonfile = fullfile(outfolder,[fname{1} '.json']);
    
            noutfname = outfname;
        
            if params.mridata(si).add_acq || params.mridata(si).add_dir
                jsondat = fileread(jsonfile);
                jsondat = jsondecode(jsondat);
        
                if params.mridata(si).add_acq
                    acq = jsondat.SeriesDescription;
        
                    noutfname=[noutfname '_acq-' acq];
                end
        
                if params.mridata(si).add_dir
                    jsonpedir = jsondat.PhaseEncodingDirection;
            
                    if contains(jsonpedir,'i'), pedir = 'LR'; end
                    if contains(jsonpedir,'j'), pedir = 'PA'; end
                    if contains(jsonpedir,'k'), pedir = 'FH'; end
                    if contains(jsonpedir,'-'), pedir = reverse(pedir); end
        
                    noutfname=[noutfname '_dir-' pedir];
                end
            end
    
            if params.mridata(si).add_echo
                if contains(fname{1},'_e')
                    esplit = split(fname{1},'_e');
                    numecho = esplit{end}(1);
                else numecho = num2str(ifile); end
    
                noutfname=[noutfname '_echo-' numecho];
            end
        
            if contains(params.mridata(si).seqtype,'T1w'), noutfname=[noutfname '_T1w']; end
            if contains(params.mridata(si).seqtype,'pepolar'), noutfname=[noutfname '_epi']; end
            if contains(params.mridata(si).seqtype,'fmri'), noutfname=[noutfname '_bold']; end
            if contains(params.mridata(si).seqtype,'aslbold'), noutfname=[noutfname '_aslbold']; end
        
            movefile(niifile,fullfile(outfolder,[noutfname '.nii']));
            movefile(jsonfile,fullfile(outfolder,[noutfname '.json']));
        end
    catch
    end
end