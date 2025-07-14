function my_spmbatch_start_vbmprocessing(sublist,nsessions,datpath,params)

save(fullfile(datpath,'params.mat'),'params')

datlist = zeros(numel(sublist)*numel(nsessions),2);

dpos = 1;
for i = 1:numel(sublist)
    for j = 1:numel(nsessions)
        datlist(dpos,1) = sublist(i);
        datlist(dpos,2) = nsessions(j);

        dpos = dpos+1;
    end
end

if params.use_parallel
    numpacks = ceil(numel(datlist(:,1))/params.maxprocesses);

    for j=1:numpacks
        if (j*params.maxprocesses)<=numel(datlist(:,1))
            maxruns = params.maxprocesses;
        else
            maxruns = params.maxprocesses-((j*params.maxprocesses)-numel(datlist(:,1)));
        end

        for is = 1:maxruns
            i = (j-1)*params.maxprocesses+is;

            fprintf(['\nStart VBM preprocessing data for subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) '\n'])
    
            logfile{i} = fullfile(datpath,['logfile_' sprintf(['%0' num2str(params.sub_digits) 'd'],datlist(i,1)) '_' sprintf('%02d',datlist(i,2)) '.txt']);
    
            if exist(logfile{i},'file'), delete(logfile{i}); end
            
            if ispc
                mtlb_cmd = sprintf("restoredefaultpath;addpath(genpath('%s'));addpath(genpath('%s'));my_spmbatch_run_vbmpreprocessing(%d,%d,'%s','%s');exit", ...
                                    params.spm_path,params.my_spmbatch_path,datlist(i,1),datlist(i,2),datpath,fullfile(datpath,'params.mat'));

                system_cmd = sprintf(['start matlab -nodesktop -nosplash -r "%s" -logfile %s'],mtlb_cmd,logfile{i});
            else
                mtlb_cmd = sprintf('"restoredefaultpath;addpath(genpath(''%s''));addpath(genpath(''%s''));my_spmbatch_run_vbmpreprocessing(%d,%d,''%s'',''%s'');exit"', ...
                                    params.spm_path,params.my_spmbatch_path,datlist(i,1),datlist(i,2),datpath,fullfile(datpath,'params.mat'));

                system_cmd = sprintf([fullfile(matlabroot,'bin') '/matlab -nosplash -r ' mtlb_cmd ' -logfile ' logfile{i} ' & ']);
            end
            [status,result]=system(system_cmd);
        end
    
        %% wait for all processing to be finnished
        isrunning = true;
        pfinnished = 0;
        while isrunning
            for is = 1:maxruns
                i = (j-1)*params.maxprocesses+is;

                if exist(logfile{i},'file')
                    FID     = fopen(logfile{i},'r');
                    txt     = textscan(FID,'%s');
                    txt     = txt{1}; 
                    test=find(cellfun('isempty',strfind(txt,'PP_Completed'))==0,1,'first');
                    errortest=find(cellfun('isempty',strfind(txt,'PP_Error'))==0,1,'first');
                    fclose(FID);

                    if ~isempty(errortest)
                        pfinnished = pfinnished+1;

                        nlogfname = fullfile(datpath,['error_vbm_logfile_' sprintf(['%0' num2str(params.sub_digits) 'd'],datlist(i,1)) '_' sprintf('%02d',datlist(i,2)) '.txt']);
                        movefile(logfile{i},nlogfname);

                        fprintf(['\nError during VBM processing data for subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) '\n'])
                    elseif ~isempty(test)
                        pfinnished = pfinnished+1;

                        if ~params.keeplogs
                            delete(logfile{i}); 
                        else
                            nlogfname = fullfile(datpath,['done_vbm_logfile_' sprintf(['%0' num2str(params.sub_digits) 'd'],datlist(i,1)) '_' sprintf('%02d',datlist(i,2)) '.txt']);
                            movefile(logfile{i},nlogfname);
                        end

                        fprintf(['\nDone VBM preprocessing data for subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) '\n'])
                    end
                end
            end
    
            if pfinnished>=maxruns 
                isrunning = false; 
            else
                pause(60);
            end
        end
    end
else
    for i=1:numel(datlist(:,1))
        itstart = tic;

        my_spmbatch_run_vbmpreprocessing(datlist(i,1),datlist(i,2),datpath,fullfile(datpath,'params.mat'));

        itstop = toc(itstart);

        fprintf(['subject ' num2str(datlist(i,1)) ' session ' num2str(datlist(i,2)) ' processed in ' datestr(duration([0,0,itstop],'InputFormat','ss'),'HH:MM:SS') '\n'])
    end
end

delete(fullfile(datpath,'params.mat'))