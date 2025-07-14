function [delfiles,keepfiles] = my_spmbatch_preprocfasl(ppparams,params,delfiles,keepfiles)

%% Make GM, WM masks
if ~isfield(ppparams.perf(1),'c1m0scanfile') || ~isfield(ppparams.perf(1),'c2m0scanfile') || ~isfield(ppparams.perf(1),'c3m0scanfile')
    [ppparams,delfiles,~] = my_spmbatch_asl_segmentation(ppparams,params,delfiles,keepfiles);
end

%% Control-Label subtraction to make deltam series
if ~isfield(ppparams.perf(1),'deltamfile')
    fprintf('Start subtraction \n')

    [ppparams,delfiles,keepfiles] = my_spmbacth_faslsubtraction(ppparams,params,delfiles,keepfiles);
end

%% Make cbf series
if ~isfield(ppparams.perf(1),'cbffile')
    fprintf('Start CBF mapping \n')
    
    [ppparams,delfiles,keepfiles] = my_spmbatch_fasl_cbfmapping(ppparams,params,delfiles,keepfiles);
end

%% Normalise CBF series
if ~isfield(ppparams.perf(1),'wcbffile')
    fprintf('Do normalization\n')

    [ppparams,delfiles,keepfiles] = my_spmbatch_aslbold_normalization(ppparams,params,delfiles,keepfiles);
end

%% Smooth CBF series        
if ~isfield(ppparams.perf(1),'scbffile')
    fprintf('Do smoothing \n')

    [ppparams,delfiles,keepfiles] = my_spmbatch_dosmoothasl(ppparams,params,delfiles,keepfiles);
end