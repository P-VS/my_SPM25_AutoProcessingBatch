function my_spmbatch_paralleldcm2nii(sublist,params)

pa=parpool(min([3,numel(sublist)])); %25 is the maximum number of workers allowed in the 'local' profile while 10 is set to avoid memory issues on my computer
parfor i = 1:numel(sublist)
    my_spmdcm2niibatch(sublist(i),params,false);
end
delete(pa)