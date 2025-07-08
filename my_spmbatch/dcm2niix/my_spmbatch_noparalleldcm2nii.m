function my_spmbatch_noparalleldcm2nii(sublist,params)

for i = sublist
    my_spmdcm2niibatch(i,params,params.use_parallel);
end