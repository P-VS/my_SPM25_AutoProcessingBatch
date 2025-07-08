function outdat=my_spmbatch_st(dat,Vin,SliceT, TR)

nslices = Vin(1).dim(3);

SliceT = SliceT;

if nslices ~= numel(SliceT)
    error('Mismatch between number of slices and length of ''Slice timings'' vector.');
end

%-Slice timing correction
%==========================================================================

nimgo = numel(Vin);
nimg  = 2^(floor(log2(nimgo))+1);
if Vin(1).dim(3) ~= nslices
    error('Number of slices differ: %d vs %d.', nslices, Vin(1).dim(3));
end

% Compute shifting amount from reference slice and slice timings
% Compute time difference between the acquisition time of the
% reference slice and the current slice by using slice times
% supplied in sliceorder vector
rtime=SliceT(1);
shiftamount = (SliceT - rtime)/TR;

% For loop to perform correction slice by slice

outdat = do_spm_hb_stc(dat,Vin,nimgo,nimg,nslices,shiftamount);

fprintf('%-40s: %30s\n','Completed',spm('time'))  %-#

%%-------------------------------------------------------------------------------------------

function nvol = do_spm_hb_stc(vol,Vin,nimgo,nimg,nslices,shiftamount)

% Set up [time x voxels] matrix for holding image info
nvol=vol;

task = sprintf('Correcting acquisition delay: session %d', 1);
spm_progress_bar('Init',nslices,task,'planes complete');
fprintf('Start Slice Time correction\n')

mask = vol(:,:,:,1)>max(vol(:,:,:,1),[],'all')*0.015;

for k=1:nslices

    slices = vol(:,:,k,:);

    slicemask = mask(:,:,k);
    tmp = find(slicemask>0);

    if numel(tmp)>0

        stack  = zeros([nimg numel(tmp)]);

        rslices = reshape(slices,[Vin(1).dim(1)*Vin(1).dim(2) nimgo]);
        mslices=rslices(tmp,:);
        
        % Set up shifting variables
        len     = size(stack,1);
        phi     = zeros(1,len);
        
        % Check if signal is odd or even -- impacts how Phi is reflected
        %  across the Nyquist frequency. Opposite to use in pvwave.
        OffSet  = 0;
        if rem(len,2) ~= 0, OffSet = 1; end
        
        % Phi represents a range of phases up to the Nyquist frequency
        % Shifted phi 1 to right.
        for f = 1:len/2
            phi(f+1) = -1*shiftamount(k)*2*pi/(len/f);
        end
        
        % Mirror phi about the center
        % 1 is added on both sides to reflect Matlab's 1 based indices
        % Offset is opposite to program in pvwave again because indices are 1 based
        phi(len/2+1+1-OffSet:len) = -fliplr(phi(1+1:len/2+OffSet));
            
        % Transform phi to the frequency domain and take the complex transpose
        shifter = [cos(phi) + sin(phi)*sqrt(-1)].';
        shifter = shifter(:,ones(size(stack,2),1)); % Tony's trick
         
        % Extract columns from slices
        stack(1:nimgo,:) = mslices';
            
        % Fill in continous function to avoid edge effects
        for g=1:size(stack,2)
            stack(nimgo+1:end,g) = linspace(stack(nimgo,g),...
            stack(1,g),nimg-nimgo)';
        end
        
        % Shift the columns
        stack = real(ifft(fft(stack,[],1).*shifter,[],1));
            
        % Re-insert shifted columns
        rslices(tmp,:) = stack(1:nimgo,:)';
        newslices = reshape(rslices,[Vin(1).dim(1:2) nimgo]);

        nvol(:,:,k,:) = newslices;
    end

    spm_progress_bar('Set',k);
end

spm_progress_bar('Clear');
