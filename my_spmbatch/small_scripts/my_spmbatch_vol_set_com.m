function Affine = my_spmbatch_vol_set_com(V)
% use center-of-mass (COM) to roughly correct for differences in the
% position between image and template
% ______________________________________________________________________
% FORMAT:  Affine = cat_vol_set_com(varargin)
%
% V      - mapped images or filenames 
% Affine - affine transformation to roughly correct origin 
% 
% Only if no input is defined the function is called interactively and the
% estimated transformation is applied to the images. Otherwise, only the 
% Affine paramter is returned.
% ______________________________________________________________________
%
% Christian Gaser, Robert Dahnke
% Structural Brain Mapping Group (https://neuro-jena.github.io)
% Departments of Neurology and Psychiatry
% Jena University Hospital
% ______________________________________________________________________
% $Id: 2170 2023-01-26 $


if nargin == 1
  if isstruct(V)
    V = V;
  else
    P = char(V);
    V = spm_vol(P);
  end
else
  P = spm_select(Inf,'image','Select images to filter');
  V = spm_vol(P);
end
n = numel(V);

MM = V(1).private.mat0;

% pre-estimated COM of MNI template
com_reference = [0 -25 -15];

%fprintf('Correct center-of-mass                                            ');
spm_progress_bar('Init',n,'Set center-of-mass','volumes completed');
for i=1:n
  V(i).mat = MM;
  Affine = eye(4);
  if isfield(V(i),'dat')
    vol(:,:,:) = V(i).dat(:,:,:);
  else
    vol = spm_read_vols(V(i));
  end
  avg = mean(vol(:));
  avg = mean(vol(vol>avg));
  
  % don't use background values
  [x,y,z] = ind2sub(size(vol),find(vol>avg));
  com = V(i).mat(1:3,:)*[mean(x) mean(y) mean(z) 1]';
  com = com';

  %M = spm_get_space(V(i).fname);
  Affine(1:3,4) = (com - com_reference)';
  
  if nargin < 1
    spm_get_space(V(i).fname,Affine\MM);
    fprintf('\n');
  end
  
  if ~nargout
    clear Affine
  end
  spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');
