% =======================================================================================
% (c) Lars Ruthotto and Siawoosh Mohammadi 2013
% http://www.eos.ubc.ca/about/researcher/L.Ruthotto.html
% 
% function acid_hysco_write(POI1,POI2,Bc,pe_direction)
%
% Applying HySCO result to other blip-up or  driver for HySCO (Hyperelastic Susceptibility COrrection of DTI)
%
% Input:
%
%  POI1         - matrix of filenames for additional blip-up volumes
%  POI2         - matrix of filenames for additional blip-down volumes
%  VB           - filename of inhomogeneity estimate produced by HySCO
%  pe_direction - phase-encoding direction, 1 for x_1, 2 for x_2, 3 for x_3
%                 (data dimensions will be flipped accordingly)
%
% Please cite one of the following works when using this software
%
% @inproceedings{Ruthotto2013,
%   author    = {Ruthotto, L and Mohammadi, S and Heck, C and Modersitzki, J and Weiskopf, N},
%   title     = {HySCO - Hyperelastic Susceptibility Artifact Correction of DTI in SPM}},
%   booktitle = {Bildverarbeitung f{\"u}r die Medizin 2013},
%   year      = {2013}
% }
%
% @article{Ruthotto2012,
%   author  = {Ruthotto, L and Kugel, H and Olesch, J and Fischer, B and Modersitzki, J and Burger, M and Wolters, CH},
%   title   = {Diffeomorphic Susceptibility Artefact Correction of Diffusion-Weighted Magnetic Resonance Images}},
%   journal = {Physics in Medicine and Biology},
%   volume  = {57},
%   number  = {18},
%   pages   = {5715--5731}
%   year    = {2012}
% }
%
% =======================================================================================

function outdat = my_spmbatch_hysco_write(PVG1,in_blipup,in_blipdown,PB,pe_direction,dummy_3dor4d)
VG1  = spm_vol(PVG1);

% default reslice parameter
if(~exist('res','var'))
    res = -4;
end

% extract data resolution and domain info 
% 
% Note that domain is assumed to be rectangular and aligned to the coordinate system,
% i.e. omega = [omega(1),omega(2)] x [omega(3),omega(4)] x [omega(5),omega(6)]

m     = VG1.dim;
Vmat  = sqrt(sum(VG1.mat(1:3,1:3).^2)); % modified by SM to make sure that voxel size is kept
    
omega = zeros(1,6);
omega(2:2:end) = Vmat(1:3).*m; % modified by SM to make sure that voxel size is kept

% find out if inhomogeneity is nodal or staggered
Bc = spm_read_vols(spm_vol(PB)); 
isNodal = numel(Bc)==prod(m+1);
mstg    = m; mstg(pe_direction) = mstg(pe_direction)+1;
isStg   = numel(Bc)==prod(mstg);
if not(isNodal) && not(isStg)
    error('resolution of inhomogeneity not compatible with data');
end


% permute data dimensions such that phase encoding is along second index
if isNodal
    switch pe_direction
        case 1
            read_data   = @(str) permute(sm_read_vols(str,VG1,res),[2 1 3]);
            write_data  = @(A,V,prefix) spm_write_image(ipermute(A,[2 1 3]),V,prefix,VG1);
            write_data1 = @(A,V,VG,prefix,Ni,Nvol,sz) spm_write_image_4d(ipermute(A,[2 1 3]),V,VG,prefix,Ni,Nvol,sz);
            omega=omega([3 4 1 2 5 6]);  m= m([2 1 3]);
            vecperm = [2 1 3];
        case 2
            read_data   = @(str) sm_read_vols(str,VG1,res);
            write_data  = @(A,V,prefix) spm_write_image(A,V,prefix,VG1);
            write_data1 = @(A,V,VG,prefix,Ni,Nvol,sz) spm_write_image_4d(A,V,VG,prefix,Ni,Nvol,sz);
            vecperm = [1 2 3];
        case 3
            read_data   = @(str) permute(sm_read_vols(str,VG1,res),[1 3 2]);
            write_data  = @(A,V,prefix) spm_write_image(ipermute(A,[1 3 2]),V,prefix,VG1);
            write_data1 = @(A,V,VG,prefix,Ni,Nvol,sz) spm_write_image_4d(ipermute(A,[1 3 2]),V,VG,prefix,Ni,Nvol,sz);
            omega=omega([1 2 5 6 3 4]);  m=m([1 3 2]);
            vecperm = [1 3 2];
    end
elseif isStg
    switch pe_direction
        case 1
            read_data   = @(str) sm_read_vols(str,VG1(1),res);
            write_data  = @(A,V,prefix) spm_write_image(A,V,prefix,VG1);
            write_data1 = @(A,V,VG,prefix,Ni,Nvol,sz) spm_write_image_4d(A,V,VG,prefix,Ni,Nvol,sz);
            vecperm = [1 2 3];
        case 2
            read_data   = @(str) permute(sm_read_vols(str,VG1(1),res),[2 1 3]);
            write_data  = @(A,V,prefix) spm_write_image(ipermute(A,[2 1 3]),V,prefix,VG1);
            write_data1 = @(A,V,VG,prefix,Ni,Nvol,sz) spm_write_image_4d(ipermute(A,[2 1 3]),V,VG,prefix,Ni,Nvol,sz);
            omega=omega([3 4 1 2 5 6]);  m= m([2 1 3]);
            vecperm = [2 1 3];
        case 3
            read_data   = @(str) permute(sm_read_vols(str,VG1(1),res),[3 1 2]);
            write_data  = @(A,V,prefix) spm_write_image(ipermute(A,[3 1 2]),V,prefix,VG1);
            write_data1 = @(A,V,VG,prefix,Ni,Nvol,sz) spm_write_image_4d(ipermute(A,[3 1 2]),V,VG,prefix,Ni,Nvol,sz);
            omega=omega([ 5 6 1:4]);  m=m([3 1 2]);
            vecperm = [3 1 2];
    end
end

% save inhomogeneity
Bc = permute(spm_read_vols(spm_vol(PB)),vecperm); 

% compute transformations and intensity modulations
if isNodal
    y1    = acid_hysco_getTrafoEPI(Bc,[0;1;0],omega,m,'matrixFree',1);
    y2    = acid_hysco_getTrafoEPI(Bc,[0;-1;0],omega,m,'matrixFree',1);
    y1    = nodal2center(y1,m);
    y2    = nodal2center(y2,m);
    pB    = acid_hysco_getPartialB(Bc,omega,m,'cc','matrixFree',1);
elseif isStg
    xc    = reshape(getCellCenteredGrid(omega,m),[],3);
    Bc    = reshape(Bc,m+[1,0,0]);                  % 1-staggered
    Bcc   = .5*(Bc(1:end-1,:,:) + Bc(2:end,:,:));   % cell-centered
    y1    = xc; y1(:,1) = y1(:,1) + Bcc(:);
    y2    = xc; y2(:,1) = y2(:,1) - Bcc(:);
    h     = (omega(2:2:end)-omega(1:2:end))./m;
    pB    = (Bc(2:end,:,:) - Bc(1:end-1,:,:))/h(1); % partial derivative
end
Jac1  = 1 + pB;
Jac2  = 1 - pB;

% save corrected data. prefix is motivated by fieldmap toolbox
prefix = 'u';

% number of image volumes does not agree. apply best-known
% field-inhomogeneity to image volumes
if ~isempty(in_blipup)
    dim = size(in_blipup);
    if numel(dim)<4
        dim=[dim 1]; 
        in_blipup = reshape(in_blipup,[dim(1),dim(2),dim(3),1]);
    end
    outdat = zeros(dim);

    for vol=1:dim(4)
        fprintf('Apply estimated field inhomogeneity to blip-up volume %d\n',vol);
        I1 = in_blipup(:,:,:,vol);
        outdat(:,:,:,vol) = reshape( linearInterMex(I1, omega,y1).*Jac1(:) ,m);
    end
end

if ~isempty(in_blipdown)
    dim = size(in_blipdown);
    if numel(dim)<4
        dim=[dim 1]; 
        in_blipdown = reshape(in_blipdown,[dim(1),dim(2),dim(3),1]);
    end
    outdat = zeros(dim);
    
    for vol=1:dim(4)
        fprintf('Apply estimated field inhomogeneity to blip-down volume %d\n',vol);
        I2 = in_blipdown(:,:,:,vol);
        outdat(:,:,:,vol) = reshape( linearInterMex(I2, omega,y1).*Jac1(:) ,m);
    end
end



%{
    (c) Lars Ruthotto and Jan Modersitzki 2013

    This file is part of HySCO (Version 1.0, 2013/03/28)
                           -  Hyperelastic Susceptibility Artefact Correction for DTI

    
    HySCO is free but copyright software, distributed under the terms of the 
    GNU General Public Licence as published by the Free Software Foundation 
    (Version 3, 29 June 2007) http://www.gnu.org/licenses/gpl.html

 
    This code is provided "as is", without any warranty of any kind, either
    expressed or implied, including but not limited to, any implied warranty
    of merchantibility or fitness for any purpose. In no event will any party
    who distributed the code be liable for damages or for any claim(s) by
    any other party, including but not limited to, any lost profits, lost
    monies, lost data or data rendered inaccurate, losses sustained by
    third parties, or any other special, incidental or consequential damages
    arising out of the use or inability to use the program, even if the
    possibility of such damages has been advised against. The entire risk
    as to the quality, the performace, and the fitness of the program for any
    particular purpose lies with the party using the code.

    This code is especially not intended for any clinical or diagnostic use. 
  
%}
