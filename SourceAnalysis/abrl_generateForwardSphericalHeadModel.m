function [gainMatrix] = abrl_generateForwardSphericalHeadModel(chanlist,sphereModel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% abrl_generateForwardSphericalHeadModel
%   - To compute a volume conduction model from the geometry of the head 
% 
% Usage: [gainMatrix] = abrl_generateForwardSphericalHeadModel(chanlist,sphereModel)
% 
% Inputs:
%       chanlist    = cell list of channel labels
%       sphereModel = structure with 
%                       dipole locations
%                       center of head sphere
%                       radii & conductivity the sphere layers
% 
% Outputs:
%       gainMatrix  = EEG forward model gain matrix [nSensors x (3*nDipoles)]
% 
% DISCLAIMER: This function is based on the "bst_eeg_sph.m" function 
%   that is part of the Brainstorm software: https://neuroimage.usc.edu/brainstorm
% 
%   DESCRIPTION:  EEG multilayer spherical forward model
%   This function computes the voltage potential forward gain matrix for an array of 
%   EEG electrodes on the outermost layer of a single/multilayer conductive sphere. 
%   Each region of the multilayer sphere is assumed to be concentric with 
%   isontropic conductivity.  EEG sensors are assumed to be located on the surface
%   of the outermost sphere. 
% 
%   Method: Series Approximiation of a Multilayer Sphere as three dipoles in a 
%           single shell using "Berg/Sherg" parameter approximation.
%   Ref:    Z. Zhang "A fast method to compute surface potentials generated by 
%           dipoles within multilayer anisotropic spheres" 
%           (Phys. Med. Biol. 40, pp335-349,1995)    
% 
%   Dipole generator(s) are assumed to be interior to the innermost "core" layer. For those 
%   dipoles external to the sphere, the dipole "image" is computed and used determine the 
%   gain function. The exception to this is the Legendre Method where all dipoles MUST be 
%   interior to the innermost "core" layer.
% 
%   The default 2000 dipole head sphere is computed pre-hand and will be
%   used in this modified function.
% 
%   Berg parameters are also computed pre-hand
%    [mu_berg, lam_berg] = bst_berg(sphereModel.radii,sphereModel.conductivity);
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original authors: 2000-18 Sylvain Baillet, John Mosher, John Ermer,
%   Francois Tadel (Brainstorm)
% Modified author: Dec 2018; Arun Sasidharan, ABRL
% 
% Copyright (C) 2000-2018 University of Southern California & McGill University
% Copyright (C) 2018 ABRL, Axxonet System Technologies Pvt Ltd., Bengaluru, India
% 
% Permission is hereby granted, free of charge, to any person obtaining a 
% copy of this software and associated documentation files (the "Software"), 
% to deal in the Software without restriction, including without limitation 
% the rights to use, copy, modify, merge, publish, distribute, sublicense, 
% and/or sell copies of the Software, and to permit persons to whom the 
% Software is furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included 
% in all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS 
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL 
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
% OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
% ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
% OTHER DEALINGS IN THE SOFTWARE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get channel locations
[chanlocs,chan_absent] = abrl_getChanLocs(chanlist);
chanlocs(chan_absent)  = [];


n_layers     = length(sphereModel.center);   % # of concentric sphere layers
n_dipoles    = size(sphereModel.dipoles,1);
n_electrodes = length(chanlocs);

% Extract x,y,z coordinates of electrodes, scale it to head radius and convert to meters
electrodes = [chanlocs.headpoints]';

% Center electrode coordinates on the center of the sphere
electrodes  = bsxfun(@minus, electrodes, sphereModel.center(:)');

% Make sure EEG sensors project on the head sphere
[elec_theta,elec_phi,elec_sph]  = cart2sph(electrodes(:,1),electrodes(:,2),electrodes(:,3));
elec_sph = sphereModel.radii(n_layers) * ones(size(elec_sph));
[electrodes(:,1),electrodes(:,2),electrodes(:,3)] = sph2cart(elec_theta,elec_phi,elec_sph);

% electrodes   = [chanlocs.X;chanlocs.Y;chanlocs.Z]';
% electrodes   = bsxfun(@times, electrodes/1000, sphereModel.radii(n_layers)/(chanlocs(1).sph_radius/1000));
% [x,y,z] = sph2cart([chanlocs.sph_theta]',[chanlocs.sph_phi]',repmat(sphereModel.radii(n_layers),n_electrodes,1));
% [x,y,z] = sph2cart([chanlocs.sph_theta]',[chanlocs.sph_phi]',[chanlocs.sph_radius]');
% electrodes = [x,y,z];

% Extract x,y,z, coordinates of dipoles
dipoles     = sphereModel.dipoles;

% Center dipole coordinates on the center of the sphere
dipoles     = bsxfun(@minus, dipoles, sphereModel.center(:)');

% Make sure EEG sensors project on the head sphere
% [elec_theta elec_phi elec_sph]  = cart2sph(electrodes(:,1),electrodes(:,2),electrodes(:,3));
% elec_sph = sphereModel.radii(n_layers) * ones(size(elec_sph));
% [electrodes(:,1),electrodes(:,2),electrodes(:,3)] = sph2cart(elec_theta,elec_phi,elec_sph);

% Pre-Allocate Gain Matrix
gainMatrix = zeros(n_electrodes,3*n_dipoles);


electrode_mag        = sphereModel.radii(n_layers);
dipole_mag           = sqrt(sum(dipoles.^2, 2));% Calculate the Euclidean norm
dipole_dot_electrode = dipoles*electrodes';

for n_berg = 1:length(sphereModel.mu_berg)
    mu = sphereModel.mu_berg(n_berg);
    % This part checks for the presence of Berg dipoles which are external to
    %  the sphere. For those dipoles external to the sphere, the dipole parameters
    %  are replaced with the electrical image (internal to sphere) of the dipole
    dipoleindx = find(mu .* dipole_mag > electrode_mag);
    if ~isempty(dipoleindx)
        warning('Check results...');
        dipoles(dipoleindx,:) = electrode_mag.^2 * bsxfun(@rdivide, dipoles(dipoleindx,:), sum(dipoles(dipoleindx,:).^2,2));
    end
    
    % Calculation of Forward Gain Matrix Contribution due to K-th Berg Dipole
    dipole_mag_sq = repmat((mu * dipole_mag).^2,1,n_electrodes);          %(PxM)
    const = 1 ./ (4.0 * pi * sphereModel.conductivity(n_layers) / sphereModel.lam_berg(n_berg) * mu.^2 * dipole_mag.^2);
    d_matrix = reshape(repmat(electrodes,1,n_dipoles)',3,n_dipoles*n_electrodes)' - mu .* repmat(dipoles,n_electrodes,1);
    d_mag = reshape(sqrt(sum(d_matrix.^2, 2)),n_dipoles,n_electrodes);% Calculate the Euclidean norm
    F_scalar = d_mag .* (electrode_mag.*d_mag + electrode_mag.^2 - mu.*dipole_dot_electrode); %(PxM)
    c1 = bsxfun(@times, (2*((mu.*dipole_dot_electrode - dipole_mag_sq) ./ d_mag.^3) + 1./d_mag - 1./electrode_mag), const);%(PxM)
    c2 = bsxfun(@times, (2./d_mag.^3) + (d_mag + electrode_mag) ./ (electrode_mag.*F_scalar), const);
    gainMatrix = gainMatrix + reshape(repmat((c1 - c2.*mu.*dipole_dot_electrode)',3,1)...
        ,n_electrodes,3*n_dipoles) .* mu .* repmat(reshape(dipoles',1,3*n_dipoles),n_electrodes,1) ...
        + reshape(repmat((c2.*dipole_mag_sq)',3,1),n_electrodes,3*n_dipoles) .* repmat(electrodes,1,n_dipoles);
end

end