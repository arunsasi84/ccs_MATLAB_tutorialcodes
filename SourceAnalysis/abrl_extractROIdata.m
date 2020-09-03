function [roi_data,roi_names] = abrl_extractROIdata(source_data,sphereModel,roi_function)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% abrl_extractROIdata
%   - To extract time series from  
% 
% Usage: [leadFieldMatrix] = abrl_sLORETA(chanlist,gainMatrix)
% 
% Inputs:
%       source_data  = dipole time series [(3*nDipoles) x timepoints]
%       sphereModel  = structure with 
%                       dipole locations
%                       center of head sphere
%                       radii & conductivity the sphere layers
%       roi_function = function to compile ROI from dipole time series
%                       'mean' [default]
%                       'pca' (first pca component)
% 
% Outputs:
%       roi_data     = ROI time series [nROIs x timepoints]
%                       "mindboggle" scout from Brainstorm is used
%       roi_names    = cell array of ROI names 
% 
% DISCLAIMER: This function is inspired from Brainstorm software: 
%               https://neuroimage.usc.edu/brainstorm
%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original authors: Jun 2019; Arun Sasidharan, ABRL
% 
% Copyright (C) 2019 ABRL, Axxonet System Technologies Pvt Ltd., Bengaluru, India
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

if ~exist('roi_function','var')
    roi_function = 'mean';
end

source_data = reshape(source_data,size(sphereModel.dipoles,2),size(sphereModel.dipoles,1),size(Souredata,2));
if length(size(source_data))>2
    n_dipoles = size(source_data,1);
else
    n_dipoles = 1;
end

%% Apply ROI
load('scout_mindboggle.mat');
n_roi = length(scout_mindboggle);
roi_data = zeros(n_dipoles,n_roi,n_time);
for scout_no = 1:n_roi
    for dipole_no = 1:n_dipoles
        switch roi_function
            case 'mean'
                roi_data(dipole_no,scout_no,:)  = ...
                    mean(source_data(dipole_no,scout_mindboggle(scout_no).Vertices,:),2);
            case 'pca'
                pc = abrl_PCA(squeeze(source_data(dipole_no,scout_mindboggle(scout_no).Vertices,:)));
                roi_data(dipole_no,scout_no,:)  = pc(1,:);
        end
    end
end
roi_names = {scout_mindboggle.Label};


%% Get the norm of 3 dipoles at each voxel
if n_dipoles == 3
    roi_data = squeeze(sqrt(roi_data(1,:,:).^2 + roi_data(2,:,:).^2 + roi_data(3,:,:).^2));
else
    roi_data = squeeze(abs(roi_data));
end

%% DC correction
roi_data     = bsxfun(@minus,roi_data,mean(roi_data,2));