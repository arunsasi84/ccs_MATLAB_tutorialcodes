function [leadFieldMatrix] = abrl_sLORETA(chanlist,gainMatrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% abrl_sLORETA
%   - To compute inverse solution usinf sLORETA 
% 
% Usage: [leadFieldMatrix] = abrl_sLORETA(chanlist,gainMatrix)
% 
% Inputs:
%       chanlist    = cell list of channel labels
%       gainMatrix  = EEG forward model gain matrix [nSensors x (3*nDipoles)]
% 
% Outputs:
%       leadMatrix  = EEG inverse model matrix [nDipoles x nSensors]
% 
% DISCLAIMER: This function is based on the "process_inverse_2018.m" function 
%   that is part of the Brainstorm software: https://neuroimage.usc.edu/brainstorm
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

SnrFixed        = 3;
NoiseReg        = 0.1; 
iW_noise        = 9.5491e+04;


n_chan = length(chanlist);
n_dipoleComponents = 3;
n_sources =  size(gainMatrix,2) / n_dipoleComponents;

NoiseCov = eye(n_chan) * eps;% No Noise covariance
NoiseVar = diag(NoiseCov); % Diagonal vector of the noise variances per channel

%% Apply average reference
avgRefMatrix = eye(n_chan) - (1/n_chan * ones(n_chan));
Gain         = avgRefMatrix * gainMatrix;
NoiseCov     = avgRefMatrix * NoiseCov * avgRefMatrix'; 

%% Compute the regularization and whitening matrix form Noise covariance matrix
% [It gives similar results to creating a scalar value of iW_noise = 9.5491e+04]
% [U_noise,S_noise2]  = svd(NoiseCov,'econ');
% S_noise             = sqrt(diag(S_noise2)); % singular values
% tol_noise           = length(S_noise) * eps(single(S_noise(1))); % single precision tolerance
% Rank_noise          = sum(S_noise > tol_noise);
% U_noise             = U_noise(:,1:Rank_noise);
% S_noise             = S_noise(1:Rank_noise);
% NoiseCov            = U_noise*diag(S_noise.^2)*U_noise';% rebuild noise covariance matrix with non-zero components
% RidgeFactor         = mean(diag(S_noise2)) * NoiseReg;% Ridge Regression
% NoiseCov            = NoiseCov + RidgeFactor * eye(size(NoiseCov,1));
% iW                  = U_noise*diag(1./sqrt(S_noise.^2 + RidgeFactor))*U_noise'; % inverse whitener
% iW_noise            = iW;

%% Apply the whitening matrix
Gain = iW_noise * Gain; 

%% Global decomposition & determine SNR
[U_gain,S_gain2]  = svd((Gain*Gain'));
S_gain2           = diag(S_gain2);
S_gain            = sqrt(S_gain2); % the singular values of the lead field matrix
% tol             = length(S_gain)*eps(single(S_gain(1))); % single precision tolerance
% Rank_Leadfield  = sum(S_gain > tol);

% Calculate Lambda multiplier to achieve desired SNR
SNR     = SnrFixed^2;
% Hamalainen definition of SNR in the min norm:
Lambda  = SNR/mean(S_gain.^2); % should be equalivent to Matti's average eigenvalue definition
[LambdaY,~,LambdaU] = engunits(sqrt(Lambda));% Getting standard units

fprintf('Assumed RMS of the sources is %.3g %sA-m\n',LambdaY,LambdaU);
fprintf('Assumed SNR is %.1f (%.1f dB)\n',SNR,10*log10(SNR));

%% Generate first the inversions (current dipole time series)
Kernel = Lambda * Gain' * (U_gain * diag(1./(Lambda * S_gain2 + 1)) * U_gain');

for spoint = 1:n_dipoleComponents:n_sources*n_dipoleComponents
    R   = Kernel(spoint:spoint+n_dipoleComponents-1,:) * Gain(:,spoint:spoint+n_dipoleComponents-1);
    [Ur,Sr,Vr] = svd(R); 
    Sr  = diag(Sr);
    RNK = sum(Sr > (length(Sr) * eps(single(Sr(1))))); % single precision Rank
    SIR = Vr(:,1:RNK) * diag(1./sqrt(Sr(1:RNK))) * Ur(:,1:RNK)'; % square root of inverse

    Kernel(spoint:spoint+n_dipoleComponents-1,:) = SIR * Kernel(spoint:spoint+n_dipoleComponents-1,:);
end
leadFieldMatrix = Kernel * iW_noise; % apply the overall whitener

