function [pc,pc_wt,pc_residualvar] = abrl_PCA(EEGdata,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% abrl_PCA
%   - To do Principal component decomposition using eigen value decomposition 
% 
% Usage: [pc,pc_wt,pc_residualvar] = abrl_PCA(EEGdata)
%        [pc,pc_wt,pc_residualvar] = abrl_PCA(EEGdata,srate,chanlist,3)
% 
% Inputs:
%       EEGdata = [nchan x nsamples]
% 
% Optional Inputs (For plotting):
%       srate       = to plot component time series
%       chanlist    = cell list of channel labels (to plot headmap of components)
%       n_plotcomps = number of PCA components to plot       
% 
% Outputs:
%       pc      = PCA component time series [ncomp x nsamples]
%       pc_wt   = PCA component weights [ncomp x nchan]
%       pc_residualvar 
%               = explained residual variance by each PCA component [ncomp x 1]
% 
% DISCLAIMER: code is written based on Mike X Cohen's book titled 
%   "Analyzing Neural Time Series Data"(MIT Press) Chapter 23
%   Mike X Cohen (mikexcohen@gmail.com)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original author: 2014 Mike X Cohen (mikexcohen@gmail.com)
% Modified author: Sep 2018; Arun Sasidharan, ABRL
% 
% Copyright (C) 2014 Mike X Cohen (mikexcohen@gmail.com)
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

npnts = size(EEGdata,2);

% subtract mean and compute covariance
EEGdata     = bsxfun(@minus,EEGdata,mean(EEGdata,2));
covar       = (EEGdata*EEGdata')./(npnts-1);

% principle components analysis via eigenvalue decomposition
[pc_wt,pc_residualvar] = eig(covar);

% components are listed in increasing order, and converted here to descending order for convenience
pc_residualvar  = diag(pc_residualvar);
pc_residualvar  = 100*pc_residualvar./sum(pc_residualvar); % convert to percent change
[pc_residualvar,indx] = sort(pc_residualvar,'descend');
pc_wt           = pc_wt(:,indx);

%% Compute the Principle components
pc = (pc_wt'*EEGdata);

%% Plot components if requested
if nargin > 1 && isnumeric(varargin{1})
    srate = varargin{1};
    plotoption = true;
    if nargin > 3 && isnumeric(varargin{3})
        n_plotcomps = varargin{3};
    else
        n_plotcomps = min(9,length(pc_residualvar));
    end
    if nargin > 2 && iscell(varargin{2})
        plotheadmap = true;
        chanlist = varargin{2};
    elseif nargin > 3 && iscell(varargin{3})
        plotheadmap = true;
        chanlist = varargin{3};
    else
        plotheadmap = false;
    end
else
    plotoption = false;
end

if plotoption
    fig_no = 1;
    timebins = (1:size(pc,2))/srate;
    figure;
    for i=1:n_plotcomps
        if plotheadmap
            subplot(n_plotcomps,2,fig_no);
            abrl_headmapdensity(double(pc_wt(:,i)),chanlist);
            title([ 'PC #' num2str(i) ', eigval=' num2str(pc_residualvar(i)) ]);
            fig_no = fig_no + 1;
            
            subplot(n_plotcomps,2,fig_no);
            plot(timebins,pc(i,:));
            title([ 'PC #' num2str(i) ', eigval=' num2str(pc_residualvar(i)) ]);
            xlim([min(timebins),max(timebins)]);
            xlabel('Time (s)');
            fig_no = fig_no + 1;
        else
            subplot(n_plotcomps,1,fig_no)
            plot(timebins,pc(i,:))
            title([ 'PC #' num2str(i) ', eigval=' num2str(pc_residualvar(i)) ]);
            xlim([min(timebins),max(timebins)]);
            xlabel('Time (s)');
            fig_no = fig_no + 1; 
        end
    end
end
  