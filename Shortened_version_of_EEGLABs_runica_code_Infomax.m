function [ic,ic_weights,ic_residualvar,sphere] = abrl_ICA_runica(EEGdata)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% abrl_ICA_runica
%   - To do Independent component decomposition 
%       using EEGLAB's runica (Infomax algorithm extended with kurtosis) 
% 
% Usage: [ic,ic_weights,ic_residualvar,sphere] = abrl_ICA_runica(EEGdata)
% 
% Inputs:
%       EEGdata     = [nchan x nsamples]
% 
% Optional Inputs (For plotting):
%       srate       = to plot component time series
%       chanlist    = cell list of channel labels (to plot headmap of components)
%       n_plotcomps = number of ICA components to plot       
% 
% Outputs:
%       ic          = ICA component time series [ncomp x nsamples]
%       ic_wt       = ICA component weights [ncomp x nchan]
%       ic_residualvar 
%                   = explained residual variance by each ICA component [ncomp x 1]
%       sphere      = sphering matrix [ncomp x nchan]
% 
% DISCLAIMER: code is written based on EEGLAB's runica.m
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original author: 1996-00 Scott Makeig, CNL / Salk Institute
%       with contributions from Tony Bell, Te-Won Lee, Tzyy-Ping Jung, 
%       Sigurd Enghoff, Michael Zibulevsky et al.
% Modified author: Jun 2019; Arun Sasidharan, ABRL,
% 		Axxonet System Technologies Pvt Ltd., Bengaluru, India
% 
% Copyright (C) 1996 Scott Makeig et al, SCCN/INC/UCSD, scott@sccn.ucsd.edu
% 
% The source of this code "runcia.m" is part of EEGLAB; http://www.eeglab.org
% Reference (please cite):
% Makeig, S., Bell, A.J., Jung, T-P and Sejnowski, T.J.,
%   "Independent component analysis of electroencephalographic data," 
%   In: D. Touretzky, M. Mozer and M. Hasselmo (Eds). Advances in Neural 
%   Information Processing Systems 8:145-151, MIT Press, Cambridge, MA (1996).
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data0       = EEGdata; % store the original data
datalength  = size(EEGdata,2);
nchans      = size(EEGdata,1);
ncomps      = size(EEGdata,1);
block       = ceil(5*log(datalength));

%% Remove mean from data
EEGdata         = bsxfun(@minus,EEGdata,mean(EEGdata,2));

%% Compute sphering matrix and use it to Decorrelate the electrode signals
sphere          = 2.0*inv(sqrtm(double(cov(EEGdata')))); %#ok<MINV>
EEGdata         = sphere*EEGdata;

%% Intialize variables
ic_weights      = eye(ncomps,nchans); % begin with the identity matrix
startweights    = ic_weights;
prevweights     = startweights; %#ok<*NASGU>
oldweights      = startweights;
onesrow         = ones(1,block);
bias            = zeros(ncomps,1);
BI              = block*eye(ncomps,ncomps);
lrate           = 0.00065/log(nchans);
maxweight       = 1e8;
annealstep      = 0.98;
degconst        = 180./pi;
wts_blowup      = 0;
step            = 0;
laststep        = 0;
maxsteps        = 512;
annealdeg       = 60;
extmomentum     = 0.5; % extra momentum for kurtosis
old_kk          = zeros(1,ncomps); % initialised kurtosis
signs           = ones(1,ncomps);    
signs(1)        = -1; % initialize signs to nsub -1, rest +1
signs           = diag(signs); % make a diagonal matrix
oldsigns        = zeros(size(signs));
signsbias       = 0.02; % bias towards super-Gaussian components
signcounts      = [];
extblocks       = 1;
urextblocks     = extblocks;    % original value, for resets
blockno         = 1; % running block counter for kurtosis interrupts

% Stopping Criteria
if ncomps > 32
    nochange    = 1E-7;
else
    nochange    = 1E-6;
end

while step < maxsteps, %%% ICA step = pass through all the data %%%%%%%%%%%
    
    %% Randomly shuffle the data's temporal order at each step
    timeperm = randperm(datalength);

    %% ICA Training = pass through each block
    for t = 1:block:datalength-block,      

        % Compute "u" and "y" matrices to derive weight change
        u   = ic_weights*double(EEGdata(:,timeperm(t:t+block-1))) + bias*onesrow;
        y   = tanh(u);     % hyperbolic tangent (Normalisation)
        
        weights_change = lrate*(BI-signs*y*u'-u*u')*ic_weights;% using "signs" from kurtosis
        bias_change    = lrate*sum((-2*y)')'; %#ok<*UDIM>
        
        % Update the weights and bias
        ic_weights = ic_weights + weights_change;        
        bias    = bias + bias_change;
        
        % Check if weight has blown out of proportion
        if max(max(abs(ic_weights))) > maxweight
            wts_blowup  = 1;
            change      = nochange;
        end
        
        %% Extended-ICA kurtosis estimation (recompute signs vector using kurtosis)
        if ~wts_blowup
            if extblocks > 0 && rem(blockno,extblocks) == 0,
                
                % Apply the weights
                partact     = ic_weights*double(EEGdata); 
                
                % Kurtosis estimate of each component
                m2          = mean(partact'.^2).^2;
                m4          = mean(partact'.^4);
                kk          = (m4./m2)-3.0;  
                if extmomentum %#ok<*BDLGI>
                    kk      = extmomentum*old_kk + (1.0-extmomentum)*kk; % use momentum
                    old_kk  = kk;
                end

                % Update the signs matrix baed on component kurtosis 
                signs       = diag(sign(kk + signsbias)); 

                % Check and update the sign changes based on kurtosis
                if signs == oldsigns,
                    signcount = signcount+1;
                else
                    signcount = 0;
                end
                oldsigns    = signs;

                signcounts  = [signcounts, signcount]; %#ok<*AGROW>
                if signcount >= 25,
                    extblocks = fix(extblocks * 2); % make kurt() estimation less frequent if sign is not changing
                    signcount = 0; 
                end
            end
        end
        
        blockno = blockno + 1;
        if wts_blowup
            break
        end
    end   
    
    %% Compute the overall change in weight and angle at the end of each step
    if ~wts_blowup
        oldwtchange     = ic_weights - oldweights;
        step            = step + 1;

        lrates(step)    = lrate;
        angledelta      = 0;
        delta           = reshape(oldwtchange,1,nchans*ncomps);
        change          = delta*delta';
    end


    %% If weights blow up, Re-initialize variables and Restart
    if wts_blowup || isnan(change) || isinf(change),
        step            = 0;                        % start again
        change          = nochange;
        wts_blowup      = 0;                    
        blockno         = 1;
        lrate           = lrate*0.9;                % with lower learning rate
        ic_weights         = startweights;          % and original weight matrix
        oldweights      = startweights;
        oldwtchange     = zeros(nchans,ncomps);
        delta           = zeros(1,nchans*ncomps);
        olddelta        = delta;
        extblocks       = urextblocks;
        prevweights     = startweights;
        prevwtchange    = zeros(nchans,ncomps);
        lrates          = zeros(1,maxsteps);
        bias            = zeros(ncomps,1);

        signs           = ones(1,ncomps);
        signs(1)        = -1;
        signs           = diag(signs);
        oldsigns        = zeros(size(signs));

    
    else % if weights in bounds
    %% Update the learning rate by annealing  at the end of each step
        %%%%%%%%%%%%% Print weight update information %%%%%%%%%%%%%%%%%%%%%
        if step > 2
            angledelta = acos((delta*olddelta')/sqrt(change*oldchange));
        end
        places = -floor(log10(nochange));
        fprintf('step %d - lrate %5f, wchange %8.8f, angledelta %4.1f deg\n', ...
                            step, lrate, change, degconst*angledelta);

        %%%%%%%%%%%%%%%%%%%% Save current values %%%%%%%%%%%%%%%%%%%%%%%%%%
        changes(step)   = change;
        oldweights      = ic_weights;

        %%%%%%%%%%%%%%%%%%%% Anneal learning rate %%%%%%%%%%%%%%%%%%%%%%%%%
        if degconst*angledelta > annealdeg,
            lrate       = lrate*annealstep;      % anneal learning rate
            olddelta    = delta;                 % accumulate angledelta until
            oldchange   = change;                % annealdeg is reached
        elseif step == 1                         % on first step only
            olddelta    = delta;                 % initialize
            oldchange   = change;
        end

        %%%%%%%%%%%%%%%%%%%% Apply stopping rule %%%%%%%%%%%%%%%%%%%%%%%%%%
        if step >2 && change < nochange,    % apply stopping rule
            laststep    = step;
            step        = maxsteps;         % stop when weights stabilize
        elseif change > 1000000000.0,       % if weights blow up,
            lrate       = lrate*0.8;        % keep trying
        end;                                % with a smaller learning rate
    end; % end if weights in bounds
end % End: ICA Step

%% After looping through all steps or reaching the limit of learning,   
%% Find mean variances and Sort components in descending order of explained variances
winv               = inv(ic_weights*sphere);
ic_residualvar     = sum(winv.^2).*sum((EEGdata').^2)/((nchans*datalength)-1);

[~, windex]        = sort(ic_residualvar);
windex             = windex(ncomps:-1:1);  % order large to small 
ic_residualvar     = ic_residualvar(windex);
  
lrates      = lrates(1,1:laststep); % truncate lrate history vector
laststep    = step;

ic_weights  = ic_weights(windex,:); % reorder the weight matrix
bias        = bias(windex);         % reorder bias
signs       = diag(signs);          % vectorize the signs matrix
signs       = signs(windex);        % reorder sings

wt_ica      = sphere / ic_weights;
ic          = wt_ica * data0;

end