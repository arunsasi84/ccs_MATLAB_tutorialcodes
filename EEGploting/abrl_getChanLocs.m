function [chanlocs,chan_absent] = abrl_getChanLocs(chanlist)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% abrl_getChanLocs
%   - To get channel locations compatible with EEGLAB 
% 
% Usage: [chan_loc,chan_absent] = abrl_getChanLocs(chan_list)
% 
% Inputs:
%       chanlist = cell list of channel labels
% 
% Outputs:
%       chanlocs    = channel location in EEGLAB format
%       chan_absent = indices of channels that did not have a valid channel location       
% 
% NOTE: the channel locations are read from "ChanLoc64.mat" file that has a
% list of all standard 64 EEG channels in 10-10 system and also 4 channels
% that are part of the older 10-20 system (i.e., T3, T4, T5, T6). Even A1,
% A2, M1 and M2 have locations and will not be marked as "chan_absent".
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original author: Apr 2017; Arun Sasidharan, ABRL
% 
% Copyright (C) 2017 ABRL, Axxonet System Technologies Pvt Ltd., Bengaluru, India
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

%% Load the default channel location file (This file should be in the path)
default_chanlocfile = 'ChanLoc64.mat';
load(which(default_chanlocfile));

%% Read channel location for each channel and store into a structure
noof_chan = length(chanlist);
chanlocs = defaultChanLocs(1:noof_chan); %#ok<NODEF>
chan_absent = [];
for chan_no = 1:noof_chan
    chan_index = find(cellfun(@(x) strncmpi(x,chanlist{chan_no},...
        length(chanlist{chan_no})), {defaultChanLocs(:).labels}), 1, 'first');
    if ~isempty(chan_index)
        chanlocs(chan_no) = defaultChanLocs(chan_index);
    else
        chan_absent = [chan_absent,chan_no]; %#ok<AGROW>
				
				% Create empty chanloc info
				chanlocs(chan_no).labels = chanlist{chan_no};
				chanlocs(chan_no).type = '';
				chanlocs(chan_no).theta = [];
				chanlocs(chan_no).radius = [];
				chanlocs(chan_no).X = [];
				chanlocs(chan_no).Y = [];
				chanlocs(chan_no).Z = [];
				chanlocs(chan_no).sph_theta = [];
				chanlocs(chan_no).sph_phi = [];
				chanlocs(chan_no).sph_radius = [];
				chanlocs(chan_no).urchan = chan_no;
				chanlocs(chan_no).ref = '';		
    end
end
