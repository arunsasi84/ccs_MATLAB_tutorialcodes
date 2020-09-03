function [fig_ID] = abrl_plotEEG_static(EEGdata,srate,chanlist,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% abrl_plotEEG_static
%   - To plot EEG data as multi-line waveform snapshot
% 
% Usage:    abrl_plotEEG_static(EEGdata,srate,chanlist,scale,fig_ID,times,cursor)
%      OR   h = abrl_plotEEG_static(EEGdata,srate,chanlist)
% 
% Inputs:
%   EEGdata  = [nchan x nsamples]
%   srate    = sampling rate
%   chanlist = string array of channel list
% 
% Optional Input
%   scale    = scalar value that determins the spread of the waveforms
%               [Default = Auto]
%   fig_ID   = figure ID
%   times    = vector with time values to plot
%   cursor   = sample point at which a vertical cursor need to be drawn
%              0 -> no cursor [Default] 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original author: Oct 2018; Arun Sasidharan, ABRL
% 
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
% 
% Modified on: 08 Nov 2018; by Arun Sasidharan
%   - Optimised the code and added channel labels as yticks
% Modified on: 19 Dec 2018; by Arun Sasidharan
%   - Added auto scale as default and displays the scaling factor used
% Modified on: 28 Jun 2019; by Arun Sasidharan
%   - Bug fix for optional arguments
% Modified on: 29 Jun 2019; by Arun Sasidharan
%   - Added optional arguments - time and cursor
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nbchan = size(EEGdata,1);
nbpts = size(EEGdata,2);

%% Add numbers as channel labels when not specified
if ~exist('chanlist','var')
    chanlist = strsplit(num2str([1:nbchan]))';
elseif isempty(chanlist)
    chanlist = strsplit(num2str([1:nbchan]))';
end

if nargin > 3 && ~isempty(varargin{1}) && isnumeric(varargin{1})
    scale = varargin{1};
else
    % Find scaling factor from data range
    [EEGmin, EEGmax] = deal(min(min(EEGdata')),max(max(EEGdata')));
    scale = abs(EEGmin-EEGmax);
end

if nargin > 4 && ~isempty(varargin{2}) && isnumeric(varargin{2})
    fig_ID = varargin{2};
else
    fig_ID = figure;
end

if nargin > 5 && ~isempty(varargin{3}) && numel(varargin{3}) == nbpts
    times = varargin{3};
else
    times = (1:nbpts)/srate;
end

if nargin > 6 && ~isempty(varargin{4}) && isnumeric(varargin{4})
    cursor = varargin{4};
else
    cursor = 0;
end

if nbchan ~= length(chanlist)
    error('Channel list do not match the EEG data')
end

if size(chanlist,1)<size(chanlist,2)% Make sure channel labels are oriented right
    chanlist = chanlist';
end


chanlabel_pos = zeros(1,nbchan+2);
chanlabel_pos(1) = (nbchan+1)*scale;
for chan = 1:nbchan
    yfactor = ones(size(EEGdata(chan,:)))*(nbchan-chan+1)*scale;
    chanlabel_pos(chan+1) = (nbchan-chan+1)*scale;
    plot(times,EEGdata(chan,:)+yfactor);
    hold on
    axis on
end
set(gca,'ytick',fliplr(chanlabel_pos));
set(gca,'yticklabel',flipud([' ';chanlist;' ']));
xlabel('Time (s)');
xlim([min(times),max(times)]);
ylim([min(chanlabel_pos),max(chanlabel_pos)]);
title(sprintf('SR: %iHz \t Scaling: %f',srate,scale));

if cursor ~= 0
    cursor_time = times(cursor);
    plot([cursor_time,cursor_time],[min(chanlabel_pos),max(chanlabel_pos)],'-.r')
end