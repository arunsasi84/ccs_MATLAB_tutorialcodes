function [figure_id,colorbar_id] = abrl_headmapdensity(channel_values,chanlist,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% abrl_headmapdensity
%       - To generate headMap with colour-coded density map
% 
% Usage: [figure_id,colorbar_id] = abrl_headmapdensity(channel_values,chanlist)
% 
% Inputs:
%		channel_values 	= vector with channel values to plot
%       chanlist		= cell array with list of channels
%
% Optional Inputs:
%		plottype		= filling of plot
%                            'filled' (Default)
%                            'empty'
%       electrodetype   = appearance of electrodes
%                            'black' (Default)
%                            'coloured'
%                            'none'
%
% Outputs:
%       figure_id 		= figure handle
%		colorbar_id		= colorbar handle
%
% Dependancy: Uses "kernel_density" function from "econometrics" Octave package
% 
% DISCLAIMER: The code is re-written from "topoplot" function of EEGLAB, with
%             modifications and optimizations to allow plotting in Octave.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original author: Aug 1996; Colin Humphries & Scott Makeig, CNL / Salk Institute
% Modified author: Apr 2017; Arun Sasidharan, ABRL
% 
% Copyright (C) 1996 Colin Humphries & Scott Makeig, CNL / Salk Institute, USA
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
% 
% Modified on 04 May 2017; by Arun Sasidharan; 
%   - added provision for colorbar position
% Modified on 14 Jul 2017; by Arun Sasidharan; 
%   - added significant electrode color based on channel diff values
% Modified on 17 Jan 2019; by Arun Sasidharan; 
%   - removed the colorbar adjustment
% Modified on 17 Mar 2019; by Arun Sasidharan; 
%   - separated fill color and electrode fill as separate options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Default parameters
rmax = 0.5;
MAPLIMITS           = 'absmax';     % Other option 'maxmin';
BACKCOLOR           = [1 1 1];    	% Figure background
GRID_SCALE          = 150;          % grids size to plot map
CIRCGRID            = 201;          % number of angles to use in drawing circles
CONTOURNUM          = 6;            % number of contour levels to plot
AXHEADFAC           = 1.3;          % head to axes scaling factor
CCOLOR              = [0.1 0.1 0.1];% default contour color
CWIDTH              = 0.5;          % default contour linewidth
BLANKINGRINGWIDTH   = .35;          % width of the blanking ring 
HEADRINGWIDTH       = .007;         % width of the cartoon head ring
HLINEWIDTH          = 2;            % default linewidth for head, nose, ears
HEADCOLOR           = [0 0 0];      % default head outline color (black)
headrad             = 0.5;          % default head outline radius
EMARKER             = '.';          % mark electrode locations with small disks
ECOLOR              = [0 0 0];      % default electrode color (black)
EMARKERLINEWIDTH    = 1;            % default edge linewidth for emarkers

if nargin <= 2
	plottype = 'filled';
else
	plottype = varargin{1};
end

if nargin <= 3
	electrodetype = 'black';
else
	electrodetype = varargin{2};
end

% Convert to column vector
if size(channel_values,1) > 1
	channel_values = channel_values';
end

%% Get electrode coordinates
[chanlocs,chan_absent] = abrl_getChanLocs(chanlist);
Theta           = [chanlocs.theta];
Radii           = [chanlocs.radius];
Theta           = pi/180*Theta;   		% convert degrees to radians
[chanX,chanY]   = pol2cart(Theta,Radii);% transform electrode locations from polar to cartesian coordinates
x               = -chanY;               % rotate cordinates
y               = chanX;                % rotate cordinates

channel_values(chan_absent) = [];

%% Squeeze electrode arc_lengths towards the vertex to plot all inside the head cartoon
plotrad = min(1.0,max(Radii)*1.02);    % default: just outside the outermost electrode location
plotrad = max(plotrad,0.5);            % default: plot out to the 0.5 head boundary
if isempty(plotrad),plotrad = 0.5; end;
squeezefac = rmax/plotrad;

%% Open a figure and set the axis
hold on
h = gca; % uses current axes
set(h,'Xlim',[-rmax rmax]*AXHEADFAC,'Ylim',[-rmax rmax]*AXHEADFAC);% specify size of head axes in gca

if strcmp(plottype,'filled') && ~isempty(y),
	
	x    = x'*squeezefac;        
	y    = y'*squeezefac;   

	%% Make sure outermost channel will be plotted just inside rmax
	xmin = min(-rmax,min(x)); xmax = max(rmax,max(x));
	ymin = min(-rmax,min(y)); ymax = max(rmax,max(y));

	%% Spread data into a square grid
	xi = linspace(xmin,xmax,GRID_SCALE);   % x-axis description (row vector)
	yi = linspace(ymin,ymax,GRID_SCALE);   % y-axis description (row vector)

	[~,~,zi] = griddata(x,y,channel_values',xi,yi'); % interpolate data
	zi(isnan(zi)) = 0;
	zi = fliplr(zi);
	try %#ok<ALIGN>
        if isOctave==1; pkg load 'image'; end % add octave package for 'fspecial' function
    catch
    end
	Zi = conv2(zi, fspecial('disk',15),'same');


	%% Create 2-d grid of coordinates and function values, suitable for 3-d plotting
	[Xi,Yi]     = meshgrid(xi,yi);

	%% Add a mask outside the plotting circle
	mask    = (sqrt(Xi.^2 + Yi.^2) <= rmax); % mask outside the plotting circle
	Zi(mask == 0)  = NaN;                    % mask non-plotting voxels with NaNs

	%% Set map limits
	if strcmp(MAPLIMITS,'absmax')        %#ok<ALIGN>
        amax = max(max(abs(Zi)));
        amin = -amax;
    elseif strcmp(MAPLIMITS,'maxmin') || strcmp(MAPLIMITS,'minmax')
        amin = min(min(Zi));
        amax = max(max(Zi));
    else
        error('unknown ''maplimits'' value.');
    end

    %% Plot the density map
	figure_id = imagesc(Xi(1,:),Yi(:,1),Zi);
	axis xy;
    
    %% Adjust the color map and add a color bar
% 	cm      = colormap(jet);
% 	n       = size(cm,1);       % size of colormap
% 	dmap    = (amax-amin)/n;    % color step	
% 	colormap([[1,1,1]; cm]);    % add nan color to colormap	
% 	caxis([amin-dmap amax]);    % change color limits	
% 	colorbar_id = colorbar;             % place a colorbar    	
% 	ylim(colorbar_id,[amin amax])       % change Y limit to hide NaN color
    colorbar_id = colorbar;
    caxis([amin,amax]);
    
	hold on

	%% Add countour lines        
	[~, chs] = contour(Xi,Yi,Zi,CONTOURNUM,'k'); % surface handle
	set(chs,'linecolor',CCOLOR,'linewidth',CWIDTH);
    
end

%% Plot filled ring to mask jagged grid boundary
hin     = squeezefac*headrad*(1- HEADRINGWIDTH/2);  % inner head ring radius
rwidth  = BLANKINGRINGWIDTH;                        % width of blanking outer ring
rin     =  rmax*(1-rwidth/2);                       % inner ring radius
if hin>rin
	rin = hin;                                      % dont blank inside the head ring
end

circ    = linspace(0,2*pi,CIRCGRID);
rx      = sin(circ); 
ry      = cos(circ); 
ringx   = [[rx(:)' rx(1) ]*(rin+rwidth)  [rx(:)' rx(1)]*rin];
ringy   = [[ry(:)' ry(1) ]*(rin+rwidth)  [ry(:)' ry(1)]*rin];

% Paint the border with background color
patch(ringx,ringy,zeros(size(ringx)),BACKCOLOR,'edgecolor','none');

hold on

%% Plot head ring
headx = [rx(:)' rx(1)]*hin;
heady = [ry(:)' ry(1)]*hin;
ringh = plot(headx,heady);
set(ringh, 'color',HEADCOLOR,'linewidth', HLINEWIDTH); hold on

%% Plot ears and nose
base  = rmax-.0046;
basex = 0.18*rmax;	% nose width
tip   = 1.15*rmax; 
tiphw = .04*rmax;	% nose tip half width
tipr  = .01*rmax;	% round the nose tip
q     = .04;        % lengthen the ear
EarX  = [.497-.005, .510, .518, .5299, .5419, .54, .547, .532, .510, .489-.005]; % rmax = 0.5
EarY  = [q+.0555, q+.0775, q+.0783, q+.0746, q+.0555, -.0055, -.0932, -.1313, -.1384, -.1199];
sf    = headrad/plotrad;   % squeeze the model ears and nose by this factor
% plot nose
plot3([basex;tiphw;0;-tiphw;-basex]*sf,[base;tip-tipr;tip;tip-tipr;base]*sf,...
     zeros(size([basex;tiphw;0;-tiphw;-basex])),...
     'Color',HEADCOLOR,'LineWidth',HLINEWIDTH);
% plot left ear 
plot3(EarX*sf,EarY*sf,zeros(size(EarX)),'color',HEADCOLOR,'LineWidth',HLINEWIDTH)    
% plot right ear
plot3(-EarX*sf,EarY*sf,zeros(size(EarY)),'color',HEADCOLOR,'LineWidth',HLINEWIDTH)   

%% Make plot square and large
plotax = gca;
axis square                             % make plot square
axis off

pos = get(gca,'position');
set(plotax,'position',pos);

xlm = get(gca,'xlim');
set(plotax,'xlim',xlm);

ylm = get(gca,'ylim');
set(plotax,'ylim',ylm);                 % copy position and axis limits again

axis equal;
set(gca, 'xlim', [-0.525 0.525]); set(plotax, 'xlim', [-0.525 0.525]);
set(gca, 'ylim', [-0.525 0.525]); set(plotax, 'ylim', [-0.525 0.525]);

%% Add electrode markers
if length(y)>=160
    EMARKERSIZE = 3;
elseif length(y)>=128
    EMARKERSIZE = 3;
elseif length(y)>=100
    EMARKERSIZE = 3;
elseif length(y)>=80
    EMARKERSIZE = 4;
elseif length(y)>=64
    EMARKERSIZE = 5;
elseif length(y)>=48
    EMARKERSIZE = 6;
elseif length(y)>=32 
    EMARKERSIZE = 8;
elseif length(y)>=20 
    EMARKERSIZE = 12;
elseif length(y)<20 
    EMARKERSIZE = 16;    
end

if ~isempty(y),
	if strcmp(electrodetype,'coloured'),
		EMARKERSIZE = EMARKERSIZE*1.5;
		plot3(-x(channel_values>0),y(channel_values>0),zeros(size(x(channel_values>0))),...
					EMARKER,'Color','red','markersize',...
					EMARKERSIZE,'linewidth',EMARKERLINEWIDTH);
		plot3(-x(channel_values<0),y(channel_values<0),zeros(size(x(channel_values<0))),...
					EMARKER,'Color','blue','markersize',...
					EMARKERSIZE,'linewidth',EMARKERLINEWIDTH);
    elseif strcmp(electrodetype,'black'),
		plot3(-x,y,zeros(size(x)),...
					EMARKER,'Color',ECOLOR,'markersize',...
					EMARKERSIZE,'linewidth',EMARKERLINEWIDTH);
	end
end

%% Set background color to match head border
try 
    set(gcf, 'color', BACKCOLOR); 
catch 
end

hold off
axis off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%