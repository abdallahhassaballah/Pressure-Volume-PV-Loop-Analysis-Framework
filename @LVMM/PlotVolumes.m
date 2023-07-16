function ax = PlotVolumes(L,varargin)
% Plot endocardial and epicardial volumes in one axis with additional information.
%
%   h = L.PlotVolumes;
%   h = L.PlotVolumes(...);
%
% Output:
%   - h contains two axis handles to hold epicardial & endocardial volume plots.
%
% Optional arguments:
%   - 'fig', figure. Default is a new figure.
%   - 'show_edes', boolean. Default is true.
%      Plot a vertical dashed lines marking ED & ES phases from the endo volume.
%   - 'endocol', color. Default is red.
%   - 'epicol', color. Default is blue.
%   - 'frames', numbers. Default is index numbers.
%
% Author: Avan Suinesiaputra - Centre for Advanced Imaging, Univ. of Auckland (2012)

% default arguments
opt.fig = [];
opt.endocol = 'r';
opt.epicol = 'b';
opt.show_edes = true;
opt.frames = 1:L.nframes;

% get arguments
for i=1:2:length(varargin)
    if( isfield(opt,varargin{i}) ), opt.(varargin{i}) = varargin{i+1};
    else error('Unknown option is found.'); end
end

% check figure
if( isempty(opt.fig) ), opt.fig = figure; end
if( ~ishandle(opt.fig) && ~strcmp(get(opt.fig,'Type'),'figure') )
    error('The argument ''ax'' must be a figure handle.'); 
end

% calculate volumes
V = L.ComputeVolumes;

% plotting
[ax,h1,h2] = plotyy(opt.frames,V.endo,opt.frames,V.epi);

% set properties
set(ax(1),'Parent',opt.fig,'YColor',opt.endocol);
set(ax(2),'Parent',opt.fig,'YColor',opt.epicol);

set(h1,'Color',opt.endocol);
set(h2,'Color',opt.epicol);

ylabel(ax(1),'Endocardial volumes (ml)');
ylabel(ax(2),'Epicardial volumes (ml)');

xlabel('Frames');
title(L.name,'Interpreter','none');

% set ED & ES lines
[~,edi] = max(V.endo); edi = opt.frames(edi);
[~,esi] = min(V.endo); esi = opt.frames(esi);

if( opt.show_edes )
    yl = ylim;
    hold on;
    plot([edi edi],yl,'--k');
    plot([esi esi],yl,'--k');
end
