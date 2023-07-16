function PlotAHASegmentsXi(L,si,varargin)
% Plot segments in the element space.
%
%   L.PlotAHASegmentsXi(si);
%   L.PlotAHASegmentsXi(si,'opt1',val1,...);
%
% Input:  - si is either 'epi' or 'endo'.
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

if( isempty(L) ), return; end

% options
opt.color_apicaltip = [1 0.69 0.39];
opt.color_apex = 'c';
opt.color_mid = 'm';
opt.color_base = [0 0.5 0];
opt.border_style = {'LineStyle','-','Marker','o'};
opt.show_points = true;
opt.point_style = {'Marker','x','LineStyle','none'};
opt.show_rvgp = true;
opt.rvgp_style = {'Marker', 's', 'Color', 'k', 'MarkerFaceColor', 'y', 'LineStyle', 'none'};
opt.aha_region_opts = {};

% get optional arguments
for i=1:2:length(varargin)
    if( isfield(opt,varargin{i}) ), opt.(varargin{i}) = varargin{i+1};
    else error('Unknown option.'); 
    end
end

% color
gridcolors = {'-b.','-r.'};

% create new figure
figure;
L.PlotElementGrid(si,gridcolors{strcmp(si,{'endo','epi'})});
hold on;

% get the regions
[R,mus,thetas] = L.GetAHAPoints(opt.aha_region_opts{:});
[~,lvl_idx] = L.GetAHASliceLevels;
lvl_idx = lvl_idx.(si);

% plot mu the boundaries
xi_thetas = L.E(L.map(1,:),1);
if( strcmpi(si,'epi') )
    tm = L.Xi2ThetaMu(si,[xi_thetas(:) mus.epi(1)*ones(numel(xi_thetas),1)],L.E(L.map(end,:),4));
    plot(tm(:,1),tm(:,2),opt.border_style{:},'Color',opt.color_apicaltip);
    
    amb = mus.epi(2:end);
else
    amb = mus.endo;
end

bs = {opt.color_apex,opt.color_mid,opt.color_base};
for i=1:numel(amb)
    if( amb(i)>=3 )
        eids = L.E(L.map(end,:),4) - 12;
    elseif( amb(i)>=2 )
        eids = L.E(L.map(end,:),4) - 8;
    elseif( amb(i)>=1 )
        eids = L.E(L.map(end,:),4) - 4;
    else
        eids = L.E(L.map(end,:),4);
    end
    
    mu = amb(i) - floor(amb(i) - 0.00001);
    tm = L.Xi2ThetaMu(si,[xi_thetas(:) mu*ones(numel(xi_thetas),1)],eids);
    plot(tm(:,1),tm(:,2),opt.border_style{:},'Color',bs{i});
end


% now the points
c = [repmat({'b','r'},1,3) repmat({'b','r'},1,3) {'b','r','b','r'}];
if( opt.show_points )
    
    if( strcmpi(si,'epi') )
        plot(L.surfacePSPoints(R{17},1),L.surfacePSPoints(R{17},2),'Color','m',opt.point_style{:});
        li = fliplr(lvl_idx(2:end));
    else
        li = fliplr(lvl_idx);
    end
    
    for i=1:6
        idx = intersect(li{1},R{i});
        plot(L.surfacePSPoints(idx,1),L.surfacePSPoints(idx,2),'Color',c{i},opt.point_style{:});
    end
    
    for i=7:12
        idx = intersect(li{2},R{i});
        plot(L.surfacePSPoints(idx,1),L.surfacePSPoints(idx,2),'Color',c{i},opt.point_style{:});
    end
    
    for i=13:16
        idx = intersect(li{3},R{i});
        plot(L.surfacePSPoints(idx,1),L.surfacePSPoints(idx,2),'Color',c{i},opt.point_style{:});
    end
    
end

% rvgp
if( opt.show_rvgp && isfield(L.data,'GP_RV_INSERTION') )
    L.PlotElementXi(si,L.data.GP_RV_INSERTION,opt.rvgp_style{:});
end

