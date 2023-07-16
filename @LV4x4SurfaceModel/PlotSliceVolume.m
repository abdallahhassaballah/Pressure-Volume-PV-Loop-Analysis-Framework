function h = PlotSliceVolume(L,V,ax,pos,varargin)
% Plot slice image from volume and the intersection points of an
% LV4x4SurfaceModel object. The intersection plane must be parallel with
% one of the x, y or z axis.
%
%   h = PlotSliceVolume(L,V,ax,pos,...);
%
% Input:
%   - L is an LV4x4SurfaceModel object.
%   - V is a 3D volume.
%   - ax is either 'x', 'y', or 'z', which means to which axis the plane is
%     parallel to.
%   - pos is the position of the intersection plane.
%
% Available options:
%   - 'axes', axes_handle
%     Default is gca
%   - 'origin', [x0 y0 z0]
%     Default is the origin from L.BoundingBox with default options.
%   - 'res', [Xr Yr Zr]
%     Resolution of the volume V. Default is [1 1 1].
%   - 'color', cell_array
%     Color for intersection lines for each surface.
%     Default = {'r', 'b', 'g', 'y', 'c', 'm'}
%   - 'int', true | false.
%     Set interactive selection of volume indices. Default is true.
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

% defaul option
opt.axes = [];
opt.origin = [];
opt.res = [1 1 1];
opt.color = {'r', 'b', 'g', 'y', 'c', 'm'};
opt.int = true;

% get options
for i=1:2:length(varargin)
    if( isfield(opt,varargin{i}) ), opt.(varargin{i}) = varargin{i+1};
    else error('Unhandled option.'); end
end

% check axes handle
if( isempty(opt.axes) ), opt.axes = gca; end
if( ~ishandle(opt.axes) ), error('Axes is not an axes handle.'); end

% check L
if( ~isa(L,'LV4x4SurfaceModel') ), error('The first argument must be an LV4x4SurfaceModel'); end

% check origin
if( isempty(opt.origin) )
    opt.origin = min(L.GetBoundingBox);
end

[P,img,pi] = GetIntersection(L,V,pos,ax,opt.origin,opt.res);

hi = guidata(gcf);
hi.fig = gcf;
hi.pos = pos;
hi.view = ax;
hi.L = L;
hi.V = V;
hi.origin = opt.origin;
hi.res = opt.res;
hi.color = opt.color;
guidata(hi.fig,hi);

% plot img
cla reset;
PlotIntersection(hi);
RefreshAxesProperty(gca,pos,ax);
colormap gray;

% set interactive
if( opt.int )
    set(hi.fig,'KeyPressFcn',{@InteractiveKeyPressFcn});
end

end

% ---- GetIntersection
function [P,img,pi] = GetIntersection(L,V,pos,view,origin,res)
    Vsz = size(V);
    
    if( strcmpi(view,'x') )
        
        N0 = [1 0 0];
        P = L.GetIntersectionWithPlane([pos 0 0],N0);
        for i=1:numel(P),
            if( isempty(P{i}) ), continue; end
            P{i} = ((P{i}(:,2:3) - ones(size(P{i},1),1)*origin(2:3)) ./ ...
                (ones(size(P{i},1),1)*res(2:3))) + 1;
        end

        pi = round((pos - origin(1)) / res(1) + 1);
        if( pi<1 || pi>Vsz(1) ), error('Position is out of the volume bound.'); end
        img = squeeze(V(pi,:,:));
        
    elseif( strcmpi(view,'y') )

        N0 = [0 1 0];
        P = L.GetIntersectionWithPlane([0 pos 0],N0);
        for i=1:numel(P),
            if( isempty(P{i}) ), continue; end
            P{i} = ((P{i}(:,[1 3]) - ones(size(P{i},1),1)*origin([1 3])) ./ ...
                (ones(size(P{i},1),1)*res([1 3]))) + 1;
        end

        pi = round((pos - origin(2)) / res(2) + 1);
        if( pi<1 || pi>Vsz(2) ), error('Position is out of the volume bound.'); end
        img = squeeze(V(:,pi,:));

    elseif( strcmpi(view,'z') )

        N0 = [0 0 1];
        P = L.GetIntersectionWithPlane([0 0 pos],N0);
        for i=1:numel(P),
            if( isempty(P{i}) ), continue; end
            P{i} = ((P{i}(:,1:2) - ones(size(P{i},1),1)*origin(1:2)) ./ ...
                (ones(size(P{i},1),1)*res(1:2))) + 1;
        end

        pi = round((pos - origin(3)) / res(3) + 1);
        if( pi<1 || pi>Vsz(3) ), error('Position is out of the volume bound.'); end
        img = squeeze(V(:,:,pi));

    else
        error('Plane must be parallel to one of the x, y or z axis. See the help.');
    end
end

% ---- Refresh axes property
function RefreshAxesProperty(ax,pos,view)
    % intersection plane
    if( strcmpi(view,'x') )

        labels = {'Y','Z'};
        titlestr = sprintf('Projection at x = %.2f',pos);

    elseif( strcmpi(view,'y') )

        labels = {'X','Z'};
        titlestr = sprintf('Projection at y = %.2f',pos);

    elseif( strcmpi(view,'z') )

        labels = {'X','Y'};
        titlestr = sprintf('Projection at z = %.2f',pos);

    else
        error('Plane must be parallel to one of the x, y or z axis. See the help.');
    end
    
    xlabel(ax,labels{2});
    ylabel(ax,labels{1});
    title(ax,titlestr);
    
end


% ---- KeyPressFcn
function InteractiveKeyPressFcn(hobj,evnt)

    h = guidata(hobj);
    dim = find(strcmpi({'x','y','z'},h.view));

    if( lower(evnt.Character) == 'x' )
        h.view = 'x';
    elseif( lower(evnt.Character) == 'y' )
        h.view = 'y';
    elseif( lower(evnt.Character) == 'z' )
        h.view = 'z';
    elseif( evnt.Character==29 ) % next
        newpos = h.pos+h.res(dim);
        pi = round((newpos - h.origin(dim)) / h.res(dim) + 1);
        Vsz = size(h.V);
        if( pi<1 || pi>Vsz(dim) ), return; end
        h.pos = newpos;
    elseif( evnt.Character==28 ) % prev
        newpos = h.pos-h.res(dim);
        pi = round((newpos - h.origin(dim)) / h.res(dim) + 1);
        Vsz = size(h.V);
        if( pi<1 || pi>Vsz(dim) ), return; end
        h.pos = newpos;
    elseif( evnt.Character=='0' ) % reset
        h.pos = 0;
    else
        return;
    end
    
    guidata(hobj,h);
    PlotIntersection(h);
    
end

% --- just plot the intersection
function PlotIntersection(gdata)
    [P,img,pi] = GetIntersection(gdata.L,gdata.V,gdata.pos,...
        gdata.view,gdata.origin,gdata.res);
    
    hold off;

    % plot img
    imagesc(img); axis image;
    hold on;

    % plot P
    for i=1:numel(P)
        if( isempty(P{i}) ), continue; end
        h.lines(i) = plot(P{i}(:,2),P{i}(:,1),'LineStyle','-','Marker','.','Color',gdata.color{i});
    end
    
    RefreshAxesProperty(gca,gdata.pos,gdata.view);
end
