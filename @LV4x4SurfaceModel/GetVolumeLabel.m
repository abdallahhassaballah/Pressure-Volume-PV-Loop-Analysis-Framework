function [V,org,res] = GetVolumeLabel(L,varargin)
% Get the volumetric pixel labels from an LV model.
%
%   V = L.GetVolumeLabel;
%   V = L.GetVolumeLabel('opt1',val1,'opt2',val2,...);
%
% Available options:
%   - 'res', [xres yres zres].
%     Define the increment steps for all axes. Default is [1 1 1].
%   - 'bb_margin', xyz_margin | [xmargin ymargin zmargin]
%     Define the margin to create the bounding box (see GetBoundingBox method).
%     Default is [5 5 5].
%   - 'limit', [xmin ymin zmin; xmax ymax zmax].
%     Define the limit. Default is to use GetBoundingBox method.
%     Note that specifying this option will ignore the 'bb_margin' option.
%
%   - 'output', 'struct' | 'volume'
%     Define how the output will be. Default is 'struct'.
%
%     If 'output' is 'struct', then V is structure with the following fields:
%       V.bounding_box = [xmin ymin zmin; xmax ymax zmax]. The bounding box.
%       V.resolution = [xres yres zres]. The resolution.
%       V.myocardium = [xi yi zi]. All voxels that are inside the myocardium.
%       V.cavity = [xi yi zi]. All voxels that are inside the cavity.
%
%     If 'output' is 'volume', then V is a 3D matrix containing the following labels:
%       0 = background
%       1 = myocardium
%       2 = LV cavity
%     The size of V is Nx x Ny x Nz.
%
% Notes:
% 1. This method only works if there are 2 surfaces.
% 2. If you use volume as the output type, then it becomes a 3D matrix.
%    There is no more information about it's origin and increment.
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

if( 2~=L.nSurfaces ), error('GetVolumeLabel only works for 2 surfaces model.'); end

% default options
opt.res = [1 1 1];
opt.bb_margin = 5;
opt.output = 'struct';
opt.limit = [];

% get options
for i=1:2:length(varargin)
    if( isfield(opt,varargin{i}) ), opt.(varargin{i}) = varargin{i+1};
    else error('Unhandled option.'); end
end

if( isempty(find(strcmpi(opt.output,{'struct','volume'}), 1)) )
    error('Unknown output type %s',opt.output);
end

% get the bounding box
if( isempty(opt.limit) )
    BB = L.GetBoundingBox('margin',opt.bb_margin);
    opt.limit = [min(BB); max(BB)];
end

% get cavity
Pcav = L.GetPointsInsideSurface(1,'limit',opt.limit,'res',opt.res);

% get myo
Pmyo = L.GetPointsInsideSurface(L.nSurfaces,'limit',opt.limit,...
    'res',opt.res,'exclude',Pcav);

% still need to remove some speckles above the base
Bp = L.GetBasalPoints;
Bp = Bp(:,:,1);
Px = Pmyo(Pmyo(:,1)<max(Bp(:,1)),:);
in = inpolygon(Px(:,2),Px(:,3),Bp(:,2),Bp(:,3));
if( any(in) )
    Pmyo = setdiff(Pmyo,Px(in,:),'rows');
end

% construct output
if( strcmpi(opt.output,'struct') )
    
    V.bounding_box = opt.limit;
    V.resolution = opt.res;
    V.myocardium = Pmyo;
    V.cavity = Pcav;
    
elseif( strcmp(opt.output,'volume') )
    V = zeros(...
        numel(opt.limit(1,1):opt.res(1):opt.limit(2,1)),...
        numel(opt.limit(1,2):opt.res(2):opt.limit(2,2)),...
        numel(opt.limit(1,3):opt.res(3):opt.limit(2,3)));
    
    P0 = opt.limit(1,:);
    Pmyo = Pmyo + ones(size(Pmyo,1),1)*(-P0+1);
    iMyo = sub2ind(size(V),Pmyo(:,1),Pmyo(:,2),Pmyo(:,3));
    V(iMyo) = 1;

    Pcav = Pcav + ones(size(Pcav,1),1)*(-P0+1);
    iCav = sub2ind(size(V),Pcav(:,1),Pcav(:,2),Pcav(:,3));
    V(iCav) = 2;
    
end

% set optional output
org = min(opt.limit);
res = opt.res;