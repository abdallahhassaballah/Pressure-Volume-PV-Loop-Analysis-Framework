function P = GetPointsInsideSurface(L,si,varargin)
% Get a set of points that are inside of an LV surface
%
%   P = L.GetPointsInsideSurface(si);
%   P = L.GetPointsInsideSurface(si,'opt1',val1,...);
%
% Input:  - si is the surface index, i.e. 1..L.nSurfaces or 'endo' or 'epi'
% Output: - P are 3D points that are inside surface si.
%
% Available options:
%   - 'res', [xres yres zres]
%     Resolution of the point P. Default is [1 1 1].
%   - 'exclude', points
%     Exclude a set of points.
%   - 'limit', [xmin ymin zmin; xmax ymax zmax].
%     Define the limit. Default is to use GetBoundingBox method.
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

P = [];
if( isempty(L) ), return; end

if( ischar(si) )
    if( strcmpi(si,'endo') ), si = 1;
    elseif( strcmpi(si,'epi') ), si = L.nSurfaces;
    else error('Unknown surface %s',si); 
    end
else
    if( si<1 || si>L.nSurfaces ), error('Surface index is out of range.'); end
end

% default options
opt.res = [1 1 1];
opt.exclude = [];
opt.limit = [];
opt.margin = 5;

% get options
for i=1:2:length(varargin)
    if( isfield(opt,varargin{i}) ), opt.(varargin{i}) = varargin{i+1};
    else error('Unhandled option.'); end
end

% get the surface faces
Faces = L.GetSurfaceFaces(si);
Vertices = L.surfacePoints(Faces,:);
Faces = reshape(1:size(Vertices,1),size(Faces));

% create triangulation, normal faces and circumcenters
Tri = TriRep(Faces,Vertices(:,1),Vertices(:,2),Vertices(:,3));
NF = Tri.faceNormals;
CF = Tri.circumcenters;

% The method uses short axes slices, because LV always stands straight on
% x-axis. So for each slice, intersection faces are calculated and points
% inside these intersection faces are computed.

% To speed up, select only possible short-axes slices.
if( isempty(opt.limit) )
    BB = L.GetBoundingBox('margin',opt.margin);
    opt.limit = [min(BB); max(BB)];
end
Xi = opt.limit(1,1):opt.res(1):opt.limit(2,1);
Xi = Xi(Xi>=min(Vertices(:,1)) & Xi<=max(Vertices(:,1)));

% iterate over Xi
for i=1:numel(Xi)
    
    % find the intersection triangles
    [~,Fi] = L.GetIntersectionWithPlane([Xi(i),0,0],[1 0 0],'surfaces',si);
    if( isempty(Fi{1}) ), continue; end

    % get corresponding normal vectors and circumcenters
    NXi = NF(Fi{1},:);
    CXi = CF(Fi{1},:);
    
    % remove NaN circumcenters
    iNum = find(~isnan(CXi(:,1)) & ~isnan(CXi(:,2)) & ~isnan(CXi(:,3)));
    NXi = NXi(iNum,:);
    CXi = CXi(iNum,:);

    % instead of the whole BB(Yi,Zi), create only small bounding box around
    % the intersection triangles, padding 1 pixel
    bbs = [floor(min(CXi(:,2:3)))-1; ceil(max(CXi(:,2:3)))+1];
    [Yi,Zi] = meshgrid(max(bbs(1,1),opt.limit(1,2)):opt.res(2):min(bbs(2,1),opt.limit(2,2)),...
        max(bbs(1,2),opt.limit(1,3)):opt.res(3):min(bbs(2,2),opt.limit(2,3)));
    Pi = [Xi(i)*ones(numel(Yi),1) Yi(:) Zi(:)];
    
    % exclude points
    if( ~isempty(opt.exclude) )
        Pi = setdiff(Pi,opt.exclude,'rows');
    end

    % find the minimum distance from each Pi to the intersection triangles
    PCi = repmat(Pi,[1 1 numel(iNum)]) - permute(repmat(CXi,[1 1 size(Pi,1)]),[3 2 1]);
    [~,imin] = min(squeeze(sum(PCi.^2,2)),[],2);

    % dot product PCi with the closest normal vector
    dotPN = zeros(numel(imin),1);
    for j=1:numel(imin)
        dotPN(j) = PCi(j,:,imin(j)) * NXi(imin(j),:)';
    end

    % points inside the surface are points with negative dot product value
    % with the closest normal vector of the intersection triangles
    if( any(dotPN>=0) )
        P = [P; Pi(dotPN>=0,:)]; 
    end
    
end

% there is still a problem in the base
% we need to correct points in the base area that should be inside the
% enclosed surface triangles

% get the basal points
Bp = reshape(permute(L.GetBasalPoints,[2 1 3]),3,[])';

% locate points inside the surface around the basal plane
Pb = P(P(:,1)<=max(Bp(:,1))+1,:);

% get the closed surface triangulation
[Fbe,iBaseFaces] = L.GetClosedSurfaceTriangulation(si);
Fbe.Vertices = Fbe.Vertices(Fbe.Faces(iBaseFaces,:),:);
Fbe.Faces = reshape(1:size(Fbe.Vertices,1),size(Fbe.Faces(iBaseFaces,:)));

% create the corresponding TriRep, normal vectors and circumcenters
Tbe = TriRep(Fbe.Faces,Fbe.Vertices(:,1),Fbe.Vertices(:,2),Fbe.Vertices(:,3));
Nbe = Tbe.faceNormals;
Cbe = Tbe.circumcenters;

% remove NaN circumcenters
iNum = find(~isnan(Cbe(:,1)) & ~isnan(Cbe(:,2)) & ~isnan(Cbe(:,3)));
Cbe = Cbe(iNum,:);
Nbe = Nbe(iNum,:);

% get the closest triangles
PCi = repmat(Pb,[1 1 numel(iNum)]) - permute(repmat(Cbe,[1 1 size(Pb,1)]),[3 2 1]);
[~,imin] = min(squeeze(sum(PCi.^2,2)),[],2);

% dot product with normal vectors
dotPN = zeros(numel(imin),1);
for j=1:numel(imin)
    dotPN(j) = PCi(j,:,imin(j)) * Nbe(imin(j),:)';
end

% find the negative to remove that from the points inside the surface
if( any(dotPN<0) ), 
    P = setdiff(P,Pb(dotPN<0,:),'rows');
end
