function [P,Fidx] = GetIntersectionWithPlane(L,P0,N0,varargin)
% Calculate points that intersection a plane with the LV model
%
%   P = L.GetIntersectiontWithPlane(P0,N0,'opt1',val1,...);
%   [P,Fidx] = L.GetIntersectionWithPlane(P0,N0,'opt1',val1,...);
%
% The plane is defined by the normal N0 and a point P0 on the plane.
% The outputs P & Fidx are cell array where each cell defines the
% intersection for each surface (see the 'surfaces' option).
% - P{i} are Nx3 coordinate points on surface i that intersect with the
%   plane.
% - Fidx{i} are indices of the surface faces indicating triangles that
%   intersect the plane.
%
% Available options:
%   - 'surfaces', <array_of_surface_indices>
%     Default are all surfaces from 1:L.nSurfaces.
%
% Source of inspiration: Triangle-Plane Intersection
% (http://softsurfer.com/Archive/algorithm_0105/algorithm_0105.htm)
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)


% default options
opt.surfaces = 1:L.nSurfaces;

% get options
for i=1:2:length(varargin)
    if( isfield(opt,varargin{i}) ), opt.(varargin{i}) = varargin{i+1};
    else error('Unhandled option.'); end
end

% adjust P0 & N0 into a column vector
P0 = P0(:);
N0 = N0(:) ./ norm(N0(:));

if( numel(P0)~=3 || numel(N0)~=3 )
    error('Invalid number of dimension for the P0 or N0.');
end

P = {}; Fidx = {};
for si=1:numel(opt.surfaces)
    
    % get the faces
    Faces = L.GetSurfaceFaces(opt.surfaces(si));

    % --- find triangles that intersect with plane: (P0,N0)

    % calculate sign distance of each vertices
    dist = N0' * (L.surfacePoints - ones(size(L.surfacePoints,1),1)*P0')';
    sgnDist = sign(dist(Faces));

    % find triangles that have points on both plane's sides
    Fidx{si} = find(any(sgnDist>0,2) & any(sgnDist<0,2));
    
    if( isempty(Fidx{si}) )
        P{si} = [];
        continue;
    end
    
    % -- find the intersection lines

    % find segments for each intersected triangles that intersects the plane
    % see http://softsurfer.com/Archive/algorithm_0104/algorithm_0104B.htm

    % pivot points
    iPos = Fidx{si}(sum(sgnDist(Fidx{si},:)>0,2)==1); 
    iNeg = Fidx{si}(sum(sgnDist(Fidx{si},:)<0,2)==1);

    p1 = [];
    u = [];
    for i=1:numel(iPos)  % triangles where only one point on positive side
        % pivot points
        p1 = [p1 L.surfacePoints(Faces(iPos(i),sgnDist(iPos(i),:)>0),:)'];

        % u vectors (2)
        u = [u L.surfacePoints(Faces(iPos(i),sgnDist(iPos(i),:)<=0),:)' - ...
            repmat(p1(:,end),1,2)];
    end 

    for i=1:numel(iNeg)  % triangles where only one point on negative side
        % pivot points
        p1 = [p1 L.surfacePoints(Faces(iNeg(i),sgnDist(iNeg(i),:)<0),:)'];

        % u vectors (2)
        u = [u L.surfacePoints(Faces(iNeg(i),sgnDist(iNeg(i),:)>=0),:)' - ...
            repmat(p1(:,end),1,2)];
    end

    % calculate the intersection point on each triangle side
    p1 = reshape(repmat(p1,2,1),3,[]);
    sI = -(N0' * (p1 - P0 * ones(1,size(p1,2)))) ./ (N0' * u);
    pts = (p1 + (repmat(sI,3,1) .* u))';
    
    % add vertices that are on the surface
    Pon = L.surfacePoints(Faces(sgnDist==0),:);
    pts = unique([pts; Pon],'rows');

    % order
    npts = size(pts,1);
    orderIdx = 1;
    while( numel(orderIdx) < npts )

        otherIdx = setdiff(1:npts,orderIdx);

        % dist of the end to the other indices
        dist = sum((ones(numel(otherIdx),1) * pts(orderIdx(end),:) - pts(otherIdx,:)).^2,2);
        [~,ccMini] = min(dist);

        % put it in the orderIdx
        orderIdx = [orderIdx otherIdx(ccMini)];

    end

    % find the end point
    pts = pts(orderIdx,:);
    i1 = 1;
    maxDist = sum((pts(1,:)-pts(end,:)).^2,2);
    for i=2:npts
        if( sum((pts(i,:)-pts(i-1,:)).^2,2) > maxDist ), i1 = i; end
    end
    pts = circshift(pts,-(i1-1));
    
    % store
    P{si} = pts;
    
end
