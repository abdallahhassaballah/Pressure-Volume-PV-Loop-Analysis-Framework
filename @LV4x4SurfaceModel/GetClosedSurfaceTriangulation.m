function [FV,iTopFaces] = GetClosedSurfaceTriangulation(L,si)
% Get the enclosed surface triangle representation.
%
%   FV = L.GetClosedSurfaceTriangulation(si);
%   [FV,iTopFaces] = L.GetClosedSurfaceTriangulation(si);
%
% Input:  - si is surface index, i.e. between 1 to L.nSurfaces
%           or 'endo' or 'epi'.
% Output: - FV is a structure with the following fields:
%           FV.Faces : triangle faces that points indices in FV.Vertices
%           FV.Vertices : list of vertices for each face
%           This representation is similar with patch (see help patch).
%         - iTopFaces are face indices of the closing surface on the base.
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

FV = struct('Faces',[],'Vertices',[]);
if( isempty(L) ), return; end

% get the surface faces & vertices
SF = L.GetSurfaceFaces(si);
SV = L.surfacePoints(SF,:);

% re-index
SF = reshape(1:size(SV,1),size(SF));

% identify the top most points in SF
topidx = SF(1:L.nCirc,1);

% add new point, the center
Pc = mean(SV(topidx,:));
SV = [SV; Pc];
iPc = size(SV,1);

% create triangle patches on top
TF = [ones(numel(topidx),1)*iPc topidx [topidx(2:end); topidx(1)]];

% combine SF & TF
FV.Faces = [SF; TF];
FV.Vertices = SV;

% the top faces
iTopFaces = (size(SF,1)+1):size(FV.Faces,1);