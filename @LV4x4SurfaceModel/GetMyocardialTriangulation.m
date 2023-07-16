function FV = GetMyocardialTriangulation(L)
% Get surface triangulation of the myocardium, enclosed by endo and
% epicardial surfaces.
%
%   FV = L.GetMyocardialTriangulation;
%
% Output is a structure with fields:
%   FV.Faces : face definition that points to vertices
%   FV.Vertices : 3D coordinate of vertices
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

FV = struct('Faces',[],'Vertices',[]);
if( isempty(L) ), return; end

if( L.nSurfaces<2 ), error('The number of surfaces must be more than 2.'); end

% get the epicardial & endocardial faces
FVendo = L.GetSurfaceFaces(1);
FVepi = L.GetSurfaceFaces(L.nSurfaces);

% calculate the top faces
FVtop = [FVendo(1:L.nCirc,1) FVendo([2:L.nCirc 1],1) FVepi([2:L.nCirc 1],1)];
FVtop = [FVtop; FVepi(1:L.nCirc,1) FVendo(1:L.nCirc,1) FVepi([2:L.nCirc 1],1)];

% endo faces must have normals forward
FVendo = [FVendo(:,2) FVendo(:,1) FVendo(:,3)];

% combine
FV.Faces = [FVendo; FVepi; FVtop];
FV.Vertices = L.surfacePoints;