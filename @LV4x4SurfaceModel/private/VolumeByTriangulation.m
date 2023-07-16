function V = VolumeByTriangulation(Faces,Vertices)
% Compute volume approximated by triangulation
%
%   V = VolumeByTriangulation(Faces,Vertices)
%
% Volume of a patch is approximated by:
%   V = | det([
%               FV1 FV1 FV1 1
%               FV2 FV2 FV2 1
%               FV3 FV3 FV3 1
%               C1  C2  C3  1
%                              ]) |
% where [C1,C2,C3] is the mean point of the patch [FV1,FV2,FV3].
% So we can approximate the volume of the heart by sum of volumes of each
% patch to the mean point. We assume that there is no hole in a heart.
%
% Author: Avan Suinesiaputra - Centre for Advanced Imaging, Univ. of Auckland (2012)

cog = mean(Vertices);

% calculate volumes

V = zeros(size(Faces,1),1);
for i=1:size(Faces,1)
    V(i) = abs(det([[Vertices(Faces(i,:),:); cog] ones(4,1)]));
end
V = V ./ 6;

% total volume in ml
V = sum(V) / 1000;