function P = GetIntersectionWithDICOMImage(L,dimg,varargin)
% Get the intersection contour points between a model with a DICOM image.
%
%   P = L.GetIntersectionWithDICOMImage(dimg,'opt1',val1,...);
%
% Input:  - dimg is a DICOM info or the DICOM filename.
%           DO NOT GIVE DICOM IMAGE MATRIX !!
% Output: - P is the intersection points or empty if there isn't
%
% Options:
%  - 'surfaces', [surface indices]
%    Default all
%  - 'image_shift', 3x1
%    Add manual image shifting.
%
%  These options below replace DICOM 3D information manually
%  - 'image_position', 3x1. Default is dimg.ImagePositionPatient
%  - 'image_orientation', 6x1. Default is dimg.ImageOrientationPatient
%  - 'pixel_spacing', 2x1. Default is dimg.PixelSpacing
%
%  These options below add offset to position and orientation
%  - 'image_position_offset', 3x1. Default is zero.
%  - 'image_orientation_offset', 6x1. Default is zero.
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

P = [];
if( isempty(L) ), return; end

if( ischar(dimg) )
    dimg = dicominfo(dimg);
end

% options
opt.surfaces = 1:L.nSurfaces;
opt.image_shift = [];
opt.image_position = dimg.ImagePositionPatient;
opt.image_orientation = dimg.ImageOrientationPatient;
opt.pixel_spacing = dimg.PixelSpacing;
opt.image_position_offset = zeros(3,1);
opt.image_orientation_offset = zeros(6,1);

% get options
for i=1:2:numel(varargin)
    if( isfield(opt,varargin{i}) ), opt.(varargin{i}) = varargin{i+1};
    else error('Unknown option'); end
end

% add offset
opt.image_position = opt.image_position + opt.image_position_offset;
opt.image_orientation = opt.image_orientation + opt.image_orientation_offset;

% get image position and the image vectors
V1 = opt.image_orientation(1:3);
V2 = opt.image_orientation(4:6);
V3 = cross(V1,V2);

% rotation matrix
R = [V1 V2 V3 zeros(3,1); 0 0 0 1];

% translation matrix
T = eye(4);
T(1:3,4) = opt.image_position;

% inverse transform
iR = inv(R);
iT = inv(T);

% set model's transformation to magnet
origT = L.T;
L.T = L.data.ModelToMagnetTransform;

% get the intersection points
Pi = L.GetIntersectionWithPlane(opt.image_position,V3);

% calculate the shifting
if( ~isempty(opt.image_shift) )
    shifting = opt.image_position(:) - opt.image_shift(:);
else
    shifting = zeros(3,1);
end

% for each Pi, apply scaling and add 1 because image starts from 1
P = cell(size(Pi));
for i=1:numel(Pi)
    if( isempty(Pi{i}) ), continue; end
    
    Q = iR * iT * [Pi{i}+ones(size(Pi{i},1),1)*shifting' ones(size(Pi{i},1),1)]';
    P{i} = 0.5+[Q(1,:)/opt.pixel_spacing(1); Q(2,:)/opt.pixel_spacing(2)]';
end

% reset transform
L.T = origT;