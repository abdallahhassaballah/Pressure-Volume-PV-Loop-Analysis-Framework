function ReorderFrames(L,frames)
% Force reording of frames
%
%   L.ReorderFrames(frames);
%
% Author: Avan Suinesiaputra - University of Auckland 2012

L.focalLengths = L.focalLengths(frames);
L.lambdas = L.lambdas(:,frames);
L.mus = L.mus(:,frames);
L.thetas = L.thetas(:,frames);

if( numel(L.phase_data)==L.nframes )
    L.phase_data = L.phase_data(frames);
end
