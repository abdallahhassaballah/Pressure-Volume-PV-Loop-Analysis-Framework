function [M,N] = CreateSamplingMap(L)
% Given the sampling elements (L.E), this function outputs the surface map M,
% where M is nMus x nThetas. The rows start from base (row=1) to apex
% (row=nMus). The columns starts from 0 (first) to 2*pi (last).
%
% The other output N is the nodal element indices: 4x16

% prepare N
N = zeros(16,4);

% get the nSamples
nSamples = L.nSamples;

% nPS is the number of total sample points per dimension
nPS = 4*nSamples - 3;

% create the first element row: 1, 2, 3, 4
M = reshape(1:(nSamples(1)*nPS(2)),nSamples(1),[]);
M(:,1:nSamples(2)) = fliplr(M(:,1:nSamples(2)));

j=(0:2)*(nSamples(2)-1)+1+nSamples(2);
k=j+(nSamples(2)-2);
for i=1:3
    M(:,j(i):k(i)) = fliplr(M(:,j(i):k(i)));
end
M = flipud(M);

% the nodal elements for element 1-4
N(1,:) = [1 nSamples(1) prod(nSamples)-nSamples(2)+1 prod(nSamples)];
N(2,1:2) = N(1,1:2) + prod(nSamples);
N(3,1:2) = N(2,1:2) + prod(nSamples) - nSamples(2);
N(4,1:2) = N(3,1:2) + prod(nSamples) - nSamples(2);
N(2:4,3:4) = N(1:3,1:2);

nI = (nSamples(2)-1)*(nSamples(1):(nSamples(1)-1):nPS(2))';

% create the last three element rows
for z=1:3
    newMap = reshape(1:(nSamples(1)-1)*nPS(2),[],nPS(2));
    newMap(:,1:nSamples(2)) = fliplr(newMap(:,1:nSamples(2)));
    for i=1:3
        newMap(:,j(i):k(i)) = fliplr(newMap(:,j(i):k(i)));
    end
    newMap = max(M(:))+flipud(newMap);
    
    % the nodal elements
    N(z*4+1:z*4+4,1) = newMap(nI);
    N(z*4+1:z*4+4,2) = N((z-1)*4+1:(z-1)*4+4,1);
    N(z*4+1,3) = newMap(nSamples(2)-1);
    N(z*4+1,4) = N((z-1)*4+1,3);
    N(z*4+2:z*4+4,3:4) = N(z*4+1:z*4+3,1:2);
    
    M = [M; newMap];
end

% transpose N
N = N';