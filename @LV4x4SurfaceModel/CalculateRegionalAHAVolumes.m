function [V,FV]= CalculateRegionalAHAVolumes(L,G)
% Compute regional volumes divided based on AHA segments
%
%   [V,FV] = L.CalculateRegionalAHAVolumes(G);
%
% The input G is the output from GetAHARegions.
%
% The output V contains three fields: 'endo', 'epi' and 'myo', each has 16, 17 and 17 volumes.
% The optional FV output contains volume patches.
%
% Author: Avan Suinesiaputra - Centre for Advanced Imaging, Univ. of Auckland (2012)

% remove any transformation
Ltmp = L.T;
L.T = [1 0 0; 0 1 0; 0 0 1; 0 0 0];

% get the surface patches
[F,B] = L.GetAHASurfacePatches(G);

% get tip and base
C = L.GetCircumferentialPoints([1 L.nAzimuth]);

% we'll take base from both surfaces
Cbase = mean(reshape(permute(C(:,1,:,:),[1 2 4 3]),[],3));
Capex = squeeze(mean(C(:,2,:,:),1))';

% from base to apex, Bd = endo, Bp = epi
Bd = [Cbase; squeeze(mean(B.endo,1))'; Capex(1,:)];
Bp = [Cbase; squeeze(mean(B.epi,1))'; Capex(2,:)];

% approximating for the Axes, we drop Y,Z coordinates of the borders
Ad = [Bd(:,1) zeros(4,2)];
Ap = [Bp(:,1) zeros(5,2)];

% prepare for Face Vertex structure
A.endo = permute(reshape([repmat(Ad(1:2,:),[6 1]); repmat(Ad(2:3,:),[6 1]); repmat(Ad(3:4,:),[4 1])],[2 16 3]),[1 3 2]);
A.epi = permute(reshape([repmat(Ap(1:2,:),[6 1]); repmat(Ap(2:3,:),[6 1]); repmat(Ap(3:4,:),[4 1]); Ap(4:5,:)],[2 17 3]),[1 3 2]);

% we'll use Ad & Ap to create regional volume patches

for si={'endo','epi'}
    for i=1:17
        if( i > numel(F.(si{1})) ), continue; end
        
        % append vertices at the end
        FV.(si{1})(i).Vertices = [F.(si{1}){i}.Vertices; A.(si{1})(:,:,i)];
        upper_fi = size(FV.(si{1})(i).Vertices,1)-1;
        bottom_fi = upper_fi + 1;

        % get surface faces & grids
        Fi = F.(si{1}){i}.Faces;
        Gi = F.(si{1}){i}.Grids;
                
        % upper & bottom faces
        UFi = [repmat(upper_fi,size(Gi,2)-1,1) Gi(end,1:end-1)' Gi(end,2:end)'];
        BFi = [repmat(bottom_fi,size(Gi,2)-1,1) Gi(1,1:end-1)' Gi(1,2:end)'];
        
        % side faces
        LFi = [repmat(upper_fi,size(Gi,1)-1,1) Gi(1:end-1,1) Gi(2:end,1); upper_fi bottom_fi Gi(1,1)];
        RFi = [repmat(upper_fi,size(Gi,1)-1,1) Gi(1:end-1,end) Gi(2:end,end); upper_fi bottom_fi Gi(1,end)];
        
        % put into RV
        FV.(si{1})(i).Faces = [Fi; UFi; BFi; LFi; RFi];
        
        % compute volume
        V.(si{1})(i) = VolumeByTriangulation(FV.(si{1})(i).Faces,FV.(si{1})(i).Vertices);
        
    end
end

% restore transformation
L.T = Ltmp;
if( ~isequal(L.T,Ltmp) && nargout>1 )
    % we need to transform the vertices
    for si={'endo','epi'}
        for i=1:17
            if( i > numel(RV.(si{1})) ), continue; end
            FV.(si{1})(i) = [FV.(si{1})(i) ones(size(FV.(si{1})(i),1),1)] * L.T;
        end
    end
end