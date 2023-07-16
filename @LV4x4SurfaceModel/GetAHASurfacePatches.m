function [FV,Pbx] = GetAHASurfacePatches(L,G)
% Get regional AHA patches
%
%   [FV,Pb] = L.AHASurfacePatches(G)
%
% Input: G must be the output from L.AHAPoints
% Output: 
%   - FV is a struct of 'endo' and 'epi' fields with 17 cells. Each cell
%     contains Faces & Vertices struct that can be fed into patch command.
%   - Pb contains border points
%
% Author: Avan Suinesiaputra - Centre for Advanced Imaging, Univ. of Auckland (2012)

% handle empty case
if( isempty(L) )
    FV = [];
    return;
end
    
% reset transformation to identity
Ttmp = L.T;
L.T = [1 0 0; 0 1 0; 0 0 1; 0 0 0];

% ------------------------------------------------------------------------------------------
% COMPUTE NEW ELEMENTS OF THE SLICE LEVEL BORDERS

% based on grid G indices, I must compute the slice level borders

% GM is mu values but they are continuous from 0 .. 1 .. 2 .. 3 .. 4
GM = L.E(:,2) + (3-floor((L.E(:,4)-1)/4));

% get fixed Xi thetas
xi_thetas = L.E(L.map(1,:),1);

% there are 5 borders: [endo:bm endo:ma epi:bm epi:ma epi:aa]
% where: bm = base-with-mid, ma = mid-with-apex, aa = apex-with-apicaltip

% so xi for lambda borders are:
xi_lambda_borders = [0 0 1 1 1]; % 2 for endo, 3 for epi

% xi mu border levels
% base-mid at endo
xi_mu_borders(1) = GM(G.endo{7}(end,1)) + 0.5*(GM(G.endo{1}(1,1)) - GM(G.endo{7}(end,1)));
% mid-apex at endo
xi_mu_borders(2) = GM(G.endo{13}(end,1)) + 0.5*(GM(G.endo{7}(1,1)) - GM(G.endo{13}(end,1)));
% base-mid at epi
xi_mu_borders(3) = GM(G.epi{7}(end,1)) + 0.5*(GM(G.epi{1}(1,1)) - GM(G.epi{7}(end,1)));
% mid-apex at epi
xi_mu_borders(4) = GM(G.epi{13}(end,1)) + 0.5*(GM(G.epi{7}(1,1)) - GM(G.epi{13}(end,1)));
% apex-tip at epi
xi_mu_borders(5) = GM(G.epi{17}(1,1)) + 0.5*(GM(G.epi{13}(1,1)) - GM(G.epi{17}(1,1)));

% the idea is to create element vectors, just like L.E but for the slice level borders.
% these are additional elements, let's call it BE
BE = [];
for i=1:numel(xi_mu_borders)
    % create element vector
    if( xi_mu_borders(i) > 3 )
        E = L.E(L.map(1,:),4);
        M = repmat(xi_mu_borders(i)-3,size(xi_thetas));
    elseif( xi_mu_borders(i) > 2 )
        E = L.E(L.map(10,:),4);
        M = repmat(xi_mu_borders(i)-2,size(xi_thetas));
    elseif( xi_mu_borders(i) > 1 )
        E = L.E(L.map(18,:),4);
        M = repmat(xi_mu_borders(i)-1,size(xi_thetas));
    else
        E = L.E(L.map(26,:),4);
        M = repmat(xi_mu_borders(i),size(xi_thetas));
    end
    
    BE = [BE; xi_thetas M repmat(xi_lambda_borders(i),size(xi_thetas)) E];
end

% now BE is 5 x L.nCirc matrix
% each row of BE is the element of slice level borders

% ------------------------------------------------------------------------------------------
% CLOSE GAPS BETWEEN SEGMENTS ON EACH SLICE LEVEL

% next is to expand grid G that closes all segments, so there is no gap between segments on each level
% we put the new grid in Gf

Gf = G;

segs = [3 4 5 6 1 9 10 11 12 7 13 14 15 16];
segs_adj = [4 5 6 1 2 10 11 12 7 8 14 15 16 13];

for si=1:numel(segs)
    
    % add the leftmost grid of the adjacent to the segments
    Gf.endo{segs(si)}(:,end+1) = Gf.endo{segs_adj(si)}(:,1);
    Gf.epi{segs(si)}(:,end+1) = Gf.epi{segs_adj(si)}(:,1);
    
end

% ------------------------------------------------------------------------------------------
% CREATING FACES & VERTICES

% from this on, we're going to work on RC space
% we need to return the space transformation first
L.T = Ttmp;

% collect the border slice points, by calculating the surface points of slice level borders BE
Pb = L.GetSurfacePointsFromElements(BE); 
Pb = permute(reshape(Pb,L.nCirc,[],3),[1 3 2]);  % Pb is now [L.nCirc x 3 x 5] matrix

% to make it easier
Pbx.endo = Pb(:,:,1:2);
Pbx.epi = Pb(:,:,3:end);

% we'll browse through surfaces and then the 17 segments
offset = [0 numel(L.map)];
segment_names = {'endo','epi'};
for si=1:numel(segment_names)
    s = segment_names{si};
    
    for cs=1:17

        if( isempty(Gf.(s){cs}) ), continue; end;

        % get faces from Gf
        [idx,~,Fi] = unique(Gf.(s){cs});
        Fi = reshape(Fi,size(Gf.(s){cs}));
        
        % get vertices from surface points
        Vi = L.surfacePoints(idx,:);
        
        % need to add slice level borders
        [mpos,~] = find((repmat(offset(si)+L.map(:),1,numel(Gf.(s){cs}(1,:))) - repmat(Gf.(s){cs}(1,:),numel(L.map),1)) == 0);
        if( cs<=6 ) 
            
            % add bottom border to the top
            B = Pbx.(s)(floor((mpos-1)/size(L.map,1))+1,:,1);
            
            Fi = [1:size(B,1); Fi+size(B,1) ];
            Vi = [B; Vi];
            
        elseif( cs<=12 )
            
            % add bottom border to the top
            B = Pbx.(s)(floor((mpos-1)/size(L.map,1))+1,:,2);
            
            Fi = [1:size(B,1); Fi+size(B,1) ];
            Vi = [B; Vi];

            % add top border to the bottom
            B = Pbx.(s)(floor((mpos-1)/size(L.map,1))+1,:,1);
            
            Fi = [Fi; size(Vi,1)+(1:size(B,1))];
            Vi = [Vi; B];
            
        elseif( cs<=16 )

            % only for epi: add bottom borders to the top
            if( si==2 )
                B = Pbx.(s)(floor((mpos-1)/size(L.map,1))+1,:,3);
                
                Fi = [1:size(B,1); Fi+size(B,1) ];
                Vi = [B; Vi];
            end
            
            % add top borders to the bottom
            B = Pbx.(s)(floor((mpos-1)/size(L.map,1))+1,:,2);
        
            Fi = [Fi; size(Vi,1)+(1:size(B,1))];
            Vi = [Vi; B];
            
        elseif( cs==17 )
            
            B = Pbx.(s)(floor((mpos-1)/size(L.map,1))+1,:,3);

            Fi = [1:size(B,1); Fi+size(B,1) ];
            Vi = [B; Vi];
                
        end
        
        % put Fi as grids
        FV.(s){cs}.Grids = Fi;

        % create half triangle faces
        F1 = reshape(Fi(1:end-1,1:end-1),[],1);
        F1(:,2) = reshape(Fi(2:end,1:end-1),[],1);
        F1(:,3) = reshape(Fi(1:end-1,2:end),[],1);

        % complete another half triangle faces
        F2 = F1(:,2:3);
        F2(:,3) = reshape(Fi(2:end,2:end),[],1);
        
        % put into FV
        FV.(s){cs}.Faces = [F1; F2];
        FV.(s){cs}.Vertices = Vi;
        
    end
end
