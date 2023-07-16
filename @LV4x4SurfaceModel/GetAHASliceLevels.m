function [mu_levels,idx,Pb] = GetAHASliceLevels(L)
% Get the AHA slice levels, i.e. apical tip (epi only), apex, mid and base
% boundaries.
%
%   [mu_levels,idx,Pb] = L.GetAHASliceLevels;
%
% Output:
%   - mu_levels is a struct with two fields: epi and endo. Each field
%     contains the slice level in the mu direction. 
%     For epi: [apical_tip, apex, mid, base]
%     For endo: [apex, mid base]
%
%   - idx is also a struct with two fields: epi and endo. Each field
%     contains cells of indices of points that belong to a particular level.
%
%   - Pb is the point samples of the AHA slice level borders.
%     Pb is a structure with two fields:
%     Pb.endo = Nx3x2 is slice level border points on endocardium
%     Pb.epi = Nx3x3 is slice level border poitns on epicardium
%     The number of points N is the same as L.nCirc (number of circumf. points)
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

% Slice level for 17-AHA models

% do not use any transformation
Ttmp = L.T;
L.T = [1 0 0; 0 1 0; 0 0 1; 0 0 0];

% the output: from apex to base
mu_levels = struct('epi',[],'endo',[]);
idx = struct('epi',{},'endo',{});
Pb = struct('epi',[],'endo',[]);

if( isempty(L) ), return; end
 
opt.tol = 1e-4;

% get the endo apical tip points
aptip = L.GetApicalPoints;
aptip = aptip(1,:,1);

% intersect endocardial apical tip with epicardium,
% the plane to intersect is perpendicular to the x-axes
P1 = L.GetIntersectionWithPlane(aptip,[1 0 0]);
P1 = P1{2};

% calculate its Xi, get the max of Xi(2)
Xi1 = L.FindSurfaceElementXi('epi',P1);
mu_levels.epi(1) = max(Xi1(:,2));

% divide the rest by 3
level_height = (4 - mu_levels.epi(1)) / 3;
mu_levels.epi(2:4) = mu_levels.epi(1) + (1:3)*level_height;

% get the indices for epi mu_level
epi_idx = {};
for i=numel(mu_levels.epi):-1:1
    xi = mu_levels.epi(i) - floor(mu_levels.epi(i)-0.0001);
    
    if( mu_levels.epi(i)-3>opt.tol ) % base takes all
        idx = union(find(L.E(:,3)==1 & L.E(:,4) >= 1 & L.E(:,4) <= 4 & L.E(:,2) <= xi ), ...
            find(L.E(:,3)==1 & L.E(:,4)>4));
    elseif( mu_levels.epi(i)-2>opt.tol ) % mid takes 5..16
        idx = union(find(L.E(:,3)==1 & L.E(:,4) >= 5 & L.E(:,4) <= 8 & L.E(:,2) <= xi ),...
            find(L.E(:,3)==1 & L.E(:,4)>8));
    elseif( mu_levels.epi(i)-1>opt.tol ) % apex takes 9..16
        idx = union(find(L.E(:,3)==1 & L.E(:,4) >= 9 & L.E(:,4) <=12 & L.E(:,2) <= xi ),...
            find(L.E(:,3)==1 & L.E(:,4)>12));
    else % the lowest one
        idx = find(L.E(:,3)==1 & L.E(:,4) >= 13 & L.E(:,2) <= xi);
    end
    
    if( isempty(epi_idx) ), epi_idx{1} = idx;
    else
        epi_idx{end} = setdiff(epi_idx{end},idx);
        epi_idx = [epi_idx idx];
    end
    
end
% reverse
epi_idx = fliplr(epi_idx);

% for the endo, the same xi level of endo must be translated into epi
% take only the top circular points of mid & apex levels
for i=2:3
    [~,imax] = max(L.surfacePSPoints(epi_idx{i},2));
    top_idx = epi_idx{i}(imax);
    
    if( L.E(top_idx,4) >= 13 )
        idx = intersect(find(L.E(:,2) == L.E(top_idx,2) & L.E(:,3)==1 & (L.E(:,4)==13 | L.E(:,4)==14)), ...
            find(L.E(:,1) == 0 | L.E(:,1)==1));
    elseif( L.E(top_idx,4) >= 9 )
        idx = intersect(find(L.E(:,2) == L.E(top_idx,2) & L.E(:,3)==1 & (L.E(:,4)==9 | L.E(:,4)==10)), ...
            find(L.E(:,1) == 0 | L.E(:,1)==1));
    elseif( L.E(top_idx,4) >= 5 )
        idx = intersect(find(L.E(:,2) == L.E(top_idx,2) & L.E(:,3)==1 & (L.E(:,4)==5 | L.E(:,4)==6)), ...
            find(L.E(:,1) == 0 | L.E(:,1)==1));
    elseif( L.E(top_idx,4) >= 1 )
        idx = intersect(find(L.E(:,2) == L.E(top_idx,2) & L.E(:,3)==1 & (L.E(:,4)==1 | L.E(:,4)==2)), ...
            find(L.E(:,1) == 0 | L.E(:,1)==1));
    else
        error('INTERNAL BUG !!');
    end

    % intersect plane with endocardium
    T = TriRep([1 2 3],L.surfacePoints(idx,1),L.surfacePoints(idx,2),L.surfacePoints(idx,3));
    Pi = L.GetIntersectionWithPlane(T.circumcenters,T.faceNormals);
    Pi = Pi{1};
    
    % calculate its Xi, get the min of Xi(2)
    Xi = L.FindSurfaceElementXi('endo',Pi);
    [mu_levels.endo(i-1),imax] = min(Xi(:,2));
    
    % add proper xi
    if( Xi(imax,4) <= 4 ), mu_levels.endo(i-1) = mu_levels.endo(i-1) + 3;
    elseif( Xi(imax,4) <= 8 ), mu_levels.endo(i-1) = mu_levels.endo(i-1) + 2;
    elseif( Xi(imax,4) <= 12 ), mu_levels.endo(i-1) = mu_levels.endo(i-1) + 1;
    end
   
end

% the last one is the top base
mu_levels.endo(3) = 4;

% get the indices for endo mu_level
endo_idx = {};
for i=numel(mu_levels.endo):-1:1
    xi = mu_levels.endo(i) - floor(mu_levels.endo(i)-0.0001);
    
    if( mu_levels.endo(i)-3>opt.tol ) % base takes all
        idx = union(find(L.E(:,3)==0 & L.E(:,4) >= 1 & L.E(:,4) <= 4 & L.E(:,2) <= xi ), ...
            find(L.E(:,3)==0 & L.E(:,4)>4));
    elseif( mu_levels.endo(i)-2>opt.tol ) % mid takes 5..16
        idx = union(find(L.E(:,3)==0 & L.E(:,4) >= 5 & L.E(:,4) <= 8 & L.E(:,2) <= xi ),...
            find(L.E(:,3)==0 & L.E(:,4)>8));
    elseif( mu_levels.endo(i)-1>opt.tol ) % apex takes 9..16
        idx = union(find(L.E(:,3)==0 & L.E(:,4) >= 9 & L.E(:,4) <=12 & L.E(:,2) <= xi ),...
            find(L.E(:,3)==0 & L.E(:,4)>12));
    else % the lowest one
        idx = find(L.E(:,3)==0 & L.E(:,4) >= 13 & L.E(:,2) <= xi);
    end
    
    if( isempty(endo_idx) ), endo_idx{1} = idx;
    else
        endo_idx{end} = setdiff(endo_idx{end},idx);
        endo_idx = [endo_idx idx];
    end
    
end
% reverse
endo_idx = fliplr(endo_idx);


% prepare output
clear idx;
idx.epi = epi_idx;
idx.endo = endo_idx;


% restore transformation
L.T = Ttmp;

% compute boundaries if required
if( nargout >= 3 )

    xi_thetas = L.E(L.map(1,:),1);

    % EPI
    eids = [];
    for i=1:3
        if( mu_levels.epi(i)>=3 )
            eids(:,i) = L.E(L.map(end,:),4) - 12;
        elseif( mu_levels.epi(i)>=2 )
            eids(:,i) = L.E(L.map(end,:),4) - 8;
        elseif( mu_levels.epi(i)>=1 )
            eids(:,i) = L.E(L.map(end,:),4) - 4;
        else
            eids(:,i) = L.E(L.map(end,:),4);
        end    
    end

    xi_mu_epi = mu_levels.epi(1:3) - floor(mu_levels.epi(1:3) - 0.00001);
    XI.epi = [repmat(xi_thetas(:),3,1) reshape(ones(numel(xi_thetas),1)*xi_mu_epi,[],1) ones(L.nCirc*3,1) eids(:)];

    % ENDO
    eids = [];
    for i=1:2
        if( mu_levels.endo(i)>=3 )
            eids(:,i) = L.E(L.map(end,:),4) - 12;
        elseif( mu_levels.endo(i)>=2 )
            eids(:,i) = L.E(L.map(end,:),4) - 8;
        elseif( mu_levels.endo(i)>=1 )
            eids(:,i) = L.E(L.map(end,:),4) - 4;
        else
            eids(:,i) = L.E(L.map(end,:),4);
        end    
    end

    xi_mu_endo = mu_levels.endo(1:2) - floor(mu_levels.endo(1:2) - 0.00001);
    XI.endo = [repmat(xi_thetas(:),2,1) reshape(ones(numel(xi_thetas),1)*xi_mu_endo,[],1) zeros(L.nCirc*2,1) eids(:)];

    % Surface points of the boundaries
    Xb = L.GetSurfacePointsFromElements([XI.endo; XI.epi]);

    % reorder 
    Pb.endo = permute(reshape(Xb(1:(2*L.nCirc),:),L.nCirc,[],3),[1 3 2]);
    Pb.epi = permute(reshape(Xb((2*L.nCirc)+1:end,:),L.nCirc,[],3),[1 3 2]);
end
