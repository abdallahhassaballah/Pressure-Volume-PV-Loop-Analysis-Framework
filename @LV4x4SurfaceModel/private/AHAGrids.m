function G = AHAGrids(L,R)
% Reorder and create grids for each AHA points
%
%   G = L.AHAGrids(R)
%
% Input R must be the output from L.AHAPoints;
% Output G is a struct with endo and epi fields, each has 17 cells.
%
% Author: Avan Suinesiaputra - Centre for Advanced Imaging, Univ. of Auckland (2012)

% do not use any transformation
Ttmp = L.T;
L.T = [1 0 0; 0 1 0; 0 0 1; 0 0 0];

% split endo & epi
Rx.epi = cell(1,17);
Rx.endo = cell(1,17);

n = size(L.surfacePoints,1)/2;
for i=1:17
    Rx.epi{i} = R{i}(R{i}>n);
    Rx.endo{i} = R{i}(R{i}<=n);
end

% reordering

% we have to make unique elements for all segments
LE = L.E;
LE(:,1) = (1-LE(:,1)) + mod(LE(:,4)-1,4);
LE(:,2) = LE(:,2) + floor((LE(:,4)-1)/4);

G.epi = cell(1,17);
G.endo = cell(1,17);
ss = {'endo','epi'};

% 1-16
for cs=1:16    
    for si=1:numel(ss)
        
        Ps = L.surfacePSPoints(Rx.(ss{si}){cs},:);
        E = LE(Rx.(ss{si}){cs},:);

        % create index grid
        thetas = unique(E(:,1));
        Vi = [];
        for i=1:numel(thetas)

            % get the points with the same theta level
            idx = find(E(:,1)==thetas(i));

            % sort based on mu
            [~,ii] = sort(Ps(idx,2));

            % reindex and put into Vi
            Vi(:,i) = idx(ii(:));

        end

        % use global index
        G.(ss{si}){cs} = Rx.(ss{si}){cs}(Vi);

    end
end

% special case for 14, because it is cut on theta = 0
theta_14 = L.surfacePSPoints(G.endo{14}(1,:),1);
idx = theta_14<pi/2;
theta_14(idx) = 2*pi + theta_14(idx)+0.00001;
[~,idx] = sort(theta_14);

G.endo{14} = G.endo{14}(:,idx);
G.epi{14} = G.epi{14}(:,idx);


% 17 only for epi

Ps = L.surfacePSPoints(Rx.epi{17},:);
E = L.E(Rx.epi{17},:);

% create index grid: column: theta, row: mu
mus = flipud(unique(E(:,2)));
Vi = [];
for i=1:numel(mus)
    
    % get points with the same mu level
    idx = find(E(:,2)==mus(i));
    
    % sort based on theta, but take the global theta value
    [~,ii] = sort(Ps(idx,1));
    
    % reindex idx and put it to Fi
    Vi(i,:) = idx(ii(:))';
    
end

G.epi{17} = Rx.epi{17}(Vi);

% restore transformation
L.T = Ttmp;