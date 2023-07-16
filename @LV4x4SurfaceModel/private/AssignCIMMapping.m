function L = AssignCIMMapping(L)
% This function assign CIM's GLM mapping to an LV object
%
% Author: Avan Suinesiaputra - Centre of Advanced Imaging, Univ. of Auckland (2012)

% read the global to element CIMM conversion (it's fixed)
% mind this weird thing for lambda
X = load('BgToBe_1.txt');

% apply this mapping for 32 nodes 
localNodes = mod(X(:,1),32);
localNodes(localNodes == 0)=32; %%%Accounts for mod(32)

%Transform the local nodes to our numbering scheme using the map below
map=[1,2,5,6, 4,3,8,7, 13,14,9,10, 16,15,12,11,...
    17,18,21,22, 20,19,24,23, 29,30,25,26, 32,31,28,27];
L.G2E.Lambda = sparse((ceil(X(:,1)/32)-1)*32+map(localNodes)',X(:,2),X(:,3));

X = load('BgToBe_2.txt'); 
L.G2E.Mu = sparse(X(:,1),X(:,2),X(:,3));

X = load('BgToBe_2.txt'); 
L.G2E.Theta = sparse(X(:,1),X(:,2),X(:,3));

