function Lmm = ReadFromCIMFolder( folder,model_name, varargin )
% LVMM/ReadFromCIMFolder is a static object factory function to create an LVMM object from a CIM 
% (Cardiac Image Modeller (tm)) software package.
%
%   L = LVMM.ReadFromCIMFolder(cim_folder);
%   L = LVMM.ReadFromCIMFolder(cim_folder,...);
%
% Available options:
%   - 'read_data', false | true
%     Read all CIM data: planes, volinfo and slices. Default is true.
%
% Author: Avan Suinesiaputra - Centre for Advanced Imaging, Univ. of Auckland (2012)
%Edited KG - to allow the file path to just be passed - removed searching
%for the models
Lmm = LVMM;

% check existence
if( ~exist(folder,'dir') )
    error('%s is not a valid folder.',folder);
end

% default options
opt.read_data = true;

% get options
for i=1:2:numel(varargin)
    if( isfield(opt,varargin{i}) ), opt.(varargin{i}) = varargin{i+1};
    else error('Unknown option %s.',varargin{i}); end
end


% set model's name
Lmm.name = model_name;

% read model files
model_files = dir(sprintf('%s/%s/*.model',folder,model_name));
model_files = cellfun(@(x) sprintf('%s/%s/%s',folder,model_name,x),{model_files.name},...
    'UniformOutput',false);

% sort frames
%[~,sidx] = sort(cellfun(@(x) str2num(x{1}),regexp(model_files,'_(\d+).model','tokens','once')));
B= regexp(model_files,'_(\d+).model','tokens');
B_ = vertcat(B{:});
B__ = B_(:,2);
[~,sidx] = sort(cellfun(@(x) str2num(x{1}),B__));
model_files = model_files(sidx);

% read one by one
for mi=1:numel(model_files)
    
    P = ReadCIMModelFileRaw(model_files{mi});
    
    Lmm.focalLengths(mi) = P.FocalLength;
    Lmm.lambdas(:,mi) = P.LambdaParameters;
    Lmm.mus(:,mi) = P.MuParameters;
    Lmm.thetas(:,mi) = P.ThetaParameters;
    
    if( isfield(P,'ModelToMagnetTransform') )
        Lmm.phase_data(mi).ModelToMagnetTransform = P.ModelToMagnetTransform;
    end
    
    if( isfield(P,'BasePlanePoints' ) )
        Lmm.phase_data(mi).BasePlanePoints = P.BasePlanePoints;

    end
    
    if( isfield(P,'BasePlaneNormal' ) )
        Lmm.phase_data(mi).BasePlaneNormal = P.BasePlaneNormal;
    end
    
end

% read system for ES & ED frames
system_file = dir(sprintf('%s%ssystem%s*.%s',folder,filesep,filesep,model_name));
if( numel(system_file)>0 )
    system_file = sprintf('%s%ssystem%s',folder,filesep,filesep,system_file(1).name);
    fprintf(1,'Reading %s\n',system_file);
    fid = fopen(system_file);
    if( fid < 3 ), error('Cannot read %s.',system_file); end
    
    tline = fgetl(fid);
    while( ischar(tline) )
        
        re = regexp(tline,'(End-systolic|End-diastolic)\s+Frame\s*:\s*(\d+)','tokens','once');
        if( ~isempty(re) && strcmpi(re{1},'end-diastolic') )
            Lmm.ed = str2num(re{2}) + 1;
        elseif( ~isempty(re) && strcmpi(re{1},'end-systolic') ) 
            Lmm.es = str2num(re{2}) + 1;
        end
        
        tline = fgetl(fid);
    end
    
    fclose(fid);
else
    Lmm.EstimateEDES;
end

% compute global LV functions
Lmm.UpdateGlobalLVFunctions;

% read cimm data
if( opt.read_data )
    data = ReadCIMData(folder,model_name);
    if( ~isempty(data) ), Lmm.data = data; end
end
