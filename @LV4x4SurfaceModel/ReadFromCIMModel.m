function L = ReadFromCIMModel(L,cim_file)
% Read LV model from CIM model
%
%   L.ReadFromCIMModel(cim_file);
%
% Input is a CIM filename, which usually ends with .model suffix
%
% Author: Avan Suinesiaputra - Auckland Bioengineering Institute (2011)

% read CIM file
% CIM model parameters
params = {
    'FocalLength',            'scalar', []; ...
    'ModelToMagnetTransform', 'invIJK',  4; ...
    'LambdaParameters',       'vector', [];...
    'MuParameters',           'vector', []; ...
    'ThetaParameters',        'vector', []; ...
    'BasePlanePoints',        'invIJK',  1; ...
    'BasePlaneNormal',        'invIJK',  1, ...
};
P = struct;

% note about the parameter types:
%   - scalar means there is only one value in the next line
%   - vector means that there are manu scalar values and the number of
%     elements is given in the next line
%   - invIJK is an inverted matrix affixed by i, j and k characters, and
%     the number of lines is given in the next element of the params

% read the file
fid = fopen(cim_file,'r');
if( fid < 3 ), error('Cannot read file %s',cim_file); end

tline = fgetl(fid);
while ischar(tline)
    % skip empty line
    if( isempty(strtrim(tline)) )
        tline = fgetl(fid);
        continue;
    end
    
    % find which parameter it is
    pi = find(strcmpi(params(:,1),tline(1:end-1)));
    if( isempty(pi) )
        error('Unhandled parameter type %s',tline);
    end
    
    %fprintf(1,'Reading %s\n',params{pi,1});
    if( strcmpi(params{pi,2},'scalar') )
        % read the next line only
        P.(params{pi,1}) = str2num(fgetl(fid));
    elseif( strcmpi(params{pi,2},'vector') )
        % read the number of elements
        nelmt = textscan(fid,'%d', 1);
        % read the data
        C = textscan(fid,'%f',nelmt);
        P.(params{pi,1}) = C{1};
    elseif( strcmpi(params{pi,2},'invIJK') )
        nrows = params{pi,3};
        M = [];
        for i=1:nrows
            tline = fgetl(fid);
            C = regexp(tline,'([0-9e.\-]+)i\s+([0-9e.\-]+)j\s+([0-9e.\-]+)k','tokens','once');
            M = [M cellfun(@(x) str2num(x),C)'];
        end
        P.(params{pi,1}) = M;
    else
        % unknown parameter type
        warning('Unknown parameter type of %s. Skip.\n',params{pi,1});
        continue; 
    end
    
    % next line
    tline = fgetl(fid); 
end

% close file
fclose(fid);


% setup the model (it's CIMM fixed node mapping)
L = AssignCIMMapping(L);

% save focal length & update model
L.focalLength = P.FocalLength;
L.lambdaParams = P.LambdaParameters;
L.muParams = P.MuParameters;
L.thetaParams = P.ThetaParameters;

% save other to data
if( isfield(P,'ModelToMagnetTransform') )
    L.data.ModelToMagnetTransform = P.ModelToMagnetTransform';
end
if( isfield(P,'BasePlanePoints') )
    L.data.BasePlanePoints = P.BasePlanePoints;
end
if( isfield(P,'BasePlaneNormal') )
    L.data.BasePlaneNormal = P.BasePlaneNormal;
end

L = UpdateModel(L);
