function WriteAsCIMModel(L,cim_file)
% Write a model into a CIM file
%
%   L.WriteAsCIMModel(cim_file);
%
% Note: CIM has .model suffix, but this is not required for this function.
%
% Author: Avan Suinesiaputra - Centre for Advanced Imaging, Univ. of Auckland (2012)

if( isempty(L) ), return; end

% check cim_file
if( exist(cim_file,'file') )
    fprintf(1,'File %s already exist.\n',cim_file);
    yn = input('Do you want to overwrite it (Y/[N]) ? ','s');
    if( ~strcmpi(yn,'y') ), return; end
end

% open file
fid = fopen(cim_file,'w+');
if( fid < 3 ), error('Cannot open %s for write.',cim_file); end

% save FocalLength
fprintf(fid,'FocalLength:\n%f\n',L.focalLength);

% save ModelToMagnetTransform
if( isfield(L.data,'ModelToMagnetTransform') )
    fprintf(fid,'\nModelToMagnetTransform:\n');
    fprintf(fid,'a1\t%fi\t%fj\t%fk\n',...
        L.data.ModelToMagnetTransform(1,2),L.data.ModelToMagnetTransform(1,2),L.data.ModelToMagnetTransform(1,3));
    fprintf(fid,'a2\t%fi\t%fj\t%fk\n',...
        L.data.ModelToMagnetTransform(2,2),L.data.ModelToMagnetTransform(2,2),L.data.ModelToMagnetTransform(2,3));
    fprintf(fid,'a3\t%fi\t%fj\t%fk\n',...
        L.data.ModelToMagnetTransform(3,2),L.data.ModelToMagnetTransform(3,2),L.data.ModelToMagnetTransform(3,3));
    fprintf(fid,'t\t%fi\t%fj\t%fk\n',...
        L.data.ModelToMagnetTransform(4,2),L.data.ModelToMagnetTransform(4,2),L.data.ModelToMagnetTransform(4,3));
else
    fprintf(fid,'\nModelToMagnetTransform:\n');
    fprintf(fid,'a1\t%fi\t%fj\t%fk\n',L.T(1,2),L.T(1,2),L.T(1,3));
    fprintf(fid,'a2\t%fi\t%fj\t%fk\n',L.T(2,2),L.T(2,2),L.T(2,3));
    fprintf(fid,'a3\t%fi\t%fj\t%fk\n',L.T(3,2),L.T(3,2),L.T(3,3));
    fprintf(fid,'t\t%fi\t%fj\t%fk\n',L.T(4,2),L.T(4,2),L.T(4,3));
end

% save LambdaParameters
fprintf(fid,'\nLambdaParameters:\n%d\n',numel(L.lambdaParams));
for i=1:numel(L.lambdaParams)
    fprintf(fid,'%f\n',L.lambdaParams(i));
end

% save MuParameters
fprintf(fid,'\nMuParameters:\n%d\n',numel(L.muParams));
for i=1:numel(L.muParams)
    fprintf(fid,'%f\n',L.muParams(i));
end

% save ThetaParameters
fprintf(fid,'\nThetaParameters:\n%d\n',numel(L.thetaParams));
for i=1:numel(L.thetaParams)
    fprintf(fid,'%f\n',L.thetaParams(i));
end

% save BasePlanePoints
if( isfield(L.data,'BasePlanePoints') )
    fprintf(fid,'\nBasePlanePoints:\n%fi\t%fj\t%fk\n',...
        L.data.BasePlanePoints(1),L.data.BasePlanePoints(2),L.data.BasePlanePoints(3));
end

% save BasePlaneNormal
if( isfield(L.data,'BasePlaneNormal') )
    fprintf(fid,'\nBasePlaneNormal:\n%fi\t%fj\t%fk\n',...
        L.data.BasePlaneNormal(1),L.data.BasePlaneNormal(2),L.data.BasePlaneNormal(3));
end

% close
fclose(fid);

fprintf(1,'Model was saved to %s successfully.\n',cim_file);