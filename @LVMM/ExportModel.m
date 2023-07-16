function ExportModel(L,model_folder,varargin)
% Save LVMM models in a folder, similar to CIM model folder
%
%   L.ExportModel(model_folder,...);
%
% Available options:
%   - 'name', string.
%     Alter the model name to save as filenames. Default is to use L.name
%
% Author: Avan Suinesiaputra - University of Auckland 2012

if( L.nframes < 1 ), return; end

% get options
opt.name = L.name;

for i=1:2:numel(varargin)
    if( isfield(opt,varargin{i}) ), opt.(varargin{i}) = varargin{i+1};
    else error('Unknown optional argument.'); end
end

lvdir = [model_folder '/' L.name];
if( ~exist(lvdir,'dir') ), mkdir(lvdir); end

% save models
for i=1:L.nframes
    Lv = L.model(i);
    Lv.WriteAsCIMModel([lvdir '/' opt.name '_' num2str(i) '.model']);
end

% save spec file
specdir = [model_folder '/system'];
if( ~exist(specdir,'dir') ), mkdir(specdir); end

fid = fopen([specdir '/' opt.name '.' L.name],'w+');
if( fid < 3 )
    fprintf(2,'Cannot create spec file.\n');
else
    
    fprintf(fid,'#EXPORTED FROM LVMM CLASS\n\n');
    
    fprintf(fid,'Num Frames : %d\n',L.nframes);
    fprintf(fid,'End-systolic Frame : %d\n',L.es-1);    % START FROM 0 !!!
    fprintf(fid,'End-diastolic Frame : %d\n',L.ed-1);
    
end
fclose(fid);

% save time information if exist
timeinfo = {'TimeIntervals','NormalisedTimeIntervals'};
if( any(isfield(L.data,timeinfo)) )
    
    fid = fopen([lvdir '/TimeIntervals.data'],'w+');
    if( fid < 3 )
        fprintf(2,'Cannot create TimeIntervals.data');
    else
        
        for ti = 1:numel(timeinfo)
            if( isfield(L.data,timeinfo{ti}) )
                
                x = L.data.(timeinfo{ti});
                fprintf(fid,'\n%s:\n%d\n',timeinfo{ti},numel(x));
                for xi=1:numel(x)
                    fprintf(fid,'%.10f\n',x(xi));
                end
                
            end
        end
        
    end
    fclose(fid);
    
end