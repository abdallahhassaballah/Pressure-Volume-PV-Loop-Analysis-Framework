function V = ExtractCIMVolumeInfo(cimfolder,model_name)
% This private class function extracts volume information from a CIM folder
% for a given model_name.
%
% Author: Avan Suinesiaputra - Centre for Advanced Imaging, Univ. of Auckland (2012)

voldir = sprintf('%s/volumes/info/',cimfolder);
if( ~exist(voldir,'dir') ), V=[]; return; end

volfile = dir(sprintf('%s/*_ca.%s',voldir,model_name));
volfile = {volfile.name};

if( isempty(volfile) ), return; end
fid = fopen(sprintf('%s/%s',voldir,volfile{1}));
if( fid<3 )
    error('Error opening %s/%s',voldir,volfile{1});
end

V = struct;
cs = 0;
sections = {'APEX','MID','BASE','TOTAL'};
vols = {'Epi','Endo','Wall','WT'};
gf = {'EDV','ESV','SV','EF','MASS'};

css = 0;
subs = {'septal','inferior','lateral','anterior','postero septal','inferior lateral',...
    'anterior lateral','anterior septal'};

tline = fgetl(fid);
while( true )
    if( ~ischar(tline) ), break; end
    
    i = find(strcmpi(sections,strtrim(tline)));
    if( ~isempty(i) ), cs = i; tline = fgetl(fid); continue; end
    if( cs==0 ), tline = fgetl(fid); continue; end
    
    i = find(strcmpi(subs,strtrim(tline)));
    if( ~isempty(i) ), css = i; tline = fgetl(fid); continue; end
    
    % check GF
    x = regexp(tline,'\s*\$(\S+)\s*:\s*(\S+)','tokens','once');
    if( ~isempty(x) && any(strcmpi(gf,x{1})) )
        if( strcmpi(sections{cs},'total') )
            V.(sections{cs}).(x{1}) = str2double(x{2});
        else
            V.(sections{cs}).(strrep(subs{css},' ','_')).(x{1}) = str2double(x{2});
        end
        tline = fgetl(fid);
        continue;
    end

    % volumes
    x = regexp(tline,'\s*(\S+)\s*:\s*(.*)','tokens','once');
    if( ~isempty(x) && any(strcmpi(vols,x{1})) )
        if( strcmpi(sections{cs},'total') )
            V.(sections{cs}).(x{1}) = str2num(x{2});
        else
            V.(sections{cs}).(strrep(subs{css},' ','_')).(x{1}) = str2num(x{2});
        end
        tline = fgetl(fid);
        continue;
    end
    
    tline = fgetl(fid);
end

fclose(fid);

end