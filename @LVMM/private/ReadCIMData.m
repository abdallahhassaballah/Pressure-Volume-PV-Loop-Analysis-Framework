function L = ReadCIMData(cim_folder,mfolder)
% This function reads a CIM folder and update the LVMotionModel object L.
%
%   L = ReadCIMData(cim_folder,model_name);
%
% Author: Avan Suinesiaputra - Centre for Advanced Imaging, Univ. of Auckland (2012)


L = [];
if( ~exist(cim_folder,'dir') )
    error('Folder %s does not exist',cim_folder);
end

% get the plane information
plane_files = dir(sprintf('%s/planes/*.planes_%s',cim_folder,mfolder));
plane_files = plane_files(cellfun(@(x) x(1)~='.',{plane_files.name}));
for i=1:numel(plane_files)

    fprintf(1,'Reading plane file %s\n',plane_files(i).name);
    fid = fopen(sprintf('%s/planes/%s',cim_folder,plane_files(i).name),'r');

    % first line the number of slices
    nslices = str2num(fgetl(fid));
    planes = struct(...
        'ImagePositionOffset',zeros(nslices,3),...
        'ImageOrientationOffset',zeros(nslices,6));

    for si=1:nslices
        plane = zeros(3,6);
        for j=1:3
            tmp = str2num(fgetl(fid));
            if( numel(tmp)==3 ), tmp = [tmp tmp]; end
            plane(j,:) = tmp;
        end

        % process
        tlc = plane(1,:);

        v1 = plane(2,:) - tlc;
        v1(1:3) = v1(1:3) ./ norm(v1(1:3));
        v1(4:6) = v1(4:6) ./ norm(v1(4:6));

        v2 = plane(3,:) - tlc;
        v2(1:3) = v2(1:3) ./ norm(v2(1:3));
        v2(4:6) = v2(4:6) ./ norm(v2(4:6));

        % store shifting
        planes.ImagePositionOffset(si,:) = tlc(1:3) - tlc(4:6);
        planes.ImageOrientationOffset(si,:) = [v1(1:3) v2(1:3)] - [v1(4:6) v2(4:6)];
        
        % this is how to use it:
        % newPos = planes.ImagePositionOffset(i,:)' + dcm.ImagePositionPatient;
        % newOri = planes.imageOrientationOffset(i,:)' + dcm.ImageOrientationPatient;

    end

    fclose(fid);

    % store
    series_name = regexp(plane_files(i).name,'(.+).planes_','tokens','once');
    series_name = series_name{1};

    L.planes.(series_name) = planes;

end

% get image info
iifiles = dir(sprintf('%s/planes/*.imageInfo',cim_folder));
iifiles = iifiles(cellfun(@(x) x(1)~='.',{iifiles.name}));

for i=1:numel(iifiles)

    fprintf(1,'Reading image info file %s\n',iifiles(i).name);
    fid = fopen(sprintf('%s/planes/%s',cim_folder,iifiles(i).name),'r');

    slices = [];
    tline = fgetl(fid);
    while( ischar(tline) )

        sinfo = regexp(tline,'Slice\s+(\d+)\s+([^:]+):\s+(\S+)','tokens','once');
        if( isempty(sinfo) )
            tline = fgetl(fid);
            continue; 
        end

        si = str2num(sinfo{1});
        sfield = strrep(lower(deblank(sinfo{2})),' ','_');
        slices(si).(sfield) = str2double(sinfo{3});

        tline = fgetl(fid);
    end

    fclose(fid);

    % store
    iinfo_name = regexp(iifiles(i).name,'(.+).imageInfo','tokens','once');
    iinfo_name = iinfo_name{1};

    L.slices.(iinfo_name) = slices;

end

% get volume info
% volinfo = ExtractCIMVolumeInfo(cim_folder,mfolder);
% if( ~isempty(volinfo) )
%     L.volinfo = volinfo;
% end

% read guide points
% guide_points = ReadCIMGuidePoints(cim_folder,mfolder);
% if( ~isempty(guide_points) ), L.guide_points = guide_points; end


% read disabled series
gp_spec_file = strcat(cim_folder,filesep,mfolder,filesep,'gp_window_spec.gpt');
if( exist(gp_spec_file,'file') )
    fprintf(1,'Guide points specification file is found.\n');

    % open spect
    fid = fopen(gp_spec_file,'r');
    if( fid < 3 ), error('Cannot read file %s.',gp_spec_file); end
    
    tline = fgetl(fid);
    while( ischar(tline) )
        
        re = regexp(tline,'Slice disabled slice\s*(\d+)\s*:\s*(FALSE|TRUE)','tokens','once');
        if( numel(re)==2 )
            disabled_series(str2num(re{1})) = strcmpi(re{2},'true');
        end
        
        tline = fgetl(fid);
    end
    
    % close spec
    fclose(fid);
    
    % add to data
    i = 1;
    series = fieldnames(L.slices);
    for si=1:numel(series)
        L.disabled.(series{si}) = disabled_series(i:(i-1)+numel(L.slices.(series{si})));
        i = i + numel(L.slices.(series{si}));
    end
    
end


% read frames (trigger time) if any
frame_files = dir(sprintf('%s%sframes%s*.frames',cim_folder,filesep,filesep));
for i=1:numel(frame_files)
    
    frame_file = sprintf('%s%sframes%s%s',cim_folder,filesep,filesep,frame_files(i).name);
    fid = fopen(frame_file,'r');
    if( fid < 3 ), continue; end
    series_name = strrep(frame_files(i).name,'.frames','');
    fprintf(1,'Reading frames for %s\n',series_name);
    
    time = [];
    tline = fgetl(fid);
    while( ischar(tline) )
        
        re = regexp(tline,'Frame\s+Time\s+(\d+)\s+:\s+([0-9.-]+)','tokens','once');
        if( ~isempty(re) )
            time(str2num(re{1})) = str2double(re{2});
        end
        
        tline = fgetl(fid);
    end
    L.frames.(series_name) = time;
    
    fclose(fid);
    
end


% read time intervals information
time_file = [cim_folder filesep mfolder filesep 'TimeIntervals.data'];
if( exist(time_file,'file') )

    fid = fopen(time_file,'r');
    if( fid >= 3 )
        
        tline = fgetl(fid);
        while( ischar(tline) )
            
            re = regexp(tline,'(\w+)\s*:','tokens','once');
            if( ~isempty(re) )
                
                n = str2num(fgetl(fid));
                time_data = textscan(fid,'%f',n);
                
                L.(re{1}) = time_data{1}';
                
            end
            tline = fgetl(fid);
            
        end
        
    end
    fclose(fid);
    
end

