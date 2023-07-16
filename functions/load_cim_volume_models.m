function results = load_cim_volume_models (image_metadata_Dir,cim_modelDir)

main_directory = pwd;

caseList = dir([cim_modelDir '/*']);
caseList(ismember({caseList.name}, {'.', '..'})) = [];
nCases = length(caseList);

results = [];
for j=1:nCases
    caseID = caseList(j).name;
    inDir = [cim_modelDir '\' caseID]; 
    % identify sequence solutions in input folder
    modelList = dir([inDir '\model_*']);
    modelList = {modelList.name};
    nModels = numel(modelList);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cd(strcat(image_metadata_Dir,caseID))
    opts = detectImportOptions('image_metadata.csv');
    T = table2struct(readtable('image_metadata',opts));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:nModels
    %for i=1:1
        cd(main_directory)
        modelID = char(modelList(i));
        scanID = extractAfter(modelID,strcat(caseID,'_'));
        % find HR of model scan id
        idx = find(ismember(cell2mat({T.scan_id}'),str2double(scanID)));
        HR = T(idx).heart_rate;
        number_of_cycles= T(idx).number_of_cycles;
        result = extract_cim_model_summary(caseID,modelID,inDir,HR,number_of_cycles);
        results = [results,result];
    end
  clearvars T opts HR  
 cd (main_directory) 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                         Internal Functions                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Summary = extract_cim_model_summary(caseID,modelID,inDir,HR,number_of_cycles)
persistent L;
Summary = struct('ID',[],'model_name',[],...
    'ed_frame',[],'es_frame',[],'nframes',[],...
    'EDV',[], 'ESV',[], 'SV',[], 'MASS',[], 'EF',[], 'GLS',[], 'HR',[],'Volumes',[]);

GLS = struct('ID',[], 'nframes',[], 'edframe',[], 'esframe',[], ...
    'A2C',[], 'A3C',[], 'A4C',[], 'AvgGLS',[], 'GLS',[], ...
    'ant',[], 'antsep',[], 'infsep',[], 'inf',[], 'inflat',[], 'antlat',[]);

circvars = {'ant', 'antsep', 'infsep', 'inf', 'inflat', 'antlat'};



% reconstruct CIM model
fprintf('Reading CIM %s... ', modelID)
evalc('L = LVMM.ReadFromCIMFolder(inDir, modelID)');
fprintf('success\n')

% store volume data
GLS.ID = caseID;
Summary.ID   = caseID;
Summary.EDV  = L.EDV;
Summary.ESV  = L.ESV;
Summary.SV   = L.SV;
Summary.MASS = L.MASS;
Summary.EF   = L.EF;
Summary.HR = HR;


% calculate reference state length (ED)
ref_frame = L.es; % TODO: set ref_frame = ES frame
lvsm = L.model(ref_frame);
[L0_a2c, L0_a3c, L0_a4c, L0_circs] = GetEndoLengths(lvsm);

GLS.nframes = L.nframes;
GLS.edframe = L.ed;
GLS.esframe = L.es;

V = L.ComputeVolumes;
Volumes = V.endo;

Summary.Volumes   = Volumes;
Summary.model_name = L.name;
Summary.ed_frame = L.ed;
Summary.es_frame = L.es;
Summary.nframes = L.nframes;

% preallocate arrays
GLS.A2C = zeros(L.nframes, 1);
GLS.A3C = zeros(L.nframes, 1);
GLS.A4C = zeros(L.nframes, 1);

GLS.ant = zeros(L.nframes, 1); % circ 1
GLS.antsep = zeros(L.nframes, 1); % circ 2
GLS.infsep = zeros(L.nframes, 1); % circ 3
GLS.inf = zeros(L.nframes, 1); % circ 4
GLS.inflat = zeros(L.nframes, 1); % circ 5
GLS.antlat = zeros(L.nframes, 1); % circ 6

% calculate strain curve from ref frame
for f = ref_frame:L.nframes
    lvsm = L.model(f);
    [a2c, a3c, a4c, circs] = GetEndoLengths(lvsm);
    
    % strain (%) as change in length / original length for views
    GLS.A2C(f) = (a2c-L0_a2c)/L0_a2c*100;
    GLS.A3C(f) = (a3c-L0_a3c)/L0_a3c*100;
    GLS.A4C(f) = (a4c-L0_a4c)/L0_a4c*100;
    
    % strain (%) as change in length / original length for
    % individual directions
    for c = 1:length(circvars)
        GLS.(circvars{c})(f) = (circs(c)-L0_circs(c))/L0_circs(c)*100;
    end
    
end

% calculate average GLS
GLS.AvgGLS = mean([GLS.A2C, GLS.A3C, GLS.A4C], 2);

% passive GLS (between ES and ED2)
GLS.GLS  = mean([GLS.A2C(end), GLS.A3C(end), GLS.A4C(end)], 2);
Summary.GLS = mean([GLS.A2C(end), GLS.A3C(end), GLS.A4C(end)], 2);

clear L; 
end



function [a2c, a3c, a4c, wsegs] = GetEndoLengths(lvsm)
% Returns endocardial longitudinal lengths (a2c, a3c, and a4c views) from 
% LVMM object of [9 9] samples.
% TODO: generalise for all sampling

% verify default CIM surface point sampling
lvsm.nSamples = [9 9];

% 2-chamber view:
[w4, ~] = lvsm.GetLongitudinalPoints(9); % inferior wall
[w1, ~] = lvsm.GetLongitudinalPoints(25); % anterior wall
line4 = w4(:,:,1); line1 = w1(:,:,1); % get endo points only
a2cpoints = [line4; flip(line1(1:end-1, :))]; % combine, discard duplicate apex
a2c = SumPointDists(a2cpoints);

% 3-chamber view:
[w5, ~] = lvsm.GetLongitudinalPoints(13); % inferolateral wall
[w2, ~] = lvsm.GetLongitudinalPoints(29); % anteroseptal wall
line5 = w5(:,:,1); line2 = w2(:,:,1); % get endo points only
a3cpoints = [line5; flip(line2(1:end-1, :))]; % combine, discard duplicate apex
a3c = SumPointDists(a3cpoints);

% 4-chamber view:
[w3, ~] = lvsm.GetLongitudinalPoints(5); % inferoseptal wall
[w6, ~] = lvsm.GetLongitudinalPoints(21); % anterolateral wall
line3 = w3(:,:,1); line6 = w6(:,:,1); % get endo points only
a4cpoints = [line3; flip(line6(1:end-1, :))]; % combine, discard duplicate apex
a4c = SumPointDists(a4cpoints);

% store wall segment lengths:
wsegs = zeros(6, 1);
wsegs(1) = SumPointDists(line1);
wsegs(2) = SumPointDists(line2);
wsegs(3) = SumPointDists(line3);
wsegs(4) = SumPointDists(line4);
wsegs(5) = SumPointDists(line5);
wsegs(6) = SumPointDists(line6);

end


function len = SumPointDists(P)

diffs = diff(P, 1);
dists = sqrt(sum(diffs .* diffs, 2));
len = sum(dists);

end

end

