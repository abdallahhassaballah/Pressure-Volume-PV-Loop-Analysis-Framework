function [PreResults,Results_new] = Pressure_Analysis_Results (result_path,study_id)
fs = 240;
main_dir = result_path;
[file_paths_ECG,file_paths_LV] = Find_all_eligible_snapshots_for_ED_analysis(result_path);
nFiles = size(file_paths_ECG,1);
N = 1;
ALL_Results = [];
for ii = 1:nFiles
    % Extract location of P, Q, R, S and T waves in ECG
    [data, ~] = loadh5(char(file_paths_ECG(ii)));
    Data      = importdata(char(file_paths_LV(ii)));

    ECG = Data(:,1); Pressure = Data(:,2);
    
    P_Offsets = data.waves_peak.ECG_P_Offsets+1;
    P_Onsets  = data.waves_peak.ECG_P_Onsets+1;
    P_Peaks   = data.waves_peak.ECG_P_Peaks+1;

    Q_Peaks = data.waves_peak.ECG_Q_Peaks+1;

    R_Offsets = data.waves_peak.ECG_R_Offsets+1;
    R_Onsets  = data.waves_peak.ECG_R_Onsets+1;
    R_Peaks   = data.waves_peak.ECG_R_Peaks+1;

    S_Peaks   = data.waves_peak.ECG_S_Peaks+1;
    
    T_Offsets = data.waves_peak.ECG_T_Offsets+1;
    T_Onsets  = data.waves_peak.ECG_T_Onsets+1;
    T_Peaks   = data.waves_peak.ECG_T_Peaks+1;

    Peaks = R_Peaks';
    nanRows = any(isnan(R_Peaks'), 2);
    Peaks(nanRows, :) = [];
    R_Peaks = Peaks';

    nCycles = size(R_Peaks,2)-1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    for jj = 1:nCycles 
        Results(N).Study_ID = study_id;
        Results(N).Snapshot_ID = str2double(extractBetween(char(file_paths_ECG(ii)), 'ECG_LV_results_', '_h5.h5'));
        Results(N).Cycle_ID = jj;

        idx_start = R_Peaks(jj);
        idx_end   = R_Peaks(jj+1);

        Indexs = (idx_start:1:idx_end)';
        Results(N).Index = (idx_start:1:idx_end)';
        Results(N).ECG = ECG((idx_start:1:idx_end)');
        Results(N).LVP = Pressure((idx_start:1:idx_end)');

        idxs = (idx_start:1:idx_end)';
        time = linspace(0,length(idxs), length(idxs))*(1/fs);

        
        Results(N).RR_Interval = time(end)*1000;
        Results(N).HR = 60/time(end);

        P_Offset_idx = find(P_Offsets>= idx_start & P_Offsets <= idx_end);
        P_Onset_idx  = find(P_Onsets>= idx_start  & P_Onsets  <= idx_end);
        P_Peak_idx   = find(P_Peaks>= idx_start   & P_Peaks   <= idx_end);

        Q_Peak_idx   = find(Q_Peaks>= idx_start   & Q_Peaks   <= idx_end);

        S_Peak_idx   = find(S_Peaks>= idx_start   & S_Peaks   <= idx_end);

        T_Offset_idx = find(T_Offsets>= idx_start & T_Offsets <= idx_end);
        T_Onset_idx  = find(T_Onsets>= idx_start  & T_Onsets  <= idx_end);
        T_Peak_idx   = find(T_Peaks>= idx_start   & T_Peaks   <= idx_end);
        


        Results(N).R_Peak_start = R_Peaks(jj);
        Results(N).R_Peak_end = R_Peaks(jj+1);

        if isempty(P_Offset_idx)==0
            Results(N).P_Offset = P_Offsets(P_Offset_idx);
        else
            Results(N).P_Offset = 0/0;
        end
         if isempty(P_Onset_idx)==0
            Results(N).P_Onset = P_Onsets(P_Onset_idx);
         else
             Results(N).P_Onset = 0/0;
         end
         if isempty(P_Peak_idx)==0
            Results(N).P_Peak = P_Peaks(P_Peak_idx);
            L = P_Peaks(P_Peak_idx):1:idx_end;
            time_PR = linspace(0,length(L), length(L))*(1/fs);
            Results(N).RP = time_PR(end);
         else
             Results(N).P_Peak = 0/0;
             Results(N).RP = 0/0;
         end

         if isempty(T_Offset_idx)==0
            Results(N).T_Offset = T_Offsets(T_Offset_idx);
         else
             Results(N).T_Offset = 0/0;
        end
         if isempty(T_Onset_idx)==0
            Results(N).T_Onset = T_Onsets(T_Onset_idx);
         else
            Results(N).T_Onset = 0/0;
         end
         if isempty(T_Peak_idx)==0
            Results(N).T_Peak = T_Peaks(T_Peak_idx);
         else
            Results(N).T_Peak = 0/0;
         end

         if isempty(Q_Peak_idx)==0
            Results(N).Q_Peak = Q_Peaks(Q_Peak_idx);
         else
             Results(N).Q_Peak = 0/0;
         end

         if isempty(S_Peak_idx)==0
            Results(N).S_Peak = S_Peaks(S_Peak_idx);
         else
             Results(N).S_Peak = 0/0;
         end

       
        [ED_foot_detection,foot_idx] = ED_idxs_methods(Pressure((idx_start:1:idx_end)'),ECG((idx_start:1:idx_end)'),0); 
        Results(N).ED_idx = Indexs(foot_idx);
        Results(N).ED_idx_local = foot_idx;
        Results(N).ED_value = ED_foot_detection;

        LVP = Pressure((idx_start:1:idx_end)');
        Results(N).ED_ECG_value = LVP(1);


        N = N+1;
        clearvars -except ii N Results nFiles file_paths_ECG file_paths_LV ...
                    fs  main_dir data Data R_Peaks P_Offsets P_Onsets P_Peaks...
                    T_Offsets T_Onsets T_Peaks Q_Peaks S_Peaks ...
                    nCycles ECG Pressure Results_new Indexs study_id BB_Cases CIM parent_dir ALL_Results
    end
    
    clearvars -except ii N Results nFiles file_paths_ECG file_paths_LV ...
                peak_curvature_setting fs  main_dir Results_new study_id BB_Cases CIM parent_dir ALL_Results

    cd(main_dir)
    Results = [ALL_Results,Results];
end
PreResults = Results;

Results_new = Extract_Pressure_Cycle_Analysis_Results(PreResults);

function [foot,foot_idx] = ED_idxs_methods(LVP,ecg,fig_id)
fs = 240;
    [foot_, ~] = ED_balmer(LVP);
    
    if ~isnan(foot_(1))==1
        foot = LVP(foot_(1));
        foot_idx = foot_(1);
    else
        foot = 0/0;
        foot_idx = 0/0;
    end


if fig_id ==1
    time = linspace(0,length(LVP), length(LVP))*(1/fs);
    figure; hold on
    yyaxis left;
    plot(time,LVP)    
    if ~isnan(foot_(1))==1
        scatter(time(foot_(1)),LVP(foot_(1)),'r')
    end
    yyaxis right;
    plot(time,ecg)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [foot, line] = ED_balmer(sl)
    % Algorithm based off Balmer paper:
    % https://iopscience.iop.org/article/10.1088/1361-6579/aada72
    
    m = find(sl == max(sl)); m= m(1)-1;
    m_half = 2 + floor(m/2);
    
    % Supporting functions allowing for visualisation
    line = ED_balmer_shear_line(sl, 0:m_half-1, m_half);
    y_dist = line - sl(1:m_half)';
    foot = find(y_dist == max(y_dist));
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function line = ED_balmer_shear_line(sl, idx, m_half)
    line = idx*((sl(m_half)-sl(1)))/(m_half-1) + sl(1);
end

end




end


function Results_all = Extract_Pressure_Cycle_Analysis_Results(data_)

%https://www.biopac.com/application/cardiovascular-hemodynamics/advanced-feature/left-ventricular-pressure-lvp-analysis/#tabs
%https://www.datasci.com/products/software/ponemah/analysis-modules/left-ventricular-pressure
%https://www.adinstruments.com/blog/introduction-pv-loops-understanding-points-pv-loop-and-measures-cardiac-function
fs = 240;


nCycles = size(unique(cell2mat({data_.Snapshot_ID}')),1);

Results_all = [];

for MM = 1:nCycles

idxs =find(cell2mat({data_.Snapshot_ID}')==MM);
data = data_(idxs);

ED_idx_all = cell2mat({data.ED_idx}');

P_Peak_final = [];
for NN = 1:size(data,2)
  P_Peak_idxs = data(NN).P_Peak;
  if length(P_Peak_idxs) >1
    P_Peak_final =   [P_Peak_final;P_Peak_idxs(end)];
  else
      P_Peak_final =  [P_Peak_final; P_Peak_idxs(1)];
  end
end


% If we based the DS time points on P_peak-R-wave Intevaral
RP_segment = cell2mat({data.R_Peak_end}')-cell2mat({data.P_Peak}');
%RP_segment = cell2mat({data.R_Peak_end}')-P_Peak_final;


Index_all = []; ECG_all = []; LVP_all = []; 
for ii = 1:size(data,2)
    Index_all = [Index_all; data(ii).Index];
    ECG_all   = [ECG_all;   data(ii).ECG];
    LVP_all   = [LVP_all;   data(ii).LVP];

end

for jj=1:size(ED_idx_all,1)-1
    data_new(jj).Study_ID = data(1).Study_ID;
    data_new(jj).Snapshot_ID = data(1).Snapshot_ID;
    data_new(jj).Cycle_ID = jj;

    idx_start = find(Index_all==ED_idx_all(jj));
    idx_end   = find(Index_all==ED_idx_all(jj+1));

    data_new(jj).Index = Index_all(idx_start:idx_end);
    time  = linspace(0,length(idx_start:1:idx_end), length(idx_start:1:idx_end))*(1/fs);

    data_new(jj).time  = time';
    data_new(jj).ECG   = ECG_all(idx_start:idx_end);
    data_new(jj).LVP   = LVP_all(idx_start:idx_end);
    
    Index_new = data_new(jj).Index;

  

    LVP = LVP_all(idx_start:idx_end);
    dPdt = diff(LVP) ./ diff(time');

   
    LVP_max = max(LVP); data_new(jj).LVP_max = LVP_max(1);
    LVP_min = min(LVP); data_new(jj).LVP_min = LVP_min(1);


    %Developed pressure is the difference between the systolic pressure and the left ventricular end diastolic pressure (SYS-LVEDP).
    data_new(jj).DP = LVP_max(1)-LVP_min(1);

    % Calculate dP/dt
    data_new(jj).dPdt = (diff(LVP) ./ diff(time'));

    % Find maximum and minimum values of dP/dt
    dPdt_max = max(dPdt); data_new(jj).dPdt_max = dPdt_max(1);
    dPdt_min = min(dPdt); data_new(jj).dPdt_min = dPdt_min(1);
    

   % Left ventricular relaxation time constant, Tau, is the best index to evaluate left ventricular diastolic function. 
   % Tau was calculated as follows: Tau = P/(−dP/dt). Both P and (−dP/dt) can be measured on the P curve and dP/dt curve.
   % Here Tau is reported in milliseconds, milliseconds, 
   % https://doi.org/10.1016/j.ultrasmedbio.2018.03.023
    %data_new(jj).Tau = LVP_max(1)./abs(max(dPdt)); 
    data_new(jj).Tau = LVP_max(1)./abs(dPdt_min(1));


    %Contractility index is +dP/dt divided by the pressure at that point.
    data_new(jj).CI = dPdt_max(1)./LVP_max(1);


    data_new(jj).PR_Interval = data(jj).RP*1000;

    data_new(jj).RR_Interval = data(jj).RR_Interval;
    data_new(jj).ED_ED_Interval = time(end)*1000;


    %The Q-A Interval is the time in milliseconds from the start of the Q-wave,
    % in the ECG trigger channel, to the start of the systolic pressure rise (LVEDP).
    % https://www.datasci.com/products/software/ponemah/analysis-modules/left-ventricular-pressure

    data_new(jj).HR = data(jj).HR;
    data_new(jj).HR_ED = 60/time(end);

    
    PR = data(jj).RP;
    if ~isnan(PR)==1 && (PR>0) ==1
        DS_time = time(end)-PR;
        [~, DS_idx] = min(abs(time - DS_time)); 
        data_new(jj).DS_idx = DS_idx(1);
        data_new(jj).LVP_DS = LVP(DS_idx(1));
    else
       data_new(jj).DS_idx  = 0/0;
       data_new(jj).LVP_DS = 0/0;
    end

    data_new(jj).ED_idx  = 1;
    data_new(jj).LVP_ED = LVP(1);

    EDP = data(jj).LVP;
    data_new(jj).LVP_ED_Rwave = EDP(1);
    
    if ismember(data(jj).R_Peak_end, Index_new)
        data_new(jj).R_Peak = find(Index_new==data(jj).R_Peak_end,1);
    else
        data_new(jj).R_Peak = 0/0;
    end
    
    if ismember(data(jj).P_Offset, Index_new)
         data_new(jj).P_Offset = find_new(Index_new,data(jj).P_Offset);
    else
        data_new(jj).P_Offset =0/0;
    end


    if ismember(data(jj).P_Onset, Index_new)
         data_new(jj).P_Onset  = find_new(Index_new,data(jj).P_Onset);
    else
        data_new(jj).P_Onset  =0/0;
    end

    if ismember(data(jj).P_Peak, Index_new)
        data_new(jj).P_Peak   = find_new(Index_new,data(jj).P_Peak);
    else
        data_new(jj).P_Peak   =0/0;
    end
    
    if ismember(data(jj).T_Offset, Index_new)
        data_new(jj).T_Offset = find_new(Index_new,data(jj).T_Offset);
    else
        data_new(jj).T_Offset = 0/0;
    end

    if ismember(data(jj).T_Onset, Index_new)
         data_new(jj).T_Onset = find_new(Index_new,data(jj).T_Onset);
    else
        data_new(jj).T_Onset   =0/0;
    end

    if ismember(data(jj).T_Peak, Index_new)
    data_new(jj).T_Peak = find_new(Index_new,data(jj).T_Peak);
    else
        data_new(jj).T_Peak =0/0;
    end
   
    if ismember(data(jj).Q_Peak, Index_new)
        data_new(jj).Q_Peak = find_new(Index_new,data(jj).Q_Peak);
    else
        data_new(jj).Q_Peak =0/0;
    end

    if ismember(data(jj).S_Peak, Index_new)
        data_new(jj).S_Peak = find_new(Index_new,data(jj).S_Peak);
    else
        data_new(jj).S_Peak=0/0;
    end
end
Results_all = [Results_all,data_new];
end


function idx_ = find_new(Index_list,idx)
    
    idx_results = find(Index_list==idx);
    if(length(idx_results))>1
        idx_ = idx_results(end);
    else
        idx_ = idx_results(1);
    end
end

end






