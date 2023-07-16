clear; close all; clc;

main_path = pwd;
participant_id = 'BB055';
study_path = append(main_path,'\analyses\',participant_id,'\');
addpath(genpath(append(main_path,'\functions\')));

cd(study_path)
cim = importdata('CIM_Volumes.mat');
data_ed = importdata('Results_Pressure.mat');


rr_match_threshold = 20;
LV_outliers_threshold = 1.96;
[Results_final,matched_volumes] = PV_loop_analysis (data_ed,cim,rr_match_threshold,LV_outliers_threshold,study_path);
save('Results_final.mat','Results_final')
save('matched_volumes.mat','matched_volumes')