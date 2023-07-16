clear; close all; clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ensure that the pressure and volume data align with the structure tree below:
% Main_Directory
%   |--->pressures (i.e., invasive raw pressure data)
%       |--->BB055 
%   |--->volumes
%       |--->beas-plsr (i.e., CIM models)
%             |--->BB055
%       |--->images (i.e., image metadata)
%             |--->BB055
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
main_path = pwd;
participant_id = 'BB055';

pressure_path = append(main_path,'\pressures\',participant_id,'\');
result_path = append(main_path,'\analyses\',participant_id,'\');
addpath(genpath(append(main_path,'\functions\')));

image_metadata_Dir = append(main_path,'\volumes\images\'); 
cim_modelDir = append(main_path,'\volumes\beas-plsr\');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                    Load CIM Volumes for Analysis                        %                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CIM_Volumes = load_cim_volume_models(image_metadata_Dir,cim_modelDir);
cd(result_path)
save('CIM_Volumes.mat','CIM_Volumes')
cd(main_path)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          Load Invasive Pressure & ECG Data for Analysis                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nPressureSnapshoots = load_invasive_raw_pressure_ecg_data(pressure_path,result_path);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                          ECG Analysis                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option A: Run ECG Anlysis through Matlab
%To call a Python function from MATLAB, you first need to ensure that you 
% have a Python distribution installed and configured correctly for MATLAB 
% to recognize it. You can do this with the pyversion command, which will 
% show you the current Python version MATLAB is using.
cd(main_path)
if count(py.sys.path,'path_to_your_python_file') == 0
    insert(py.sys.path,int32(0),'path_to_your_python_file');
end
% Import the Python module
mod = py.importlib.import_module('ECG_Analysis_Main');
% Call the function
python_path = strrep(result_path, '\', '/');
mod.ECG_analysis(python_path,participant_id ,int32(nPressureSnapshoots));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Option B: ECG_Anlysis_Main_Python in a native way inside Python
% Please ensure that the following libraries are installed:
% os, neurokit2, numpy, matplotlib, pandas, h5py, csv
% Command: python ECG_Anlysis_Main_Python.py "python_path" "participant_id" nPressureSnapshoots
%Example: python ECG_Analysis_Main.py "Z:/Pressure_Analysis/Main_Directory/analyses/" "BB055" 4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                       Pressure-volume Analysis                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[PreResults_Pressure,Results_Pressure] = Pressure_Analysis_Results (result_path,participant_id);
cd(result_path)
save('PreResults_Pressure.mat','PreResults_Pressure')
save('Results_Pressure.mat','Results_Pressure')

rr_match_threshold = 20;
LV_outliers_threshold = 1.6;
[Results_final,matched_volumes] = PV_loop_analysis (Results_Pressure,CIM_Volumes,rr_match_threshold,LV_outliers_threshold,result_path);
save('Results_final.mat','Results_final')
save('matched_volumes.mat','matched_volumes')


