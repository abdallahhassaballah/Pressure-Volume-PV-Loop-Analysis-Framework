# Load NeuroKit and other useful packages
import os
import neurokit2 as nk
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from neurokit2.epochs import epochs_to_df
import h5py
import csv
import matplotlib
matplotlib.use('Agg')

def ECG_analysis(Path,participant_id,loop_idxs):
    ########################################################################################################################
    error_study_ID_indices = []
    error_snapshot_indices = []
    error_study_ID = []

    os.chdir(Path)
    print(Path)

    loop_idxs = int(loop_idxs)

    # load All ECG data
    for n in range(0, loop_idxs):
        # for n in [1,2]:
        try:
            print(participant_id)
            print('File = ', n)
            os.chdir(Path)
            csv_filepath = os.path.join('ECG_LV' + str(n + 1) + '.txt')
            ECG_LV_results = os.path.join('ECG_LV_results_' + str(n + 1))
            ECG = pd.read_csv(csv_filepath)['ECG']
            fs = 240
            ECG_PQRST_Processing(ECG, fs, Path, ECG_LV_results)
        except Exception as e:
            # handle the exception
            print("An exception occurred in Study ID = ", participant_id, " at snapshot = ", n, ":", e)
            error_study_ID.append(participant_id)
            error_snapshot_indices.append(n + 1)
            continue

    # name of the file
    os.chdir(Path)
    filename = "error_cases.csv"

    Cases_error = [error_study_ID_indices, error_study_ID, error_snapshot_indices]
    data = [list(row) for row in zip(*Cases_error)]
    header = ["Study_idx", "Study_ID", "Snapshot_idx"]

    # write the data to the file
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(data)
########################################################################################################################
########################################################################################################################
def ECG_PQRST_Processing(ecg_signal_org, fs, Path, filename):
    os.chdir(Path)
    Method = "neurokit"

    # Clean ECG signal
    ecg_cleaned = nk.ecg_clean(ecg_signal_org,  sampling_rate=fs, method=Method)
    ecg_signal = ecg_cleaned


    # All the heart beats
    signals, info = nk.ecg_process(ecg_signal, sampling_rate=fs, method=Method)
    epochs = nk.ecg_segment(ecg_signal, rpeaks=None, sampling_rate=fs, show=True)
    # plt.show()

    ALL_INFO = nk.ecg_process(ecg_signal, sampling_rate=fs, method=Method)

    #plt.savefig(os.path.join(Path, filename + '_ecg_segment.png'))
    #plt.close()
    #plt.clf()

    # ecg_segments = [Signal, Index, Label, Time]
    ecg_segments = np.array(epochs_to_df(epochs))
    # ecg_segments_IDs = epochs_to_df(epochs).columns.values

    # Extract R-peaks locations
    _, rpeaks = nk.ecg_peaks(ecg_signal, sampling_rate=fs)
    #_, rpeaks_org = nk.ecg_peaks(ecg_signal_org, sampling_rate=fs)

    # Visualize R-peaks in ECG signal
    nk.events_plot(rpeaks['ECG_R_Peaks'], ecg_signal)
    # plt.show()
    #plt.savefig(os.path.join(Path, filename + '_ECG_R_Peaks.png'))
    #plt.close()
    #plt.clf()

    # Delineate the ECG signal and visualizing all peaks of ECG complexes
    #_, waves_peak1 = nk.ecg_delineate(ecg_signal, rpeaks, sampling_rate=fs, method="peaks", show=True, show_type='all')
    _, waves_peak = nk.ecg_delineate(ecg_signal, rpeaks, sampling_rate=fs, method="dwt", show=True, show_type='all')

    # plt.show()
    #plt.savefig(os.path.join(Path, filename + '_waves_peak.png'))
    #plt.close()
    #plt.clf()
    ####################################################################################################################
    #                                         Save ECG_Analysis results                                                #
    ####################################################################################################################
    filepath = os.path.join(Path, filename + '_h5' + '.h5')
    print(filepath)
    with h5py.File(filepath, 'w') as hdf:
        G1 = hdf.create_group('ecg_segments')
        G1.create_dataset('Signal', data=np.array(ecg_segments[:, 0], dtype=np.float))
        G1.create_dataset('Index',  data=np.array(ecg_segments[:, 1], dtype=np.float))
        G1.create_dataset('Label',  data=np.array(ecg_segments[:, 2], dtype=np.float))
        G1.create_dataset('Time',   data=np.array(ecg_segments[:, 3], dtype=np.float))

        G2 = hdf.create_group('waves_peak')
        G2.create_dataset('ECG_P_Peaks',
                          data=np.array([list_of_values[:] for list_of_values in waves_peak.values()][0],
                                        dtype=np.float))
        G2.create_dataset('ECG_P_Onsets',
                          data=np.array([list_of_values[:] for list_of_values in waves_peak.values()][1],
                                        dtype=np.float))
        G2.create_dataset('ECG_P_Offsets',
                          data=np.array([list_of_values[:] for list_of_values in waves_peak.values()][2],
                                        dtype=np.float))
        G2.create_dataset('ECG_Q_Peaks',
                          data=np.array([list_of_values[:] for list_of_values in waves_peak.values()][3],
                                        dtype=np.float))
        G2.create_dataset('ECG_R_Peaks', data=np.array(rpeaks['ECG_R_Peaks'], dtype=np.float))
        G2.create_dataset('ECG_R_Onsets',
                          data=np.array([list_of_values[:] for list_of_values in waves_peak.values()][4],
                                        dtype=np.float))
        G2.create_dataset('ECG_R_Offsets',
                          data=np.array([list_of_values[:] for list_of_values in waves_peak.values()][5],
                                        dtype=np.float))
        G2.create_dataset('ECG_S_Peaks',
                          data=np.array([list_of_values[:] for list_of_values in waves_peak.values()][6],
                                        dtype=np.float))
        G2.create_dataset('ECG_T_Peaks',
                          data=np.array([list_of_values[:] for list_of_values in waves_peak.values()][7],
                                        dtype=np.float))
        G2.create_dataset('ECG_T_Onsets',
                          data=np.array([list_of_values[:] for list_of_values in waves_peak.values()][8],
                                        dtype=np.float))
        G2.create_dataset('ECG_T_Offsets',
                          data=np.array([list_of_values[:] for list_of_values in waves_peak.values()][9],
                                        dtype=np.float))

        G3 = hdf.create_group('ecg_quality')
        G3.create_dataset('ECG_Quality', data=np.array(signals.ECG_Quality, dtype=np.float))

        G4 = hdf.create_group('ECG_info')

        G4.create_dataset('ECG_Phase_Atrial', data=np.array(ALL_INFO[0].ECG_Phase_Atrial), dtype=np.float)
        G4.create_dataset('ECG_Phase_Completion_Atrial', data=np.array(ALL_INFO[0].ECG_Phase_Completion_Atrial), dtype=np.float)
        G4.create_dataset('ECG_Phase_Completion_Ventricular', data=np.array(ALL_INFO[0].ECG_Phase_Completion_Ventricular), dtype=np.float)
        G4.create_dataset('ECG_Phase_Ventricular', data=np.array(ALL_INFO[0].ECG_Phase_Ventricular), dtype=np.float)
        G4.create_dataset('ECG_Rate', data=np.array(ALL_INFO[0].ECG_Rate), dtype=np.float)

        G5 = hdf.create_group('ECG')
        G5.create_dataset('ECG_Raw', data=np.array(ecg_signal_org, dtype=np.float))
        G5.create_dataset('ECG_Cleaned', data=np.array(ecg_cleaned, dtype=np.float))



########################################################################################################################
#Path_= os.path.join('../Main_Directory/analyses')
#participant_id = 'BB055'
#loop_idxs = 4
#ECG_analysis(Path_,participant_id,loop_idxs)

