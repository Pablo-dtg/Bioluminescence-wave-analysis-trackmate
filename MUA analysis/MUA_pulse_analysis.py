
"""
Created on Wed Jul 12 14:33:26 2023

@author: pdetorresgutierrez
"""

import matplotlib.pyplot as plt
from matplotlib.figure import figaspect
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from scipy.signal import butter, filtfilt
from scipy.signal import medfilt
from scipy.ndimage import gaussian_filter
import datetime



#location of the file
dirc = r'\\vf-mcb-circadian-data\mcb-circadian-data$\Neurofys\Karlijn\MUA data'
#name of the file
exp = '2025_08_19_10_50_19__df7_' 
experiment_start_time = '2025-08-19 10:50:19' 
#output files location
output_dirc = r'//vf-mcb-circadian-data/mcb-circadian-data$/Neurofys/Karlijn/Script/Output' 
#for data slicing
cut_point_start =2500
# cut_point_end = 18500


df = pd.read_excel(dirc + '/' + exp + '.xlsx')
pulsetimes = pd.read_csv(dirc + '/' + exp + '.txt')
pulsetimes = pulsetimes['Pulses'].tolist()
#data refinement
df = df.iloc[cut_point_start:]
# df = df.iloc[:cut_point_end]

# y1 = df.iloc[:,]                
# y1 = y1.values.tolist()

time = df.iloc[:,0].tolist()
ST1 = df.iloc[:,2]
ST2 = df.iloc[:,3]



experiment_start_dt = datetime.datetime.strptime(experiment_start_time, "%Y-%m-%d %H:%M:%S")

time_list = [experiment_start_dt + datetime.timedelta(seconds=(t.hour * 3600 + t.minute * 60 + t.second)) for t in time]

time_s =[
    t.hour * 3600 + t.minute * 60 + t.second for t in time
]


# =============================================================================
# # FILTERING
# =============================================================================
# def bandpass_filter(data, lowcut=1, highcut=400, fs=10000, order=2):
#     nyquist = 0.5 * fs
#     low = lowcut / nyquist
#     high = highcut / nyquist
#     b, a = butter(order, [low, high], btype='band')
#     return filtfilt(b, a, data)
# ST1_filter = bandpass_filter(ST1)

# =============================================================================
# #LOW PASS FILTER
# =============================================================================
# Filter parameters
cutoff_freq = 0.05  # Cutoff frequency in Hz
fs = 1  # Sampling rate in Hz
order = 4  # Filter order

# Design the low-pass filter
nyq = 0.5 * fs
normal_cutoff = cutoff_freq / nyq
b, a = butter(order, normal_cutoff, btype='lowpass')

#smoothing 
ST1_smoothed = medfilt(ST1, kernel_size=3)
ST2_smoothed = medfilt(ST2, kernel_size=3)

ST1_filter = filtfilt(b, a, ST1_smoothed)
ST2_filter = filtfilt(b, a, ST2_smoothed)

# =============================================================================
# #GAUSSIAN FILTER
# =============================================================================
# ST1_filter = gaussian_filter(ST1, sigma=4)
# ST2_filter = gaussian_filter(ST2, sigma=4)



time_s =[
    t.hour * 3600 + t.minute * 60 + t.second for t in time
]

#PULSE INFORMATION

hik1_s = pulsetimes[0]
hik1_e = pulsetimes[1]


hik8_s = pulsetimes[2]
hik8_e = pulsetimes[3]

hik11_s = pulsetimes[4]
hik11_e = pulsetimes[5]

hik15_s = pulsetimes[6]
hik15_e = pulsetimes[7]

# hik21_s = pulsetimes[8]
# hik21_e = pulsetimes[9]

# hik15_start_s = pulsetimes[10]
# hik15_start_e = pulsetimes[11]

#nim_s = pulsetimes[12]
#nim_e = pulsetimes[13]

# =============================================================================
# QUANTIFICATION
# =============================================================================

#Baseline calculation

# hik15_start_s_dt = datetime.datetime.strptime(hik15_start_s, '%Y-%m-%d %H:%M')
# if hik15_start_s_dt in time_list:
#     index_hik15_start_s_dt = time_list.index(hik15_start_s_dt)
#     print("Index of event time:", index_hik15_start_s_dt)
# else:
#     print("Event time not found in the list")

# bl_21b = np.mean(ST1_filter[index_hik15_start_s_dt-600:index_hik15_start_s_dt]) 

  

#first peak 15mM (pulsetimes[10 and 11])
# int_hik15_start_s = ST1_filter[index_hik15_start_s_dt:index_hik15_start_s_dt+2100] #35 minutes after the start of the pulse
# max_hik15_start_s = max(int_hik15_start_s)    #find the maximum in that interval
# rel_max_hik15_start_s = max_hik15_start_s/bl_21b    #dividing maximum and baseline

# int_hik15_start_s = ST1_filter[index_hik15_start_s_dt:index_hik15_start_s_dt+2100] #35 minutes after the start of the pulse
# min_hik15_start_s = min(int_hik15_start_s)    #find the maximum in that interval
# rel_min_hik15_start_s = min_hik15_start_s/bl_21b    #dividing maximum and baseline


#peak 1mM
hik1_s_dt = datetime.datetime.strptime(hik1_s, '%Y-%m-%d %H:%M')  #convert str to datetime

if hik1_s_dt in time_list:              #calculating the index for when the pulse starts
    index_hik1_s_dt = time_list.index(hik1_s_dt)
    print("Index of event time:", index_hik1_s_dt)
else:
    print("Event time not found in the list")
    
bl_1 = np.mean(ST1_filter[index_hik1_s_dt-600:index_hik1_s_dt]) 

int_hik1_s = ST1_filter[index_hik1_s_dt:index_hik1_s_dt+2100]
min_int_hik1_s = min(int_hik1_s)    #find the maximum in that interval
rel_min_int_hik1_s = min_int_hik1_s/bl_1


#peak 8mM
hik8_s_dt = datetime.datetime.strptime(hik8_s, '%Y-%m-%d %H:%M')  #convert str to datetime

if hik8_s_dt in time_list:              #calculating the index for when the pulse starts
    index_hik8_s_dt = time_list.index(hik8_s_dt)
    print("Index of event time:", index_hik8_s_dt)
else:
    print("Event time not found in the list")

bl_8 = np.mean(ST1_filter[index_hik8_s_dt-600:index_hik8_s_dt]) 

int_hik8_s = ST1_filter[index_hik8_s_dt:index_hik8_s_dt+2100]
max_int_hik8_s = max(int_hik8_s)    #find the maximum in that interval
rel_max_int_hik8_s = max_int_hik8_s/bl_8

#peak 11mM
hik11_s_dt = datetime.datetime.strptime(hik11_s, '%Y-%m-%d %H:%M')  #convert str to datetime

if hik11_s_dt in time_list:              #calculating the index for when the pulse starts
    index_hik11_s_dt = time_list.index(hik11_s_dt)
    print("Index of event time:", index_hik11_s_dt)
else:
    print("Event time not found in the list")

bl_11 = np.mean(ST1_filter[index_hik11_s_dt-600:index_hik11_s_dt])

int_hik11_s = ST1_filter[index_hik11_s_dt:index_hik11_s_dt+2100]
max_int_hik11_s = max(int_hik11_s)    #find the maximum in that interval
rel_max_int_hik11_s = max_int_hik11_s/bl_11

#peak 15mM
hik15_s_dt = datetime.datetime.strptime(hik15_s, '%Y-%m-%d %H:%M')  #convert str to datetime

if hik15_s_dt in time_list:              #calculating the index for when the pulse starts
    index_hik15_s_dt = time_list.index(hik15_s_dt)
    print("Index of event time:", index_hik11_s_dt)
else:
    print("Event time not found in the list")
    
bl_15 = np.mean(ST1_filter[index_hik15_s_dt-600:index_hik15_s_dt])

int_hik15_s = ST1_filter[index_hik15_s_dt:index_hik15_s_dt+2100]
max_int_hik15_s = max(int_hik15_s)    #find the maximum in that interval
rel_max_int_hik15_s = max_int_hik15_s/bl_15


min_int_hik15_s = min(int_hik15_s)    #find the minimum in that interval
rel_min_int_hik15_s = min_int_hik15_s/bl_15

#peak 21mM trough
# hik21_s_dt = datetime.datetime.strptime(hik21_s, '%Y-%m-%d %H:%M')  #convert str to datetime

# if hik21_s_dt in time_list:              #calculating the index for when the pulse starts
#     index_hik21_s_dt = time_list.index(hik21_s_dt)
#     print("Index of event time:", index_hik21_s_dt)
# else:
#     print("Event time not found in the list")

# bl_21 = np.mean(ST1_filter[index_hik21_s_dt-600:index_hik21_s_dt])

# int_hik21_s = ST1_filter[index_hik21_s_dt:index_hik21_s_dt+2100]
# min_int_hik21_s = min(int_hik21_s)    #find the maximum in that interval
# rel_min_int_hik21_s = min_int_hik21_s/bl_21

# int_hik21_s = ST1_filter[index_hik21_s_dt:index_hik21_s_dt+2100]
# max_int_hik21_s = max(int_hik21_s)    #find the maximum in that interval
# rel_max_int_hik21_s = max_int_hik21_s/bl_21

#Saving the values in an excel file
df_excel = pd.DataFrame()

pulselist = ['1','11','15','21']
baselinelist = [bl_1, bl_8, bl_11, bl_15]
maxlist = ['-', max_int_hik8_s, max_int_hik11_s, max_int_hik15_s]
minlist = [ min_int_hik1_s, '-', '-', '-']
relmaxlist = ['-', rel_max_int_hik8_s, rel_max_int_hik11_s, rel_max_int_hik15_s]
relminlist = [ rel_min_int_hik1_s, '-', '-', rel_min_int_hik15_s]

df_excel['Pulse K+ (mM)'] = pulselist
df_excel['Baseline'] = baselinelist
df_excel['Maximum'] = maxlist
df_excel['Minimum'] = minlist
df_excel['Relative maximum'] = relmaxlist
df_excel['Relative minimum'] = relminlist


df_excel.to_string(index=False)

with pd.ExcelWriter(output_dirc + '/peaks_'+ exp + '.xlsx') as writer:
    df_excel.to_excel(writer, sheet_name='Peaks')





# =============================================================================
# #FIGURES
# =============================================================================

leg_pos = 1700   #position of the legend for the K+ concentrations


w, h = figaspect(1/3)
fig, ax = plt.subplots(figsize=(w,h))
ax.plot(time_list, ST1, c='gray', linewidth=1)
ax.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M:%S'))
# ax.set_ylim([0, 200])
ax.set_xlabel('Time')
ax.set_ylabel('MUA counts', color='black')
ax.axvspan(hik1_s, hik1_e, alpha=0.5, color='red')
ax.text(datetime.datetime.strptime(hik1_s, '%Y-%m-%d %H:%M'), leg_pos,'1 mM', ha='left', va='bottom', fontsize= '14')
ax.axvspan(hik8_s, hik8_e, alpha=0.5, color='red')
ax.text(datetime.datetime.strptime(hik8_s, '%Y-%m-%d %H:%M'), leg_pos,'8 mM', ha='left', va='bottom', fontsize= '14')
ax.axvspan(hik11_s, hik11_e, alpha=0.5, color='red')
ax.text(datetime.datetime.strptime(hik11_s, '%Y-%m-%d %H:%M'), leg_pos,'11 mM', ha='left', va='bottom', fontsize= '14')
ax.axvspan(hik15_s, hik15_e, alpha=0.5, color='red')
ax.text(datetime.datetime.strptime(hik15_s, '%Y-%m-%d %H:%M'), leg_pos,'15 mM', ha='left', va='bottom', fontsize= '14')
# ax.axvspan(hik21_s, hik21_e, alpha=0.5, color='red')
# ax.text(datetime.datetime.strptime(hik21_s, '%Y-%m-%d %H:%M'), leg_pos,'21 mM', ha='left', va='bottom', fontsize= '14')
# ax.axvspan(hik15_start_s, hik15_start_e, alpha=0.5, color='red')
# ax.text(datetime.datetime.strptime(hik15_start_s, '%Y-%m-%d %H:%M'), leg_pos,'15 mM', ha='left', va='bottom', fontsize= '14')
# ax.axvspan(nim_s, nim_e, alpha=0.2, color='blue')
# ax.text(datetime.datetime.strptime(nim_s, '%Y-%m-%d %H:%M'), leg_pos,'nim', ha='left', va='bottom', fontsize= '14')
plt.title("ST1")
plt.savefig(output_dirc + '\ST1',dpi=300)
plt.show()

fig, ax = plt.subplots(figsize=(w,h))
ax.plot(time_list, ST1_filter, c='black')
ax.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M'))
# ax.set_ylim([0, 200])
ax.set_xlabel('Time (h)')
ax.set_ylabel('MUA counts', color='black')
# plt.axhline(y=bl, color='b', linestyle='--')    #plotting the baseline throughout the experiment
ax.axvspan(hik1_s, hik1_e, alpha=0.5, color='red')
ax.text(datetime.datetime.strptime(hik1_s, '%Y-%m-%d %H:%M'), leg_pos,'1 mM', ha='left', va='bottom', fontsize= '14')
ax.axvspan(hik8_s, hik8_e, alpha=0.5, color='red')
ax.text(datetime.datetime.strptime(hik8_s, '%Y-%m-%d %H:%M'), leg_pos,'8 mM', ha='left', va='bottom', fontsize= '14')
ax.axvspan(hik11_s, hik11_e, alpha=0.5, color='red')
ax.text(datetime.datetime.strptime(hik11_s, '%Y-%m-%d %H:%M'), leg_pos,'11 mM', ha='left', va='bottom', fontsize= '14')
ax.axvspan(hik15_s, hik15_e, alpha=0.5, color='red')
ax.text(datetime.datetime.strptime(hik15_s, '%Y-%m-%d %H:%M'), leg_pos,'15 mM', ha='left', va='bottom', fontsize= '14')
# ax.axvspan(hik21_s, hik21_e, alpha=0.5, color='red')
# ax.text(datetime.datetime.strptime(hik21_s, '%Y-%m-%d %H:%M'), leg_pos,'21 mM', ha='left', va='bottom', fontsize= '14')
# ax.axvspan(hik15_start_s, hik15_start_e, alpha=0.5, color='red')
# ax.text(datetime.datetime.strptime(hik15_start_s, '%Y-%m-%d %H:%M'), leg_pos,'15 mM', ha='left', va='bottom', fontsize= '14')

# ax.axvspan(nim_s, nim_e, alpha=0.2, color='blue')
# ax.text(datetime.datetime.strptime(nim_s, '%Y-%m-%d %H:%M'), leg_pos,'nim', ha='left', va='bottom', fontsize= '14')


##ax.axhline(y=bl_1, color='blue', linestyle='--', label='Baseline 1 mM')
##ax.axhline(y=bl_8, color='green', linestyle='--', label='Baseline 8 mM')
# # Plot 10-minute baseline segment before 1 mM pulse
# baseline1_times = time_list[index_hik1_s_dt - 600:index_hik1_s_dt]
# baseline1_values = [bl_1] * len(baseline1_times)
# ax.plot(baseline1_times, baseline1_values, color='blue', linestyle='solid',  linewidth = 3, alpha = 0.5)

# # Plot 10-minute baseline segment before 8 mM pulse
# baseline8_times = time_list[index_hik8_s_dt - 600:index_hik8_s_dt]
# baseline8_values = [bl_8] * len(baseline8_times)
# ax.plot(baseline8_times, baseline8_values, color='blue', linestyle='solid',  linewidth = 3, alpha = 0.5)

# # Plot 10-minute baseline segment before 11 mM pulse
# baseline11_times = time_list[index_hik11_s_dt - 600:index_hik11_s_dt]
# baseline11_values = [bl_11] * len(baseline11_times)
# ax.plot(baseline11_times, baseline11_values, color='blue', linestyle='solid',  linewidth = 3, alpha = 0.5)

# # Plot 10-minute baseline segment before 15 mM pulse
# baseline15_times = time_list[index_hik15_s_dt - 600:index_hik15_s_dt]
# baseline15_values = [bl_15] * len(baseline15_times)
# ax.plot(baseline15_times, baseline15_values, color='blue', linestyle='solid',  linewidth = 3, alpha = 0.5)

# Plot 10-minute baseline segment before 21 mM pulse
# baseline21_times = time_list[index_hik21_s_dt - 600:index_hik21_s_dt]
# baseline21_values = [bl_21] * len(baseline21_times)
# ax.plot(baseline21_times, baseline21_values,color='blue', linestyle='solid',  linewidth = 3, alpha = 0.5)

# Plot 10-minute baseline segment before starting 15/21 mM pulse
# baseline21b_times = time_list[index_hik15_start_s_dt - 600:index_hik15_start_s_dt]
# baseline21b_values = [bl_21b] * len(baseline21b_times)
# ax.plot(baseline21b_times, baseline21b_values, color='blue', linestyle='solid',  linewidth = 3, alpha = 0.5)

ax.legend(loc="best")
# plt.title("ST1 filtered")
plt.savefig(output_dirc+'\ST1 filter.eps',dpi=300)
plt.show()

fig, ax = plt.subplots(figsize=(w,h))
ax.plot(time_list, ST2, c='gray')
ax.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M:%S'))
ax.set_xlabel('Time (s)')
ax.set_ylabel('MUA counts', color='grey')
ax.axvspan(hik1_s, hik1_e, alpha=0.5, color='red')
ax.text(datetime.datetime.strptime(hik1_s, '%Y-%m-%d %H:%M'), leg_pos,'1 mM', ha='left', va='bottom', fontsize= '14')
ax.axvspan(hik8_s, hik8_e, alpha=0.5, color='red')
ax.text(datetime.datetime.strptime(hik8_s, '%Y-%m-%d %H:%M'), leg_pos,'8 mM', ha='left', va='bottom', fontsize= '14')
ax.axvspan(hik11_s, hik11_e, alpha=0.5, color='red')
ax.text(datetime.datetime.strptime(hik11_s, '%Y-%m-%d %H:%M'), leg_pos,'11 mM', ha='left', va='bottom', fontsize= '14')
ax.axvspan(hik15_s, hik15_e, alpha=0.5, color='red')
ax.text(datetime.datetime.strptime(hik15_s, '%Y-%m-%d %H:%M'), leg_pos,'15 mM', ha='left', va='bottom', fontsize= '14')
# ax.axvspan(hik21_s, hik21_e, alpha=0.5, color='red')
# ax.text(datetime.datetime.strptime(hik21_s, '%Y-%m-%d %H:%M'), leg_pos,'21 mM', ha='left', va='bottom', fontsize= '14')
# ax.axvspan(hik15_start_s, hik15_start_e, alpha=0.5, color='red')
# ax.text(datetime.datetime.strptime(hik15_start_s, '%Y-%m-%d %H:%M'), leg_pos,'21 mM', ha='left', va='bottom', fontsize= '14')
# ax.axvspan(nim_s, nim_e, alpha=0.2, color='blue')
# ax.text(datetime.datetime.strptime(nim_s, '%Y-%m-%d %H:%M'), leg_pos,'nim', ha='left', va='bottom', fontsize= '14')
plt.legend(loc="best")
plt.title("ST2")
plt.savefig(output_dirc + '\ST2',dpi=300)
plt.show()

fig, ax = plt.subplots(figsize=(w,h))
ax.plot(time_list, ST2_filter, c='black')
ax.xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%H:%M:%S'))
ax.set_xlabel('Time (s)')
ax.set_ylabel('MUA counts', color='black')
ax.axvspan(hik1_s, hik1_e, alpha=0.5, color='red')
ax.text(datetime.datetime.strptime(hik1_s, '%Y-%m-%d %H:%M'), leg_pos,'1 mM', ha='left', va='bottom', fontsize= '14')
ax.axvspan(hik8_s, hik8_e, alpha=0.5, color='red')
ax.text(datetime.datetime.strptime(hik8_s, '%Y-%m-%d %H:%M'), leg_pos,'8 mM', ha='left', va='bottom', fontsize= '14')
ax.axvspan(hik11_s, hik11_e, alpha=0.5, color='red')
ax.text(datetime.datetime.strptime(hik11_s, '%Y-%m-%d %H:%M'), leg_pos,'11 mM', ha='left', va='bottom', fontsize= '14')
ax.axvspan(hik15_s, hik15_e, alpha=0.5, color='red')
ax.text(datetime.datetime.strptime(hik15_s, '%Y-%m-%d %H:%M'), leg_pos,'15 mM', ha='left', va='bottom', fontsize= '14')
# ax.axvspan(hik21_s, hik21_e, alpha=0.5, color='red')
# ax.text(datetime.datetime.strptime(hik21_s, '%Y-%m-%d %H:%M'), leg_pos,'21 mM', ha='left', va='bottom', fontsize= '14')
# ax.axvspan(hik15_start_s, hik15_start_e, alpha=0.5, color='red')
# ax.text(datetime.datetime.strptime(hik15_start_s, '%Y-%m-%d %H:%M'), leg_pos,'21 mM', ha='left', va='bottom', fontsize= '14')
# ax.axvspan(nim_s, nim_e, alpha=0.2, color='blue')
# ax.text(datetime.datetime.strptime(nim_s, '%Y-%m-%d %H:%M'), leg_pos,'nim', ha='left', va='bottom', fontsize= '14')
plt.legend(loc="best")
plt.title("ST2 filtered")
plt.savefig(output_dirc + '\ST2 filter.eps',dpi=300)
plt.show()

# with pd.ExcelWriter(r"\\vf-mcb-circadian-data\mcb-circadian-data$\Neurofys\Pablo\EI and period\EI balance MUA\Script\output.xlsx") as writer:        #saving in a new document in different data sheets
#     res.to_excel(writer) #corrected 340 (minus bg)