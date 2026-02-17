# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 15:39:15 2023

@author: pdetorresgutierrez
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Dec 18 16:58:22 2023

@author: pdetorresgutierrez
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.figure import figaspect
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['pdf.use14corefonts'] = True
plt.rcParams['ps.useafm'] = True
plt.rcParams['text.usetex'] = False
import numpy
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
from scipy import signal
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
from scipy.signal import savgol_filter
from datetime import datetime, timedelta
import math
import statistics
from matplotlib import cm
from matplotlib.ticker import  MultipleLocator
import os

path = os.getcwd() #reads directory
print(path)

figures_on = False  #toggle to visualise processing of individual traces
exp = 'export_260211_2'
experiment_start_time = '2024-08-12 12:39:00'  # Example format: 'YYYY-MM-DD HH:MM:SS'



dic_traces = {}   #x axis
dic_times = {}      #y axis
dic_t_length = {}    #not used yet
dic_t_id= {}
dic_pos_x = {}
dic_pos_x_up = {}  # updated later to have the correct keys
dic_pos_y = {}
dic_pos_y_up = {}
dic_disp_x = {}  #displacement along the x axis from first to last spots of the trace
dic_disp_y = {} #displacement along the y axis from first to last spots of the trace
dic_tot_dist = {} #dictionary containing the distance for every trace
dic_dy_dx = {} #firds derivative
dic_d2y_dx = {}
dic_x_fine ={} #first derivative x dat for plotting
df = pd.read_csv(r'\\vf-mcb-circadian-data\mcb-circadian-data$\Neurofys\Pablo\Membrane excitability and period\Bioluminiscence PER2 data\Script' + '/'  + exp + '.csv')



df = df.iloc[3:]

del df['MANUAL_SPOT_COLOR']
df = df.dropna(axis = 0)
df['TRACK_ID'] = df['TRACK_ID'].astype(float)
df_tracks = df.groupby('TRACK_ID')

dfs = [df_tracks.get_group(x) for x in df_tracks.groups]

i=0
for a in dfs:
    if a.shape[0] <= 72:                  #important to filter traces that arent long enough
        continue
    else:
    
        trace = a['MEAN_INTENSITY_CH1'].astype(float).to_numpy()
        t = a['POSITION_T'].astype(float).to_numpy()
        pos_x = (a['POSITION_X'].astype(float).to_numpy())/96         #I divided by 96 because 1px is approx 1/96th of an inch
        pos_y = (-1*a['POSITION_Y'].astype(float).to_numpy())/96
        dic_t_id[i] = a.iloc[0,2]
        # if t[0] > 12:                                                        #Maitreyi script
        #     continue
        dic_traces[i] = trace
        dic_times[i] = t
        dic_pos_x[i] = pos_x
        dic_pos_y[i] = pos_y
        if figures_on == True:
            fig, ax = plt.subplots()
            ax.plot(t, trace, 'b-')
            ax.set(title='Track '+str(a.iloc[0,2]))
            print(str(a.iloc[0,2]))
        
        plt.show()
        i = i+1

for key in dic_pos_x.keys():
    pos_x_start = dic_pos_x[key][0]
    pos_x_end = dic_pos_x[key] [-1]
    pos_y_start = dic_pos_y[key][0]
    pos_y_end = dic_pos_y[key] [-1]  
    dic_disp_x[key] = pos_x_end - pos_x_start
    dic_disp_y[key] = pos_y_end - pos_y_start
    
df_c_y = pd.DataFrame()   #all cells y values
df_c_x = pd.DataFrame()   #all cells x values

if figures_on == True:
    fig, ax2 = plt.subplots()
    for e in range(0,i):
        ax2.plot(dic_times[e],dic_traces[e],linewidth= 0.25)
    
    plt.savefig(r'\\vf-mcb-circadian-data\mcb-circadian-data$\Neurofys\Pablo\Membrane excitability and period\Bioluminiscence PER2 data\Script\traces', dpi=300) 
    plt.show() 


for trac in range(0,i):
    # df_c['Track_'+str(dic_t_id[trac])] = dic_times[trac]
    df_c_y = pd.concat([df_c_y,pd.Series(dic_traces[trac])], ignore_index = True, axis =1)
    df_c_x = pd.concat([df_c_x,pd.Series(dic_times[trac])], ignore_index = True, axis =1) 

#trace
def calculate_distance(x1, y1, x2, y2):
    return math.sqrt((x2 - x1)**2 + (y2 - y1)**2)

def calculate_total_distance(x_coordinates, y_coordinates):
    total_distance = 0.0
    for i in range(1, len(x_coordinates)):
        distance = calculate_distance(x_coordinates[i-1], y_coordinates[i-1], x_coordinates[i], y_coordinates[i])
        total_distance += distance
    return total_distance

for p in dic_pos_x.keys():
    dic_tot_dist[p] = calculate_total_distance(dic_pos_x[p], dic_pos_y[p])
    


# =============================================================================
# 1 by 1 analysis for each cell trace
# =============================================================================
df_details = pd.read_excel(r'\\vf-mcb-circadian-data\mcb-circadian-data$\Neurofys\Pablo\Membrane excitability and period\Bioluminiscence PER2 data\Script\exp_details.xlsx')

kcl = df_details['KCl_conc'].tolist()
times = df_details['start_time'].tolist()
year = df_details['year'].tolist()
month = df_details['month'].tolist()
day = df_details['day'].tolist()
hour = df_details['hour'].tolist()
minute = df_details['minute'].tolist()

dic_y = {}      #making a dictionary with single points
dic_fun = {}       #dictionary for functions
dic_x = {}      #dictionary for x axes with the correct times
dic_x_h = {}    #same information as dic_x but instead of having the times in datetime format its in hours
dic_x_h_lin = {}    #x in hours linspace
dic_y_raw = {}
dic_y_func = {}       #functions dictionary
dic_y_func_n = {}    #same but normalised values
dic_y_filter = {}   #filtered y data
dic_x_smooth = {}
dic_y_smooth = {}     #smoothed data
dic_param_all = {}         #all parameters from the fitting
dic_cyc_per = {}   #dictionary with each trace's cycle periods
dic_cyc_per_avrg = {}
dic_index = {}   #dic with all the indexes of cells not excluded
#dictionaries for the raw data
dic_r_y_filter = {}
dic_r_x_smooth = {}
dic_r_y_smooth = {}
dic_r_y = {}
dic_extrema= {}     #Contains x and y coordinates of maximum and minimum for each trace
dic_amplitude = {}
dic_amplitude_cyc = {} #amplitudes per cycle
dic_peaks = {}
dic_peaks_x ={}
dic_troughs = {}
dic_single_peaktrough_x = {}
dic_r_points_x = {} #to index variables to their correct cycle
dic_cyc_index = {}   #contains info about which cycle the parameters calculated refer to
dic_SE = {}  #contains the standard error of the regression for every cell



for index in df_c_y:
    y1 = df_c_y[index].values.astype(float)
    x1 = df_c_x[index].values.astype(float)
    #remove Nan
    y1 = y1[~np.isnan(y1)]
    x1 = x1[~np.isnan(x1)]
      

    
    # Parse the experiment start time into a datetime object
    experiment_start_dt = datetime.strptime(experiment_start_time, '%Y-%m-%d %H:%M:%S')
    
    
    timepoints = []
    timepoints_h = []
    time_interval = timedelta(hours=1)
    
    # Starting point for time
    time_from_start = timedelta(hours=x1[0] - 1)  # Adjust for the first time interval
    timepoint = experiment_start_dt + time_from_start  # Start from the experiment start time
    
    # Generate timepoints based on the intervals in x1 (this replaces the need for exp_time)
    for n in range(len(x1)):        
        timepoint = timepoint + time_interval  # Add 1 hour for each step
        timepoints.append(timepoint)
    
    dic_x[index] = timepoints
    
    # Calculate time in hours from the experiment start time
    for o in timepoints:
        time_h = (o - experiment_start_dt).total_seconds() / 3600  # Convert time difference to hours
        timepoints_h.append(time_h)
    
    dic_x_h[index] = np.array(timepoints_h)
    
    timepoints_h = np.array(timepoints_h)


    # =============================================================================
    #    BINNING
    # =============================================================================
    # n_points = len(y1)
    # bin_size = 6      #in per2 no need to bin
    # # remainder = n_points % bin_size
    # # n_points = n_points - remainder


    # n_bin = int(n_points/bin_size)
    # binned_y = []
    # split = np.array_split(y1, n_bin)   #some splits longer than bin size e.g. 61 (61)

    # for e in split:
    #     bin_mean = np.mean(e)
    #     binned_y.append(bin_mean)
    
    #for per2 no binning and x1 determined by the dictionary
    
    # # x1 = np.arange(0,n_points,bin_size)
    # x1 = np.arange(0,len(binned_y))

    # #delete initial parts
    # first_value =24                      #GIVE IN SECONDS, points of data you wish to exclude
    # first_value = int(first_value/bin_size)                       #transformation of that to hours
    # x1 = (x1[first_value:]-first_value+1)                         #slices

    # # in per2 already hours no need to change
    # # x1= x1 *(bin_size/3600)                               # standardize so it shows the results in hours otherwise the parameters change


    binned_y= y1



    # =============================================================================
    #     DETRENDING
    # =============================================================================
    # det_y = signal.detrend( , type='linear')
    #detrend2             How to know which one is best?? Use Akaike Information Criterion

    k = 4
    model = np.polyfit(x1, binned_y, k)
    predicted = np.polyval(model, x1)
    if math.isnan(predicted[0]) :
        k = 5
        model = np.polyfit(x1, binned_y, k)
        predicted = np.polyval(model, x1)


    if figures_on == True:
        fig, axes = plt.subplots(nrows=2, sharex=True)
        axes[0].plot(x1, binned_y, 'ro', markersize='4')
        axes[0].plot(x1, predicted, 'k-')
        axes[0].set(title='Original Data and Polynomial Trend')
    
        axes[1].plot(x1, binned_y - predicted, 'ro', markersize='4')
        axes[1].set(title='Detrended Residual')
        plt.savefig(r'\\vf-mcb-circadian-data\mcb-circadian-data$\Neurofys\Pablo\Membrane excitability and period\EI balance MUA\Script\Ouput graphs\Detrend',dpi=300)
        plt.show()

    det_y = binned_y-predicted
    



    # Need to transform x to hours (standardize) when subtracting the fitting, because the parameters of the fit are dependent on the number of timepoints
          
    # =============================================================================
    #       CURVE FITTING
    # =============================================================================


    A_guess= (abs(np.max(det_y))+abs(np.min(det_y)))/2
    guess_std = 3*np.std(det_y)/(2**0.5)/(2**0.5)
    phi_guess=np.pi*(0.2)
    y0_guess=np.mean(det_y)
    f_guess= 0.26
    guess = np.array([y0_guess, A_guess, phi_guess])
    lam_guess = -0.01


    data_first_guess = guess_std*np.exp(timepoints_h*lam_guess)*np.cos(f_guess*timepoints_h+phi_guess) + y0_guess

    # Define the function to optimize, in this case, we want to minimize the difference
    # between the actual data and our "guessed" parameters
    optimize_func = lambda n: n[0]*np.exp(timepoints_h*n[1])*np.cos(n[2]*timepoints_h+n[3]) + n[4] - det_y
    A, lam, f, phi, y0 = leastsq(optimize_func, [A_guess, lam_guess, f_guess, phi_guess, y0_guess])[0]

    T = abs(2*np.pi/(f*1))   #change the number times f depending on bin size


   
    # # recreate the fitted curve using the optimized parameters
    # data_fit = est_amp*np.sin(est_freq*t+est_phase) + est_mean
    
    dic_param = {}
    dic_param['A'] = A
    dic_param['lambda'] = lam
    dic_param['f'] = f
    dic_param['phi'] = phi
    dic_param['y0'] = y0
    dic_param['T'] = T
    dic_param_all[index] = dic_param    #one single dictionary with all stuff 
    print(index)
    
   
        
    # check the displacement if beyond threshold then eliminate trace
    # disp_threshold = 0.000001   #0.000001
    # if abs(dic_disp_x[index]) < disp_threshold and abs(dic_disp_y[index]) < disp_threshold:
    #     dic_param.popitem()
    #     dic_param_all.popitem()
    #     dic_times.popitem()
    #     dic_traces.popitem()
    #     dic_x.popitem()
    #     dic_x_h.popitem()
    #     # dic_t_id.popitem()
    #     dic_pos_x.pop(index)
    #     dic_pos_y.pop(index)
    #     dic_tot_dist.pop(index)
    #     print('insufficient displacement')
        # continue
    # if abs(dic_tot_dist[index]) <= 0.00004: #0.00004 #To remove cells that do not move too much (possible correction for dead pixels)
    #     dic_param.popitem()
    #     dic_param_all.popitem()
    #     dic_times.popitem()
    #     dic_traces.popitem()
    #     dic_x.popitem()
    #     dic_x_h.popitem()
    #     # dic_t_id.popitem()
    #     dic_pos_x.pop(index)
    #     dic_pos_y.pop(index)
    #     dic_tot_dist.pop(index)
    #     print('insufficient total distance')
        # continue
    if T < 20 or T > 28:           #if it cant fit a cosinor skip this trace
        dic_param.popitem()
        dic_param_all.popitem()
        dic_times.popitem()
        dic_traces.popitem()
        dic_x.popitem()
        dic_x_h.popitem()
        # dic_t_id.popitem()
        dic_pos_x.pop(index)
        dic_pos_y.pop(index)
        dic_tot_dist.pop(index)
        print('T < 20 or > 28')
        continue
    
 
    
    #after the checkpoint save values
    dic_y_raw[index] = y1
    
    # =============================================================================
    #     Smoothing (detrended line)
    # =============================================================================
    #first Savitzky-Golay filter
    filter_y = savgol_filter(det_y, 15, 3)
    dic_y_filter[index] = filter_y
    
    
    if figures_on == True:
        fig, ax = plt.subplots()
        ax.plot(timepoints_h,filter_y , 'b-', label='Smoothed function')
        ax.scatter(timepoints_h, det_y, s=2, c='gray', label ='Raw data')
        ax.set(title='Savitzky-Golay filter')
        plt.show()
    
    
    #Smoothing. It uses data after Savitzky-Golay filter. Change smoothed_y to det_y if you wish to not use this prior filter
    from scipy.interpolate import splrep, BSpline
    tck = splrep(timepoints_h, filter_y, s=150)
    tck_s = splrep(timepoints_h, filter_y, s=len(timepoints_h))
    xnew = np.arange(timepoints_h[0], timepoints_h[-1], 0.5) 
       
    
    if figures_on == True:
        # plt.plot(xnew, np.sin(xnew), '-.', label='sin(x)')
        fig, ax = plt.subplots()
        plt.plot(xnew, BSpline(*tck)(xnew), '-', label='s=0')
        plt.plot(xnew, BSpline(*tck_s)(xnew), '-', label=f's={len(timepoints_h)}')
        plt.plot(timepoints_h, det_y, 'o', markersize=2)
        ax.set(title='Smoothing')
        plt.legend()
        plt.show()
    
    
    #Dividing every hour in 60min intervals to obtain smoothed data, increasing the resolution to one point per minute
    x_extended = np.arange(timepoints_h[0], timepoints_h[-1], np.divide(1,60)) 
    y_extended = BSpline(*tck)(x_extended)
    dic_x_smooth[index] = x_extended
    dic_y_smooth[index] = y_extended
    
    
    #finding interpolated values corresponding corresponding to filter_y
    corr_y = y0 + A*np.exp(timepoints_h*lam)*np.cos(f*timepoints_h + phi)
    #Standard error of the regression, a metric that can be used to exclude fittings
    sq_dif_sum = 0
    for l in range(len(filter_y)):
        sq_dif = np.square(filter_y[l]-corr_y[l])
        sq_dif_sum = sq_dif_sum + sq_dif
    SE = np.sqrt(sq_dif_sum/(len(corr_y)-5))  #5 because number of parameters
    SE = abs(SE/np.std(filter_y)) #Correcting the SE dividing by the standard deviation to be able to compare between traces
    dic_SE[index] = SE
    #NEXT use this to find maxima and minima and make a new halfmax value
   
    if figures_on == True:
        fig, ax = plt.subplots()
        ax.scatter(timepoints_h, filter_y, s=2, c='gray', label ='Raw data')
        ax.scatter(timepoints_h, corr_y, s=2, c='blue', label='Fitted function')
        # # ax.plot(x_fitted_lin, y3,'r-', label='First guess')
        ax.set_xlabel('Time (h)')
        ax.set_ylabel('PER2 bioluminiscence', color='black')
        ax.legend(loc="best")
        ax.set(title='Cosinor fit '+ str(index+1))
        ax.text(0,0,'SE = '+ str(SE) ,verticalalignment='center')
        # plt.xticks(np.arange(0, 165, 24)) #possible to adjust for actual time where the recording starts
        # plt.savefig(r'\\vf-mcb-circadian-data\mcb-circadian-data$\Neurofys\Pablo\Membrane excitability and period\EI balance MUA\Script\Ouput graphs\Fitting',dpi=300)
        plt.show()
        
    if SE >= 0.8:           #if fit is too poor
        dic_param.popitem()
        dic_param_all.popitem()
        dic_times.popitem()
        dic_traces.popitem()
        dic_x.popitem()
        dic_x_h.popitem()
        # dic_t_id.popitem()
        dic_pos_x.pop(index)
        dic_pos_y.pop(index)
        dic_y_filter.popitem()
        dic_x_smooth.popitem()
        dic_y_smooth.popitem()
        dic_tot_dist.pop(index)
        dic_SE.popitem()
        print('Bad fit')
        continue        
   
    # =============================================================================
    #     Smoothing (raw line)
    # =============================================================================
    #first part is more like low pass filter
    r_filter_y = savgol_filter(binned_y, 15, 3)
    dic_r_y_filter[index] = r_filter_y
    
    
    if figures_on == True:
        fig, ax = plt.subplots()
        ax.plot(timepoints_h, r_filter_y , 'b-', label='Raw smoothed function')
        ax.scatter(timepoints_h, binned_y, s=2, c='gray', label ='Raw data')
        ax.set(title='Low pass filter (raw)')
        plt.show()
    
    
    #smoothing. Im using the data from the low pass but change smoothed_y to det_y if you wish to not use the filter
    r_tck = splrep(timepoints_h, r_filter_y, s=150)
    r_tck_s = splrep(timepoints_h, r_filter_y, s=len(timepoints_h))
    xnew = np.arange(timepoints_h[0], timepoints_h[-1], 0.5) 
    
    if figures_on == True:
        # plt.plot(xnew, np.sin(xnew), '-.', label='sin(x)')
        fig, ax = plt.subplots()
        plt.plot(xnew, BSpline(*r_tck)(xnew), '-', label='s=0')
        plt.plot(xnew, BSpline(*r_tck_s)(xnew), '-', label=f's={len(timepoints_h)}')
        plt.plot(timepoints_h, binned_y, 'o', markersize=2)
        ax.set(title='Smoothing (raw)')
        plt.legend()
        plt.show()
    
    
    #Dividing every hour in 60h intervals to obtain smoothed data
    r_x_extended = np.arange(timepoints_h[0], timepoints_h[-1], np.divide(1,60)) 
    r_y_extended = BSpline(*r_tck)(x_extended)
    dic_r_x_smooth[index] = r_x_extended
    dic_r_y_smooth[index] = r_y_extended
    
    #NEXT use this to find maxima and minima and make a new halfmax value
    
    # =============================================================================
    #     FINDING Maxima and minima (detrended)
    # =============================================================================
    loc_max = signal.argrelextrema(y_extended, np.greater, order = 400)  #high order better. GIves index of the extreme values
    loc_max = loc_max[0]
    loc_max_y = np.empty(loc_max.size)
    loc_max_y = []
    
    for e in loc_max:
        val = y_extended[e]
        loc_max_y.append(val)
        
    
    loc_min = signal.argrelextrema(y_extended, np.less, order = 400)
    loc_min = loc_min[0]
    loc_min_y = np.empty(loc_min.size)
    loc_min_y = []
    
    for e in loc_min:
        val = y_extended[e]
        loc_min_y.append(val)
        
        
    loc_max_h = []     #local maxima smoothed y values
    loc_min_h = []     # local minimum smoothed y values
    for i in loc_max:
        loc_max_h.append(x_extended[i])
    for i in loc_min:
        loc_min_h.append(x_extended[i])
    

    if len(loc_min_h) <= 1:
        dic_param.popitem()
        dic_param_all.popitem()
        dic_times.popitem()
        dic_traces.popitem()
        dic_x.popitem()
        dic_x_h.popitem()
        # dic_t_id.popitem()
        dic_pos_x.pop(index)
        dic_pos_y.pop(index)
        dic_y_filter.popitem()
        dic_x_smooth.popitem()
        dic_y_smooth.popitem()
        dic_r_y_filter.popitem()
        dic_r_x_smooth.popitem()
        dic_r_y_smooth.popitem()
        dic_tot_dist.pop(index)
        dic_SE.popitem()
        print('only one trough or no trough')
        continue
             
    
    if loc_max_h[0] < loc_min_h[0]:
        loc_max_h = loc_max_h[1:]
        loc_max = loc_max[1:]
        loc_max_y = loc_max_y[1:]

    if loc_max_h[-1] < loc_min_h[-1]:
        loc_min_y = loc_min_y[:-1]
        loc_min = loc_min[:-1]
        loc_min_h = loc_min_h[:-1]
    if len(loc_min_h) == 1:
        dic_param.popitem()
        dic_param_all.popitem()
        dic_times.popitem()
        dic_traces.popitem()
        dic_x.popitem()
        dic_x_h.popitem()
        # dic_t_id.popitem()
        dic_pos_x.pop(index)
        dic_pos_y.pop(index)
        dic_y_filter.popitem()
        dic_x_smooth.popitem()
        dic_y_smooth.popitem()
        dic_r_y_filter.popitem()
        dic_r_x_smooth.popitem()
        dic_r_y_smooth.popitem()
        dic_tot_dist.pop(index)
        dic_SE.popitem()
        print('only one trough2')
        continue
    if len(loc_min_h) != len(loc_max_h):
        dic_param.popitem()
        dic_param_all.popitem()
        dic_times.popitem()
        dic_traces.popitem()
        dic_x.popitem()
        dic_x_h.popitem()
        # dic_t_id.popitem()
        dic_pos_x.pop(index)
        dic_pos_y.pop(index)
        dic_y_filter.popitem()
        dic_x_smooth.popitem()
        dic_y_smooth.popitem()
        dic_r_y_filter.popitem()
        dic_r_x_smooth.popitem()
        dic_r_y_smooth.popitem()
        dic_tot_dist.pop(index)
        dic_SE.popitem()
        print('Peak and trough number does not match')
        continue

    if loc_min_h[0] and loc_min_h[1] < loc_max_h[0]:
        loc_min_y = loc_min_y[1:]
        loc_min = loc_min[1:]
        loc_min_h = loc_min_h[1:]
        
    
    points_max =  list(map(lambda x, y:[x,y], loc_max_h, loc_max_y))   #making x y points with max and min
    points_min =  list(map(lambda x, y:[x,y], loc_min_h, loc_min_y))
    points_all = points_max + points_min
    points_all =  sorted(points_all, key=lambda x: x[0])
    
    
    if points_all[0][1] > 0:  #trimming, the first trought is always detected at the beginning of the trace and last the peak at the end
        del points_all[0]
    if points_all[-1][1] <0:
        del points_all [-1]
    
    if len(points_all)%2 != 0:
        del points_all[-1]
    
    if figures_on == True:
        fig, ax = plt.subplots()
        for e in range(len(points_all)):   
            if points_all[e][1] < 0: 
                ax.scatter(points_all[e][0],points_all[e][1], s=10, c='red')
            if points_all[e][1] > 0: 
                ax.scatter(points_all[e][0],points_all[e][1], s=10,  c='blue')
   
    def pairwise(iterable):
        "s -> (s0, s1), (s2, s3), (s4, s5), ..."
        aa = iter(iterable)
        return zip(aa, aa)
    
    # def line(slo,off):
    halfmax_points_x=[]
        
    for k, l in pairwise(points_all):
        midpoint_y = (l[1]+k[1])/2
        
        # interpolation from the y_extended trace
        def find_nearest(array, value):
            array = np.asarray(array)
            idx = (np.abs(array - value)).argmin()
            return array[idx]
        split = y_extended[x_extended.tolist().index(k[0]):x_extended.tolist().index(l[0])]
        
        midpoint_y_near = find_nearest(split,midpoint_y)
        midpoint_x = x_extended[y_extended.tolist().index(midpoint_y_near)]
        if figures_on == True:    
            ax.scatter(midpoint_x,midpoint_y, s=10, c='green')
       
        halfmax_points_x.append(midpoint_x)

    
    cyc_periods = []
    #extracts periods subtracting each half max value from the previous
    for i in halfmax_points_x:
        j = i+1
        for j in halfmax_points_x:
            if j==i or j<i or halfmax_points_x.index(j)-halfmax_points_x.index(i) > 1:
                continue
            else:
                cyc_per = j-i 
            cyc_periods.append(cyc_per)
    
    
    dic_cyc_per[index] = cyc_periods
    dic_cyc_per_avrg[index] = np.mean(cyc_periods)
    
    
    if figures_on == True:
        ax.plot(x_extended, y_extended,c='gray', label ='Raw data')
        ax.set(title='Relevant points '+str(index+1))
        # plt.savefig(r'\\vf-mcb-circadian-data\mcb-circadian-data$\Neurofys\Pablo\Membrane excitability and period\EI balance MUA\Script\Ouput graphs\Fitting',dpi=300)
        plt.show()
        
        
    # =============================================================================
    #     FINDING Maxima and minima (raw)   everything again but variables have r_ in front ( for raw)
    # =============================================================================
    #No need to calculate the x values again, they should be the same as before

    r_loc_max_y = []
    
    for e in loc_max:
        val = r_y_extended[e]
        r_loc_max_y.append(val)
        
    

    r_loc_min_y = []
    
    for e in loc_min:
        val = r_y_extended[e]
        r_loc_min_y.append(val)
        
     
    
     
    r_points_max =  list(map(lambda x, y:[x,y], loc_max_h, r_loc_max_y))   #making x y points with max and min
    r_points_min =  list(map(lambda x, y:[x,y], loc_min_h, r_loc_min_y))
    r_points_all = r_points_max + r_points_min
    r_points_all =  sorted(r_points_all, key=lambda x: x[0])
    dic_extrema[index] = r_points_all
    

    if figures_on == True:
        fig, ax = plt.subplots()
        for e in range(len(points_all)):   
            if points_all[e][1] < 0: 
                ax.scatter(r_points_all[e][0],r_points_all[e][1], s=10, c='red')
            if points_all[e][1] > 0: 
                ax.scatter(r_points_all[e][0],r_points_all[e][1], s=10,  c='blue')
   
    r_points_x = []   #variable so I can trace back which cycle the amplitude and peaks correspond to
    r_amplitudes = []   
    r_peaks = []
    r_peaks_x = []
    r_troughs = []
    r_single_peaktroughs_x = []
    
    
    for k, l in pairwise(r_points_all):
        r_midpoint_y = (l[1]+k[1])/2
        

        r_amplitude =  l[1]- k[1]
        r_peak = l[1]
        r_trough = k[1]
        r_single_peaktrough_x = (k[0],l[0])
        r_peak_x = l[0]
        # r_halfmax_points_x.append(r_midpoint_x)
        r_amplitudes.append(r_amplitude)
        r_peaks.append(r_peak)
        r_peaks_x.append(r_peak_x)
        r_troughs.append(r_trough)
        r_points_x.append(k[0])
        r_single_peaktroughs_x.append(r_single_peaktrough_x)
        
    dic_amplitude[index] = r_amplitudes
    dic_peaks[index] = r_peaks
    dic_troughs[index] = r_troughs
    dic_r_points_x[index] = r_points_x
    dic_peaks_x[index]= r_peaks_x   #will be used for phase distribution
    dic_single_peaktrough_x[index] = np.array(r_single_peaktroughs_x)
     # for i in half_maxpoints_x:
     
    #testing if big amplitude values lead to a different length of trough-peak or peak-trough (as seen in the mean of all traces)     
    amp_spread = {k: v[0] for k, v in dic_amplitude.items()}
    # print(amp_spread) 
    sorted_amp_spread = sorted(amp_spread, key=amp_spread.get)
    
    amp_mid = len(sorted_amp_spread)//2
    amp_lower_keys = sorted_amp_spread[:amp_mid]
    amp_upper_keys = sorted_amp_spread[amp_mid:]
    
    dic_low_amp = {k: dic_amplitude[k] for k in amp_lower_keys}  #dictionary with the lower amplitude values
    dic_upper_amp = {k: dic_amplitude[k] for k in amp_upper_keys} #dictionary with upper amplitude values (both for the first cycle)
        

     
        
    #not necessary to extract periods from here  
    
    #To exclude based on low amplitude (to select a big region around the scn and only remove those). Only select if amplitude threshold is clear
    
    # if statistics.mean(r_amplitudes) <= 200: # removing spots with very low amplitude (which are selected outside the scn)
    #     dic_param.popitem()
    #     dic_param_all.popitem()
    #     dic_times.popitem()
    #     dic_traces.popitem()
    #     dic_x.popitem()
    #     dic_x_h.popitem()
    #     # dic_t_id.popitem()
    #     dic_pos_x.pop(index)
    #     dic_pos_y.pop(index)
    #     dic_y_filter.popitem()
    #     dic_x_smooth.popitem()
    #     dic_y_smooth.popitem()
    #     dic_r_y_filter.popitem()
    #     dic_r_x_smooth.popitem()
    #     dic_r_y_smooth.popitem()
    #     dic_tot_dist.pop(index)
    #     dic_cyc_per.popitem()
    #     dic_cyc_per_avrg.popitem()
    #     dic_extrema.popitem()
    #     dic_amplitude.popitem()
    #     dic_peaks.popitem()
    #     dic_troughs.popitem()
    #     dic_r_points_x.popitem()
    #     dic_peaks_x.popitem()
    #     print('Low amplitude')
    #     continue
        
    
    
    if figures_on == True:
        ax.plot(r_x_extended, r_y_extended,c='gray', label ='Raw data')
        ax.set(title='Relevant points (raw) '+str(index+1))
        # plt.savefig(r'\\vf-mcb-circadian-data\mcb-circadian-data$\Neurofys\Pablo\Membrane excitability and period\EI balance MUA\Script\Ouput graphs\Fitting',dpi=300)
        plt.show()
     
     

    # #representation of the raw fitted curve  
    
    # # =============================================================================
    # representation of the fitted curve
    # =============================================================================
    
    x_fitted_lin = np.linspace(np.min(timepoints_h), np.max(timepoints_h), 1000)  #extra points for fitted curve
    dic_x_h_lin[index] = x_fitted_lin
    
    # y2 = y0 + A*np.cos((2*np.pi*x_fitted_lin*(1/24)+phi))
    y2 = y0 + A*np.exp(x_fitted_lin*lam)*np.cos(f*x_fitted_lin + phi)
    dic_y_func[index] = y2
    
    y2dif = y0 + A*np.cos(f*x1 + phi)
    if A > 0:        
        y2_norm = y0 + 100*np.cos(f*x_fitted_lin + phi)
    else:
        y2_norm = y0 + -100*np.cos(f*x_fitted_lin + phi)
    dic_y_func_n[index] = y2_norm
    
    #representation of manually approximated values function
    # y3 = cos(225, 140, x1, 0)
    y3= y0_guess + guess_std*np.exp(x_fitted_lin*lam_guess)*np.cos(f_guess*x_fitted_lin+phi_guess)


    if figures_on == True:
        fig, ax = plt.subplots()
        ax.scatter(timepoints_h, det_y, s=2, c='gray', label ='Raw data')
        ax.plot(x_fitted_lin, y2, 'b-', label='Fitted function')
        # # ax.plot(x_fitted_lin, y3,'r-', label='First guess')
        ax.set_xlabel('Time (h)')
        ax.set_ylabel('PER2 bioluminiscence', color='black')
        ax.legend(loc="best")
        ax.set(title='Cosinor fit '+ str(index+1))
        ax.text(0,0,'SE = '+ str(SE) ,verticalalignment='center')
        # plt.xticks(np.arange(0, 165, 24)) #possible to adjust for actual time where the recording starts
        # plt.savefig(r'\\vf-mcb-circadian-data\mcb-circadian-data$\Neurofys\Pablo\Membrane excitability and period\EI balance MUA\Script\Ouput graphs\Fitting',dpi=300)
        plt.show()


    # plt.show()

    res = pd.DataFrame()
    res['output'] = T

    with pd.ExcelWriter(r"\\vf-mcb-circadian-data\mcb-circadian-data$\Neurofys\Pablo\Membrane excitability and period\EI balance MUA\Script\output.xlsx") as writer:        #saving in a new document in different data sheets
        res.to_excel(writer) #corrected 340 (minus bg)
        
        
    # =============================================================================
    #   Analysis of curvature. First derivative of smoothed data
    # =============================================================================
    
    # #increasing resolution
    # interp_func_r = interp1d(dic_r_x_smooth[index], dic_r_y_smooth[index], kind="cubic")
    # interp_func = interp1d(dic_x_smooth[index], dic_y_smooth[index], kind="cubic")
    
    # x_fine_r = np.linspace(dic_r_x_smooth[index].min(), dic_r_x_smooth[index].max(), 5 * len(dic_r_x_smooth[index]))
    # y_fine_r = interp_func_r(x_fine_r)
    
    # x_fine = np.linspace(dic_x_smooth[index].min(), dic_x_smooth[index].max(), 5 * len(dic_x_smooth[index]))
    # y_fine = interp_func(x_fine)
    
    
    dy_dx_r = np.gradient(dic_r_y_smooth[index], dic_r_x_smooth[index])  #derivative of raw
    dy_dx = np.gradient(dic_y_smooth[index], dic_x_smooth[index])   #derivative of detrended
    #smooth first derivatives
    dy_dx_r = savgol_filter(dy_dx_r, window_length=500, polyorder=3)
    dy_dx = savgol_filter(dy_dx, window_length=500, polyorder=3)
    

    # dic_x_fine[index] =x_fine
    dic_dy_dx[index] = dy_dx
  
    
    
    d2y_dx2_r = np.gradient(dy_dx_r, dic_r_x_smooth[index]) #second derivative of raw
    d2y_dx2 = np.gradient (dy_dx, dic_x_smooth[index])       #second derivative fo raw   
    
    dic_d2y_dx[index] = d2y_dx2
    
    #plotting raw
    
    if figures_on == True:
        fig, axs = plt.subplots(3, 1, figsize=(8, 6), sharex=True)
    
        axs[0].plot(dic_r_x_smooth[index], dic_r_y_smooth[index], label="Wave", color="blue")
        axs[0].set_ylabel("y")
        axs[0].legend()
        
        axs[1].plot( dic_r_x_smooth[index], dy_dx_r, label="First Derivative", color="red")
        axs[1].axhline(0, color='black', linestyle='--', alpha=0.5) 
        axs[1].set_ylabel("dy/dx")
        axs[1].legend()
        
        axs[2].plot( dic_r_x_smooth[index], d2y_dx2_r, label="Second Derivative", color="green")
        axs[2].axhline(0, color='black', linestyle='--', alpha=0.5) 
        axs[2].set_ylabel("d²y/dx²")
        axs[2].set_xlabel("x")
        axs[2].legend()
        
        plt.tight_layout()
        plt.show()
        
     #plotting detrended
             
        fig, axs2 = plt.subplots(3, 1, figsize=(8, 6), sharex=True)
        
        axs2[0].plot(dic_x_smooth[index],dic_y_smooth[index], label="Wave", color="blue")
        axs2[0].axhline(0, color='black', linestyle='--', alpha=0.5) 
        axs2[0].set_ylabel("y")
        axs2[0].legend()
        
        axs2[1].plot(dic_x_smooth[index], dy_dx, label="First Derivative", color="red")
        axs2[1].axhline(0, color='black', linestyle='--', alpha=0.5) 
        axs2[1].set_ylabel("dy/dx")
        axs2[1].legend()
        
        axs2[2].plot(dic_x_smooth[index], d2y_dx2, label="Second Derivative", color="green")
        axs2[2].axhline(0, color='black', linestyle='--', alpha=0.5) 
        axs2[2].set_ylabel("d²y/dx²")
        axs2[2].set_xlabel("x")
        axs2[2].legend()
        
        plt.tight_layout()
        plt.show()
    ##############################################
    dic_y[index] = det_y
    dic_fun[index] = y2
    dic_r_y[index] = binned_y
    
    
    dic_index[index] = index
    index = index +1
    

for t in dic_index.values():
    dic_pos_x_up[t] = dic_pos_x[t]   
    dic_pos_y_up[t] = dic_pos_y[t]
   


   

period_cyc_array = np.array(list(dic_cyc_per_avrg.values()), dtype='float32')
period_cyc_array_average = np.mean(period_cyc_array)   #average of averages of cycle to cycle period

period_cos_array = []
for value in dic_param_all.values():
    period_cos = value['T']
    period_cos_array.append(period_cos)
period_cos_array = np.array(period_cos_array, dtype='float32')
period_cos_array_average = np.mean(period_cos_array)    #average cosinor fit T


# =============================================================================
# calculating other parameters (amplitude, trough/peaks, phase)
# =============================================================================


smooth_xy = []  #contains the values of the smoothed raw x and y dictionaries into tuples
raw_values = list(dic_r_y_smooth.values())
for e in list(dic_index.values()):
    raw_traces_xy =[]
    for point in range(len(dic_r_x_smooth[e])):     
        tup = (dic_r_x_smooth[e][point],dic_r_y_smooth[e][point])
        raw_traces_xy.append(tup)
    smooth_xy.append(raw_traces_xy)
    
smallest_first_value = None
largest_last_value = None

for l in smooth_xy:
    current_first_value = l[0][0]
    if smallest_first_value is None or current_first_value < smallest_first_value:
        smallest_first_value = current_first_value
    
    current_last_value =l[-1][0]
    if largest_last_value is None or current_last_value > largest_last_value:
        largest_last_value = current_last_value    

step_size = 1/60
global_x = np.round(np.arange((smallest_first_value), largest_last_value, step_size),2)



sum_y_array = np.zeros_like(global_x, dtype=float)
count_array = np.zeros_like(global_x, dtype=int)
for sublist in smooth_xy:
    for x, y in sublist:
        ind = np.where(global_x == np.round((x),2))[0]
        # Update sum and count arrays
        sum_y_array[ind] += y
        count_array[ind] += 1

# Calculate mean values
mean_list = np.divide(sum_y_array, count_array, out=np.zeros_like(sum_y_array), where=count_array != 0)   
# calculate phase distribution, deviation from the maximum and minimums from the average trace

maxima_indices = signal.argrelextrema(np.array(mean_list), np.greater, order = 300)[0]
minima_indices = signal.argrelextrema(np.array(mean_list), np.less, order = 300)[0]
extrema_indices = np.sort(np.concatenate([maxima_indices, minima_indices]))
x_intervals = [global_x[start:end + 1] for start, end in zip(extrema_indices, extrema_indices[0:])]   

if x_intervals[0] < 20: #might lead to errors when changing time
    x_intervals = x_intervals[1:]     #0 = no skips, use 1 to skip the first extreme value. Make sure the first detected extreme value is always a trough!!
    # if x_intervals[0] <20:
    #     x_intervals = x_intervals[1:]
    
x_intervals_troughs = []
x_intervals_peaks = []
for e in range(len(x_intervals)):   #this loop is shortening the intervals to omit those that include peaks. So only troughs are detected
    if e == 0 or e%2 ==0:
        x_intervals_troughs.append(x_intervals[e][0])
    else:
        x_intervals_peaks.append(x_intervals[e][0])

#now calculating the distance between successive throughs and peaks
def subtract_even_odd(lst):
    # Subtract even - odd values
    trough_peak_diff = [lst[i+1] - lst[i] for i in range(0, len(lst)-1, 2)]
    
    # Subtract odd - even values (excluding first element)
    peak_trough_diff = [lst[i] - lst[i-1] for i in range(2, len(lst), 2)]
    
    return trough_peak_diff, peak_trough_diff

trough_peak_diff, peak_trough_diff = subtract_even_odd(x_intervals)
# Make both lists the same length by padding with NaN
max_len = max(len(trough_peak_diff), len(peak_trough_diff))

trough_peak_diff = [x.item() for x in trough_peak_diff]
peak_trough_diff = [x.item() for x in peak_trough_diff]



trough_peak_diff += [np.nan] * (max_len - len(trough_peak_diff))
peak_trough_diff += [np.nan] * (max_len - len(peak_trough_diff))


# x_intervals_troughs = np.asarray(x_intervals_troughs)
# half_max_values = [(mean_list[start] + mean_list[end]) / 2 for start, end in zip(extrema_indices, extrema_indices[1:])]
# # crossings = np.where(np.diff(np.sign(np.array(mean_list) - half_max_values[0])))[0]
# # adjusted_crossings = np.clip(crossings, 0, len(mean_list) - 1)

# # x_intervals = [global_x[start:end + 1] for start, end in zip(adjusted_crossings, adjusted_crossings[1:])]
# # y_intervals = [mean_list[start:end + 1] for start, end in zip(adjusted_crossings, adjusted_crossings[1:])]
# crossings = [np.where(np.diff(np.sign(np.array(mean_list[start:end + 1]) - half_max_values))[0])[0][0] + start
#               for half, start, end in zip(half_max_values, extrema_indices[:-1], extrema_indices[1:])]

# # Adjust the crossings indices to ensure they don't go beyond the length of the data
# adjusted_crossings = np.clip(crossings, 0, len(mean_list) - 1)

# # Split the x-values and mean_list into intervals based on crossings
# x_intervals = [global_x[start:end + 1] for start, end in zip(adjusted_crossings[:-1], adjusted_crossings[1:])]

#making dictionary with the correct values of cycles that can be matched 
for q in dic_index.values():
    lis = dic_r_points_x[q]
    cyc_index_list = []
    for e in lis:
        near =find_nearest(x_intervals_troughs, e)
        cyc_index = x_intervals_troughs.index(near)
        cyc_index_list.append(cyc_index)
    dic_cyc_index[q]= cyc_index_list


#calculating mean amplitude per cycle

# find all unique values in the values of dic_index
unique_values = set(value for values in dic_cyc_index.values() for value in values)

# calculate mean for each unique value
dic_amplitude_cyc = {}
dic_amplitude_list = {}
for unique_value in unique_values:
    # Collect all positions of unique_value in dic_index
    positions = [(key, pos) for key, values in dic_cyc_index.items() for pos, value in enumerate(values) if value == unique_value]

    # Collect amplitude values based on positions
    amplitude_values = [dic_amplitude[key][pos] for key, pos in positions]  #can be used to extract ALL means and standar deviation
    
    # Calculate the mean for the collected values
    mean_value = np.mean(amplitude_values)
    dic_amplitude_list[unique_value] = amplitude_values
    dic_amplitude_cyc[unique_value] = mean_value

#calculating mean peaks per cycle
dic_peaks_cyc = {}
dic_peaks_list = {}
for unique_value in unique_values:
    positions = [(key, pos) for key, values in dic_cyc_index.items() for pos, value in enumerate(values) if value == unique_value]
    peak_values = [dic_peaks[key][pos] for key, pos in positions]  #can be used to extract ALL means and standar deviation
    mean_value = np.mean(peak_values)
    dic_peaks_list[unique_value] = peak_values
    dic_peaks_cyc[unique_value] = mean_value

#calculating mean troughs per cycle
dic_troughs_cyc = {}
dic_troughs_list = {}
for unique_value in unique_values:
    positions = [(key, pos) for key, values in dic_cyc_index.items() for pos, value in enumerate(values) if value == unique_value]
    trough_values = [dic_troughs[key][pos] for key, pos in positions]  #can be used to extract ALL means and standar deviation
    mean_value = np.mean(trough_values)
    dic_troughs_list[unique_value] = trough_values
    dic_troughs_cyc[unique_value] = mean_value
    
#Calculating peak trough distance
dic_single_peaktrough_x_cyc = {}
dic_single_peaktrough_x_list = {}
# dic_single_peaktrough_values = np.array()
for unique_value in unique_values:
    positions = [(key, pos) for key, values in dic_cyc_index.items() for pos, value in enumerate(values) if value == unique_value]
    single_peaktrough_values = [dic_single_peaktrough_x[key][pos] for key, pos in positions]  #can be used to extract ALL means and standar deviation
    print(single_peaktrough_values)
    mean_value = np.mean(single_peaktrough_values, axis =0)
    dic_single_peaktrough_x_list[unique_value] = single_peaktrough_values
    dic_single_peaktrough_x_cyc[unique_value] = mean_value

#calculating peak/trough distance for those cells with higher/lower amplitude values
# dic_low_amp_x_cyc ={}
# dic_low_amp_x_list = {}
# for unique_value in unique_values:
#     positions = [(key, pos) for key, values in dic_cyc_index.items() for pos, value in enumerate(values) if value == unique_value]
#     single_peaktrough_low_values = [dic_low_amp[key][pos] for key, pos in positions]  #can be used to extract ALL means and standar deviation
#     print(single_peaktrough_values)
#     mean_value = np.mean(single_peaktrough_low_values, axis =0)
#     dic_low_amp_x_list[unique_value] = single_peaktrough_values
#     dic_low_amp_x_cyc[unique_value] = mean_value




#phase distribution
dic_peaks_x_cyc = {}
dic_peaks_x_list = {}
dic_phase_dif = {}
for unique_value in unique_values:
    positions = [(key, pos) for key, values in dic_cyc_index.items() for pos, value in enumerate(values) if value == unique_value]
    peak_x_values = [dic_peaks_x[key][pos] for key, pos in positions]  #can be used to extract ALL means and standar deviation
    mean_value = np.mean(peak_x_values)
    dic_peaks_x_cyc[unique_value] = mean_value
    phase_dif = []
    for g in peak_x_values:
        phase_dif.append(g-mean_value)
    dic_phase_dif[unique_value] = phase_dif
    dic_peaks_x_list[unique_value] = peak_x_values

#Period (structured the same way as the other variables) 
dic_cyc_per_or = dic_cyc_per
dic_cyc_per_or_cyc = {}
dic_cyc_per_list = {}
for x in dic_cyc_per_or.values():       #First adds 0s when the size doesnt match the other dictionaries. 
    x.append(0)
for h in dic_amplitude.keys():
    if len(dic_cyc_per_or[h]) < len(dic_amplitude[h]):
        dic_cyc_per_or[h].append(0)

for unique_value in unique_values:
    positions = [(key, pos) for key, values in dic_cyc_index.items() for pos, value in enumerate(values) if value == unique_value]
    period_values = [dic_cyc_per_or[key][pos] for key, pos in positions]  #can be used to extract ALL means and standar deviation
    period_values_without_zeros = np.where(np.array(period_values) == 0, np.nan, period_values)
    mean_value = np.nanmean(period_values_without_zeros)  #0s to Nan so the mean is only of the values that are not 0 (0s arent real data)
    dic_cyc_per_list[unique_value] = period_values
    dic_cyc_per_or_cyc[unique_value] = mean_value 
    


#other way to display the data: padding with NaN so each cycle can be visualised with NaN
dic_peaklist_percycle = {}
dic_amplitude_percycle = {}
dic_period_percycle = {}
for unique_value in unique_values:
    # Initialize the list with NaNs
    # values_for_unique_value_peak = [np.nan] * len(dic_cyc_index)
    # values_for_unique_value_amp = [np.nan] * len(dic_cyc_index)
    # values_for_unique_value_per = [np.nan] * len(dic_cyc_index)
    values_for_unique_value_peak_dic = {key: np.nan for key in dic_index.values()}
    values_for_unique_value_amp_dic = {key: np.nan for key in dic_index.values()}
    values_for_unique_value_per_dic = {key: np.nan for key in dic_index.values()}
    # Update specific positions with values from dic_peaks_x
    for key, values in dic_cyc_index.items():
        if unique_value in values:
            pos = values.index(unique_value)
            # Ensure the index is within bounds and use len(dic_peaks_x[key]) for indexing
            # if 0 <= int(key) < len(dic_peaks_x) and 0 <= pos < len(dic_peaks_x[key]):
                #Update no need to ensure because structure changed to NaN dictionaries instead of lists. Fixes indexing
            values_for_unique_value_peak_dic[int(key)] = np.float64(dic_peaks_x[key][pos])
            values_for_unique_value_amp_dic[int(key)] = np.float64(dic_amplitude[key][pos])
            values_for_unique_value_per_dic[int(key)] = np.float64(dic_cyc_per_or[key][pos])
            
    dic_peaklist_percycle[str(unique_value)] = values_for_unique_value_peak_dic
    dic_amplitude_percycle[str(unique_value)] = values_for_unique_value_amp_dic
    dic_period_percycle[str(unique_value)] = values_for_unique_value_per_dic

# new_dict now contains the desired mapping with NaNs and values from dic_peaks_x at the correct positions


#Now finish the peaktrough analysis but for single cell averages
x_intervals_single = [value for sublist in dic_single_peaktrough_x_cyc.values() for value in sublist] #bad name but ok
#now calculating the distance between successive throughs and peaks
def subtract_even_odd(lst):
    # Subtract even - odd values
    single_trough_peak_diff = [lst[i+1] - lst[i] for i in range(0, len(lst)-1, 2)]
    
    # Subtract odd - even values (excluding first element)
    single_peak_trough_diff = [lst[i] - lst[i-1] for i in range(2, len(lst), 2)]
    
    return single_trough_peak_diff, single_peak_trough_diff

single_trough_peak_diff, single_peak_trough_diff = subtract_even_odd(x_intervals_single)
# Make both lists the same length by padding with NaN
max_len = max(len(single_trough_peak_diff), len(single_peak_trough_diff))

single_trough_peak_diff = [x.item() for x in single_trough_peak_diff]
single_peak_trough_diff = [x.item() for x in single_peak_trough_diff]



single_trough_peak_diff += [np.nan] * (max_len - len(single_trough_peak_diff))
single_peak_trough_diff += [np.nan] * (max_len - len(single_peak_trough_diff))



# =============================================================================
# FIGURES
# =============================================================================
def create_folder(folder_name, location):
    # Join the location and folder name to create the full path
    folder_path = os.path.join(location, folder_name)
    
    # Check if the folder exists
    if not os.path.exists(folder_path):
        # If not, create the folder
        os.makedirs(folder_path)
        print(f"Folder '{folder_path}' created successfully.")
    else:
        print(f"Folder '{folder_path}' already exists.")

folder_name = exp
folder_location = path+ "\\Script output\\"
create_folder(folder_name, folder_location)       
       

create_folder('Figures',path+ "\\Script output\\" +exp+ "\\" )    
fig_folder = path+ "\\Script output\\" +exp+ "\\"+"Figures"
print('figures')
list_index = dic_index.values()



w, h = figaspect(2/7)
fig, ax = plt.subplots(figsize=(w,h))
for t in list(dic_index.values()):   
    ax.scatter(dic_x[t],dic_y[t], s=2)
ax.set_xlabel('Time (days)')
ax.set_ylabel('PER2 bioluminescence', color='black')
ax.set(title='Detrended (date)')
ax.legend(loc="best")
plt.savefig(fig_folder + '\Raw_detrended_all',dpi=300)
plt.show()



fig, ax = plt.subplots(figsize=(w,h))
# ax.scatter(x1, det_y, s=2, c='gray', label ='Raw data')
for t in list(dic_index.values()):   
    ax.plot(dic_x[t],dic_y_raw[t],linewidth=0.5)
ax.set_xlabel('Time (days)')
ax.set_ylabel('PER2 bioluminescence', color='black')
ax.set(title='Raw data')
ax.legend(loc="best")
plt.savefig(fig_folder + '\Raw data all',dpi=300)
plt.show()

fig, ax = plt.subplots(figsize=(w,h))
# ax.scatter(x1, det_y, s=2, c='gray', label ='Raw data')
for t in list(dic_index.values()):   
    ax.plot(dic_x_h[t],dic_y_raw[t],linewidth=0.5)
ax.set_xlabel('Time (days)')
ax.set_ylabel('PER2 bioluminescence', color='black')
ax.set(title='Raw data h')
ax.legend(loc="best")
# Set x-axis ticks to show every 24 points
ax.xaxis.set_major_locator(MultipleLocator(24))
# Remove the upper and right spines (borders)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.savefig(fig_folder + '\Raw data h',dpi=300)
plt.show()

fig, ax = plt.subplots(figsize=(w,h))
# ax.scatter(x1, det_y, s=2, c='gray', label ='Raw data')
for t in list(dic_index.values()):   
    ax.plot(dic_x_smooth[t],dic_dy_dx[t],linewidth=0.5)
ax.set_xlabel('Time (days)')
ax.set_ylabel('PER2 bioluminescence', color='black')
ax.set(title='dy/dx')
ax.legend(loc="best")
# Set x-axis ticks to show every 24 points
ax.xaxis.set_major_locator(MultipleLocator(24))
# Remove the upper and right spines (borders)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.savefig(fig_folder + '\Raw data h',dpi=300)
plt.show()

fig, ax = plt.subplots(figsize=(w,h))
# ax.scatter(x1, det_y, s=2, c='gray', label ='Raw data')
for t in list(dic_index.values()):   
    ax.plot(dic_x_smooth[t],dic_d2y_dx[t],linewidth=0.5)
ax.set_xlabel('Time (days)')
ax.set_ylabel('PER2 bioluminescence', color='black')
ax.set(title='d2y/dx')
ax.legend(loc="best")
# Set x-axis ticks to show every 24 points
ax.xaxis.set_major_locator(MultipleLocator(24))
# Remove the upper and right spines (borders)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.savefig(fig_folder + '\Raw data h',dpi=300)
plt.show()

fig, ax = plt.subplots(figsize=(w,h))
# ax.scatter(x1, det_y, s=2, c='gray', label ='Raw data')
for t in list(dic_index.values()):   
    ax.plot(dic_r_x_smooth[t],dic_r_y_smooth[t],linewidth=0.5)
ax.plot(global_x,mean_list, linewidth =2, color = 'black')
for element in x_intervals:
    ax.vlines(element[0],min(mean_list),max(mean_list))
ax.set_xlabel('Time (h)')
ax.set_ylabel('PER2 bioluminescence', color='black')
ax.set(title='Raw smoothed data')
ax.legend(loc="best")
plt.show()



afont = {'fontname':'Arial'}



w, h = figaspect(2/7)
fig, ax = plt.subplots(figsize=(w,h))
# ax.scatter(x1, det_y, s=2, c='gray', label ='Raw data')
for t in list(dic_index.values()):   
    x_data = dic_r_x_smooth[t]  # Get the corresponding x data
    y_data = dic_r_y_smooth[t]  # Get the corresponding y data   
    ax.plot(x_data[750:],y_data[750:],linewidth=0.75)
# ax.plot(global_x[750:],mean_list[750:], linewidth =2, color = 'black')
# ax.plot(global_x,mean_list, linewidth =2, color = 'black')
ax.set_xlabel('Time (h)', **afont)
ax.set_ylabel('PER2 bioluminescence (AU)', color='black', **afont)
ax.set(title='Smoothed traces')

# Set x-axis ticks to show every 24 points
ax.xaxis.set_major_locator(MultipleLocator(24))
# Remove the upper and right spines (borders)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.savefig(fig_folder + '\Smoothed data manuscript 075',dpi=600)
plt.savefig(fig_folder + '\Smoothed data manuscript 075.eps')
plt.show()



fig, ax = plt.subplots(figsize=(w,h))
# ax.scatter(x1, det_y, s=2, c='gray', label ='Raw data')

for t in list(dic_index.values()):   
    x_data = dic_x_h[t]
    y_data = dic_y_raw[t]
    ax.plot(x_data[12:],y_data[12:],linewidth=0.75)
# ax.plot(global_x,mean_list, linewidth =2, color = 'black')
ax.set_xlabel('Time (h)', **afont)
ax.set_ylabel('PER2 bioluminescence (AU)', color='black', **afont)

# Set x-axis ticks to show every 24 points
ax.xaxis.set_major_locator(MultipleLocator(24))
# Remove the upper and right spines (borders)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.savefig(fig_folder + '\Raw data manuscript 075',dpi=1200)
plt.savefig(fig_folder + '\Raw data manuscript 075.eps')
plt.show()



fig, ax = plt.subplots(figsize=(w,h))
# ax.scatter(x1, det_y, s=2, c='gray', label ='Raw data')
for t in list(dic_index.values()):   
    ax.plot(dic_x_smooth[t],dic_y_smooth[t],linewidth=0.5)
ax.set_xlabel('Time (h)')
ax.set_ylabel('PER2 bioluminescence', color='black')
ax.set(title='Smoothed detrended data (date)')
ax.xaxis.set_major_locator(MultipleLocator(3))
ax.legend(loc="best")
plt.savefig(fig_folder + '\Smoothed detrended all',dpi=300)
plt.show()


fig, ax = plt.subplots(figsize=(w,h))
for t in list(dic_index.values()):   
    ax.scatter(dic_x_h[t],dic_y[t], s=2)
ax.set_xlabel('Time (h)')
ax.set_ylabel('PER2 bioluminescence', color='black')
ax.set(title='Detrended data(hours)')
plt.savefig(fig_folder + '\Detrended (hours)',dpi=300)
ax.legend(loc="best")
plt.show()


fig, ax = plt.subplots(figsize=(w,h))
for t in list(dic_index.values()):   
    ax.plot(dic_x_h_lin[t],dic_y_func[t],linewidth=0.5)
ax.set_xlabel('Time (h)')
ax.set_ylabel('PER2 bioluminescence', color='black')
ax.set(title='Cosinor fit detrended data')
ax.legend(loc="best")
plt.savefig(fig_folder + '\Cosinor fit_all',dpi=300)
plt.show()


fig, ax = plt.subplots(figsize=(w,h))
for t in list(dic_index.values()):    
    ax.plot(dic_x_h_lin[t],dic_y_func_n[t],linewidth=0.5)
ax.set_xlabel('Time (h)')
ax.set_ylabel('PER2 bioluminescence', color='black')
ax.set(title='Normalised data')
ax.legend(loc="best")
plt.savefig(fig_folder + '\\Normalised_all',dpi=300)
plt.show()

#number of detected cells
n_det_cells = len(dic_index)


#colour gradient for period

fig, ax = plt.subplots()
# els = list(dic_x_h.items())
els_list = list(dic_x_h.values())
els_list2 =[]
for v in els_list:
    els_list2.append(v[0])
lat_first = max(els_list2)    #from the trace that gets detected the latest, the first point
# lat_first = els[-1][1][0] #from the trace that gets detected the latest, the first point

vmin, vmax = 20,28

cmap = 'rainbow'

colors = []
for t in list(dic_index.values()):
    lat_first_index = np.where(dic_x_h[t] == lat_first)[0]  #gives index from the lat_first for this trace. SO every trace has its position assessed at the same timepoint in the recording
    
    ax.scatter(dic_pos_x_up[t][lat_first_index], dic_pos_y_up[t][lat_first_index], s=5, c=dic_cyc_per_avrg[t], cmap=cmap, vmin = vmin, vmax = vmax)
    colors.append(dic_cyc_per_avrg[t])
    # ax.annotate(t, (dic_pos_x_up[t][lat_first_index], dic_pos_y_up[t][lat_first_index]), fontsize= 5)        #if you want the numbers  (from this analysis)
    # ax.annotate(dic_t_id[t], (dic_pos_x_up[t][lat_first_index], dic_pos_y_up[t][lat_first_index]), fontsize= 2)     #if you want the numbers (from trackmate)

colors = np.array(colors)

norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # empty array for the scalar mappable
cbar = plt.colorbar(sm, ax=ax, label='Period (h)')


ax.set_aspect('equal', adjustable='box')    
plt.xticks(visible=False)
plt.yticks(visible=False)

plt.savefig(fig_folder + '\Map_SCN_period',dpi=600)
plt.show()



fig, ax = plt.subplots()
# els = list(dic_x_h.items())
els_list = list(dic_x_h.values())
els_list2 =[]
for v in els_list:
    els_list2.append(v[0])
lat_first = max(els_list2)    #from the trace that gets detected the latest, the first point
# lat_first = els[-1][1][0] #from the trace that gets detected the latest, the first point

vmin, vmax = 22,26

cmap = 'rainbow'
#rotation transformation
def rotate_points(x, y, angle_deg):
    x = np.array(x)
    y= np.array(y)
    angle_rad = np.deg2rad(angle_deg)
    x_rot = x * np.cos(angle_rad) - y * np.sin(angle_rad)
    y_rot = x * np.sin(angle_rad) + y * np.cos(angle_rad)
    return x_rot, y_rot

rot_angle =-122


colors = []
for t in list(dic_index.values()):
    lat_first_index = np.where(dic_x_h[t] == lat_first)[0]   #gives index from the lat_first for this trace. SO every trace has its position assessed at the same timepoint in the recording
    x_rot, y_rot = rotate_points(dic_pos_x_up[t][lat_first_index], dic_pos_y_up[t][lat_first_index], angle_deg=rot_angle)
    cl=dic_cyc_per_avrg[t]
    if cl < 22.5 :
        continue
    else:
        ax.scatter(x_rot, y_rot, s=50, alpha=0.6, c=dic_cyc_per_avrg[t], cmap=cmap, vmin = vmin, vmax = vmax,  edgecolors='none')
        colors.append(dic_cyc_per_avrg[t])
        # ax.annotate(t, (dic_pos_x_up[t][lat_first_index], dic_pos_y_up[t][lat_first_index]), fontsize= 2)        #if you want the numbers  (from this analysis)
        # ax.annotate(dic_t_id[t], (dic_pos_x_up[t][lat_first_index], dic_pos_y_up[t][lat_first_index]), fontsize= 2)     #if you want the numbers (from trackmate)

colors = np.array(colors)

norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
sm = cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])  # empty array for the scalar mappable
cbar = plt.colorbar(sm, ax=ax, label='Period (h)')

for spine in ax.spines.values():        #hiding axes and ticks
    spine.set_visible(False)
ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

ax.set_aspect('equal', adjustable='box')    
plt.xticks(visible=False)
plt.yticks(visible=False)

plt.savefig(fig_folder + '\Map_SCN_period_22-25',dpi=600)
plt.savefig(fig_folder + '\Map_SCN_period_22-25.eps' )
plt.show()


# =============================================================================
# OUTPUT
# =============================================================================

output = pd.DataFrame()
df_raw_h = pd.DataFrame() 
df_raw = pd.DataFrame()
df_amplitude = pd.DataFrame()
df_troughs =pd.DataFrame()
df_peaks = pd.DataFrame()
df_phase = pd.DataFrame()
df_per = pd.DataFrame()
df_cos_fit = pd.DataFrame()
df_all_cycle = pd.DataFrame()
df_meantrace = pd.DataFrame(mean_list)
df_peak_trough_diff = pd.DataFramedf = pd.DataFrame({
    "Trough_peak_diff": trough_peak_diff,
    "Peak_through_diff": peak_trough_diff
})
df_peak_trough_diff_single = pd.DataFramedf = pd.DataFrame({
    "Trough_peak_diff_single": single_trough_peak_diff,
    "Peak_through_diff_single": single_peak_trough_diff
})
# unique_x = set(value for values in dic_r_x_smooth.values() for value in values)


#Raw traces
for e in dic_index.values():
    df_raw_h[str(e)+'_x'] = pd.Series(dic_x_h[e])
    df_raw_h[str(e)+'_y'] = pd.Series(dic_y_raw[e])

# Raw traces smoothed
for e in dic_index.values():
    # df_raw.insert(-1,str(e), dic_r_x_smooth[e])
    df_raw[str(e)+'_x'] = pd.Series(dic_r_x_smooth[e])
    df_raw[str(e)+'_y'] = pd.Series(dic_r_y_smooth[e])


for v in unique_values:
    #Amplitude
    df_amplitude = pd.DataFrame.from_dict(dic_amplitude_list, orient = 'index').transpose()
    
    #troughs
    df_troughs = pd.DataFrame.from_dict(dic_troughs_list, orient = 'index').transpose()

    #peaks
    df_peaks = pd.DataFrame.from_dict(dic_peaks_list, orient = 'index').transpose()

    #phase dif
    df_phase = pd.DataFrame.from_dict(dic_peaks_x_list, orient = 'index').transpose()

    #phase dif
    df_phase_dif = pd.DataFrame.from_dict(dic_phase_dif, orient = 'index').transpose()

    #period
    df_per = pd.DataFrame.from_dict(dic_cyc_per_list, orient = 'index').transpose()

# #Cosinor fit for every trace


#Now one for averages per cycle
# for x, y in dic_param:
#     df_all_cycle.
    
# for keys in dic_amplitude_cyc.keys():
df_all_amp = pd.DataFrame(dic_amplitude_cyc.values()).transpose()
df_all_troughs = pd.DataFrame(dic_troughs_cyc.values()).transpose()
df_all_peaks = pd.DataFrame(dic_peaks_cyc.values()).transpose()
df_all_phase = pd.DataFrame(dic_peaks_x_cyc.values()).transpose()
df_all_per = pd.DataFrame(dic_cyc_per_or_cyc.values()).transpose()

df_all_cycle = pd.concat([df_all_amp, df_all_peaks, df_all_troughs, df_all_phase, df_all_per])
df_all_cycle.index = ['Amplitude', 'Peaks', 'Troughs', 'Phase', 'Period']
    
with pd.ExcelWriter(path+ "\\Script output\\" +exp+ "\\"+ 'Output_'+exp +'.xlsx') as writer:        #saving in a new document in different data sheets
    df_raw_h.to_excel(writer, sheet_name='Raw') 
    df_raw.to_excel(writer, sheet_name='Raw_smoothed') 
    df_meantrace.to_excel(writer, sheet_name='Mean_trace')
    df_amplitude.to_excel(writer, sheet_name='Amplitude')
    df_peaks.to_excel(writer, sheet_name='Peaks')
    df_troughs.to_excel(writer, sheet_name='Troughs')
    df_phase.to_excel(writer, sheet_name='Phase')
    df_phase_dif.to_excel(writer, sheet_name='Phase dif')
    df_per.to_excel(writer, sheet_name='Period')
    df_peak_trough_diff.to_excel(writer, sheet_name='peak_trough_dif')
    df_peak_trough_diff_single.to_excel(writer, sheet_name='peak_trough_dif_single')
    df_all_cycle.to_excel(writer, sheet_name='Averages per cycle')


# =============================================================================
# Additional figures
# =============================================================================
all_phases =  [item for sublist in dic_peaks_x_list.values() for item in sublist]

fig, ax = plt.subplots(figsize=(w,h))
ax.hist(all_phases, bins=150, facecolor='black')
plt.show()



# colors = []

for value in range(len(unique_values)): 
    fig, ax2= plt.subplots()
    i=0    
    cycle = str(value)      #choose which cycle to look at. if the cell has no data for that cycle it shows in grey
    col_list = list(dic_peaklist_percycle[cycle].values())
    vmin, vmax = np.nanmean(col_list)-3, np.nanmean(col_list)+3
    
    cmap = 'viridis'
    for t in list(dic_index.values()):
        lat_first_index = np.where(dic_x_h[t] == lat_first)[0]   #gives index from the lat_first for this trace. SO every trace has its position assessed at the same timepoint in the recording
        
        
        if np.isnan(col_list[i]):
            ax2.scatter(dic_pos_x_up[t][lat_first_index], dic_pos_y_up[t][lat_first_index], s=40, alpha=0.6, c='grey')
        else:
            ax2.scatter(dic_pos_x_up[t][lat_first_index], dic_pos_y_up[t][lat_first_index], s=40, alpha=0.6, c=col_list[i], cmap=cmap, vmin = vmin, vmax = vmax)
        # colors.append(col_list[i])
        # ax.annotate(t, (dic_pos_x_up[t][lat_first_index], dic_pos_y_up[t][lat_first_index]), fontsize= 5)        #if you want the numbers  (from this analysis)
        # ax.annotate(dic_t_id[t], (dic_pos_x_up[t][lat_first_index], dic_pos_y_up[t][lat_first_index]), fontsize= 2)     #if you want the numbers (from trackmate)
        i = i+1
    # colors = np.array(colors)
    
    ax2.set_aspect('equal', adjustable='box')    
    plt.xticks(visible=False)
    plt.yticks(visible=False)
    
    
    norm2 = mcolors.Normalize(vmin=vmin, vmax=vmax)
    sm2 = cm.ScalarMappable(cmap=cmap, norm=norm2)
    sm2.set_array([])  # empty array for the scalar mappable
    cbar2 = plt.colorbar(sm2,  ax=ax, label='Phase (h)')
    
    plt.savefig(fig_folder + '\Map_SCN_phase_cycle_'+str(value+1),dpi=600)
    plt.show()

#map amplitudes

for value in range(len(unique_values)): 
    fig, ax2= plt.subplots()
    i=0    
    cycle = str(value)      #choose which cycle to look at. if the cell has no data for that cycle it shows in grey
    col_list = list(dic_amplitude_percycle[cycle].values())
    vmin, vmax = 0, np.nanmax(col_list)  
    if vmax < vmin:
        vmax = 0
    cmap = 'plasma'
    for t in list(dic_index.values()):
        lat_first_index = np.where(dic_x_h[t] == lat_first)[0]   #gives index from the lat_first for this trace. SO every trace has its position assessed at the same timepoint in the recording
        
        
        if np.isnan(col_list[i]):
            ax2.scatter(dic_pos_x_up[t][lat_first_index], dic_pos_y_up[t][lat_first_index], s=40, alpha=0.6, c='grey')
        else:
            ax2.scatter(dic_pos_x_up[t][lat_first_index], dic_pos_y_up[t][lat_first_index], s=40, alpha=0.6, c=col_list[i], cmap=cmap, vmin = vmin, vmax = vmax)
            # ax2.scatter(dic_pos_x_up[t][lat_first_index], dic_pos_y_up[t][lat_first_index], s=40, alpha=0.6, c='green', cmap=cmap, vmin = vmin, vmax = vmax)
            # ax2.annotate(, (dic_pos_x_up[t][lat_first_index], dic_pos_y_up[t][lat_first_index]), fontsize= 5)        #if you want the numbers  (from this analysis)
        # colors.append(col_list[i])
        # ax2.annotate(t, (dic_pos_x_up[t][lat_first_index], dic_pos_y_up[t][lat_first_index]), fontsize= 5)        #if you want the numbers  (from this analysis)
        # ax.annotate(dic_t_id[t], (dic_pos_x_up[t][lat_first_index], dic_pos_y_up[t][lat_first_index]), fontsize= 2)     #if you want the numbers (from trackmate)
        i = i+1
    # colors = np.array(colors)
    
    ax2.set_aspect('equal', adjustable='box')    
    plt.xticks(visible=False)
    plt.yticks(visible=False)
    
    
    norm2 = mcolors.Normalize(vmin=vmin, vmax=vmax)
    sm2 = cm.ScalarMappable(cmap=cmap, norm=norm2)
    sm2.set_array([])  # empty array for the scalar mappable
    cbar2 = plt.colorbar(sm2, ax=ax, label='Amplitude')
    
    plt.savefig(fig_folder + '\Map_SCN_amplitude_cycle_'+str(value+1),dpi=600)
    plt.show()


#map periods

for value in range(len(unique_values)): 
    fig, ax2= plt.subplots()
    i=0    
    cycle = str(value)      #choose which cycle to look at. if the cell has no data for that cycle it shows in grey
    col_list = list(dic_period_percycle[cycle].values())
    vmin, vmax = 22,26
    
    cmap = 'rainbow'
    for t in list(dic_index.values()):
        lat_first_index = np.where(dic_x_h[t] == lat_first)[0]   #gives index from the lat_first for this trace. SO every trace has its position assessed at the same timepoint in the recording
        
        
        if np.isnan(col_list[i]) or col_list[i] == 0:  #removing values of 0 to visualise better
            ax2.scatter(dic_pos_x_up[t][lat_first_index], dic_pos_y_up[t][lat_first_index], s=40, alpha=0.6, c='grey')
        else:
            ax2.scatter(dic_pos_x_up[t][lat_first_index], dic_pos_y_up[t][lat_first_index], s=40, alpha=0.6, c=col_list[i], cmap=cmap, vmin = vmin, vmax = vmax)
        # colors.append(col_list[i])
        # ax.annotate(t, (dic_pos_x_up[t][lat_first_index], dic_pos_y_up[t][lat_first_index]), fontsize= 5)        #if you want the numbers  (from this analysis)
        # ax.annotate(dic_t_id[t], (dic_pos_x_up[t][lat_first_index], dic_pos_y_up[t][lat_first_index]), fontsize= 2)     #if you want the numbers (from trackmate)
        i = i+1
    # colors = np.array(colors)
    
    ax2.set_aspect('equal', adjustable='box')    
    plt.xticks(visible=False)
    plt.yticks(visible=False)
    
    
    norm2 = mcolors.Normalize(vmin=vmin, vmax=vmax)
    sm2 = cm.ScalarMappable(cmap=cmap, norm=norm2)
    sm2.set_array([])  # empty array for the scalar mappable
    cbar2 = plt.colorbar(sm2, ax=ax, label='Period')
    
    plt.savefig(fig_folder + '\Map_SCN_period_cycle_'+str(value+1),dpi=600)
    plt.show()
