# -*- coding: utf-8 -*-
"""

@author: jessy
"""
import numpy as np 
import pandas as pd

import matplotlib.pyplot as plt

"""Plot given accelerometer data"""
def plot_acc(acc):
    plt.figure(figsize=(10, 7))
    plt.title("Acc data")
    plt.xlabel("Time (s)")
    plt.ylabel("acceleration (g)")
    plt.plot(np.arange(len(acc)), acc, label="acc (g)", color="#A651D8")

    plt.show()

"""allow 3 secs calibration time to determine cane position when placed upright, 
which marks the beginning of taking a step forward for the user""" 
def calibrate_cane_upright(timestamp, acc_x, acc_y, acc_z):
    ini_x = []
    ini_y = []
    ini_z = []
    for i in range(0,len(timestamp)): 
        if timestamp[i] <= 3: 
            ini_x.append(acc_x[i])
            ini_y.append(acc_y[i])
            ini_z.append(acc_z[i])
    mean_x = np.mean(ini_x)
    mean_y = np.mean(ini_y)
    mean_z = np.mean(ini_z)
    
    return mean_x, mean_y, mean_z

""""determine the range of accelerometer data where cane 
position is considered upright when within this range"""
def activity_start_range(acc_x, acc_y, acc_z, upright_positions):
    low_x = upright_positions[0]-max(acc_x)/4
    high_x = upright_positions[0]+max(acc_x)/4
    low_y = upright_positions[1]-max(acc_y)/4
    high_y = upright_positions[1]+max(acc_y)/4
    low_z = upright_positions[2]-max(acc_z)/4
    high_z = upright_positions[2]+max(acc_z)/4
    return (low_x, high_x), (low_y, high_y), (low_z, high_z)

"""returns the timestamps that should start/refresh the scan for obstacles"""
def scan_time_range(timestamp, acc_x, acc_y, acc_z, thresholds): 
    stand_upright_times = []
    for i in range(0,len(timestamp)): 
        # start scanning after the initial 3s calibration time
        if timestamp[i] > 3: 
            # scan start when position of cane is upright (within the error range)
            if (acc_x[i]>thresholds[0][0] and acc_x[i]<thresholds[0][1]) and \
                (acc_y[i]>thresholds[1][0] and acc_y[i]<thresholds[1][1]) and \
                (acc_z[i]>thresholds[2][0] and acc_z[i]<thresholds[2][1]): 
                stand_upright_times.append(timestamp[i])
                
    # extract the initial timestamps during a range of time when cane is 
    # upright, which marks the start times of taking a step forward where scanning for obstacles
    scan_times = [] 
    scan_times.append(stand_upright_times[0])
    for i in range(0, len(stand_upright_times)-1): 
        ## assume at least 0.5s between steps 
        if stand_upright_times[i+1]-stand_upright_times[i]>0.5:
            scan_times.append(stand_upright_times[i])

    return scan_times

# Instantaneous real-time scanning
def scan_time_range_instant(timestamps, acc_x, acc_y, acc_z, thresholds,curr_time): 
    
    scan_start = False
    i=0
    for t in range(0,len(timestamps)):
        if curr_time == timestamps[t]:
            i=t
        elif t> 0:
            if curr_time < timestamps[t] and curr_time > timestamps[t-1]:
                i = t

    if (acc_x[i]>thresholds[0][0] and acc_x[i]<thresholds[0][1]) and \
                (acc_y[i]>thresholds[1][0] and acc_y[i]<thresholds[1][1]) and \
                (acc_z[i]>thresholds[2][0] and acc_z[i]<thresholds[2][1]):
        if (acc_x[i-1]>thresholds[0][0] and acc_x[i-1]<thresholds[0][1]) and \
                (acc_y[i-1]>thresholds[1][0] and acc_y[i-1]<thresholds[1][1]) and \
                (acc_z[i-1]>thresholds[2][0] and acc_z[i-1]<thresholds[2][1]): 
                scan_start = True                 
    

    return scan_start
"""Execute main"""

def cane_position_main():
    acc_dataFile = pd.read_csv('acc_cane.csv')
    timestamp = acc_dataFile['time'] 
    acc_x = acc_dataFile['gFx'] 
    acc_y = acc_dataFile['gFy'] 
    acc_z = acc_dataFile['gFz'] 

    # plot_acc(acc_x)
    # plot_acc(acc_y)
    # plot_acc(acc_z)

    # perform calculations of upright position and thresholds
    upright_positions = calibrate_cane_upright(timestamp, acc_x, acc_y, acc_z)
    thresholds = activity_start_range(acc_x, acc_y, acc_z, upright_positions)

    # time points that represent a new step and to refresh scanning for obstacles
    scan_time = scan_time_range(timestamp, acc_x, acc_y, acc_z, thresholds)
