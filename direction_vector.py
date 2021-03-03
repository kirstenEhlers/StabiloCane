import sys
import csv
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.interpolate.fitpack2 import InterpolatedUnivariateSpline
from scipy.signal import butter, lfilter
import scipy.integrate as integrate
from scipy import signal
from scipy import fft
import cane_position as cane
import Obstacle_Detection as obstacle
import pandas as pd

#Settings
num_skip = 2000
num_still = 1000
initial_yaw = 0
initial_roll = (90+87.3)/180*math.pi
initial_pitch = 0 #-8/180 * math.pi
method_trapezoid = False
threshold_factor = 0
num_ave =1000
f_high = 2
f_low = 10
height = 0.511 #m
step_length = 0.695 #m

class gyro_measure:
    def __init__(self,timestamp, angX, angY, angZ):
        self.time = timestamp
        self.roll = angX
        self.pitch= angY
        self.yaw = angZ

    def print(self):
        sys.stdout.write(f"Roll: {self.roll}\t Pitch: {self.pitch}\t Yaw: {self.yaw}\n")

class velocity:
    def __init__(self, timestamp, x_, y_, z_):
        self.time = timestamp
        self.x = x_
        self.y = y_
        self.z = z_
    def print(self):
        sys.stdout.write(f"{self.time}:\tx: {self.x}\ty: {self.y}\tz: {self.z}\n")
class acc_measure:
    def __init__(self,timestamp, aX, aY, aZ):
        self.time = timestamp
        self.accX = aX
        self.accY= aY
        self.accZ = aZ

class threshold_data:
    def __init__(self, q3, q1, average):
        self.min = average - (q3 - q1)*threshold_factor
        self.max = average + (q3 - q1)*threshold_factor

class threshold:

    def __init__(self, roll_list, pitch_list, yaw_list, aX_list, aY_list, aZ_list):
        self.roll= threshold_data(np.percentile(roll_list, 75), np.percentile(roll_list, 25), np.average(roll_list))
        self.pitch= threshold_data(np.percentile(pitch_list, 75), np.percentile(pitch_list, 25), np.average(pitch_list))
        self.yaw =threshold_data(np.percentile(yaw_list, 75), np.percentile(yaw_list, 25), np.average(yaw_list))
        self.X = threshold_data(np.percentile(aX_list, 75), np.percentile(aX_list, 25), np.average(aX_list))
        self.Y= threshold_data(np.percentile(aY_list, 75), np.percentile(aY_list, 25), np.average(aY_list))
        self.Z = threshold_data(np.percentile(aZ_list, 75), np.percentile(aZ_list, 25), np.average(aZ_list))
        self.X_ave = np.average(aX_list)
        self.Y_ave = np.average(aY_list)
        self.Z_ave = np.average(aZ_list)
        
def read_gyroscope_data(filename):
    with open(filename) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        list_gyro_data =[]
        time_subtract = 0
        for row in csv_reader:
            if line_count > num_skip:
                
                if(len(list_gyro_data) == 0):
                    time_subtract = float(row[0])

                data_point = gyro_measure(float(row[0])- time_subtract, float(row[1]), float(row[2]), float(row[3]))
                list_gyro_data.append(data_point)
                line_count += 1
            else:
                line_count += 1
        print(f'Processed {line_count} lines.')

        return [list_gyro_data, time_subtract]

def read_acc_data(filename, time_subtract):
    with open(filename) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        list_acc_data =[]
        
        for row in csv_reader:

            if line_count > 0:
                if (float(row[0]) - time_subtract) >=0:
                    
                    data_point = acc_measure(float(row[0])-time_subtract, float(row[1]), float(row[2]), float(row[3]))
                    list_acc_data.append(data_point)
                    line_count += 1
                else:
                    line_count += 1
            else:
                line_count += 1
        print(f'Processed {line_count} lines.')

        return list_acc_data

def read_data():
    [list_gyro_data, time_subtract] = read_gyroscope_data("./Data Files/hipstraught_gyrm.csv")
    sys.stdout.write("Read gyroscope data\n")
    list_acc_data = read_acc_data("./Data Files/hipstraught_lacm.csv", time_subtract)
    sys.stdout.write("Read acceleration data\n")

    return [list_gyro_data, list_acc_data]

def get_thresholds(gyro_data, acc_data):
    yaw_list = []
    roll_list = []
    pitch_list = []
    accX_list = []
    accY_list = []
    accZ_list = []

    for i in range(0, num_still):
        yaw_list.append(gyro_data[i].yaw)
        roll_list.append(gyro_data[i].roll)
        pitch_list.append(gyro_data[i].pitch)
        accX_list.append(acc_data[i].accX)
        accY_list.append(acc_data[i].accY)
        accZ_list.append(acc_data[i].accZ)
    
    return threshold(np.array(roll_list), np.array(pitch_list), np.array(yaw_list), np.array(accX_list), np.array(accY_list), np.array(accZ_list))
    
def filter_rates(acc, threshold):

    acc_out = 0
    if (acc < threshold.min or acc > threshold.max):
        acc_out = acc

    return acc_out

def track_angle(time_interval, thresholds, last_gyro, gyro, current_angles, curr_time):
    

    current_rate_yaw = filter_rates(gyro[2][0], thresholds.yaw)
    prev_rate_yaw = filter_rates(last_gyro[2][0], thresholds.yaw)
    current_yaw = current_angles.yaw + 0.5*(current_rate_yaw + prev_rate_yaw) * time_interval

    current_rate_pitch = filter_rates(gyro[1][0], thresholds.pitch)
    prev_rate_pitch = filter_rates(last_gyro[1][0],  thresholds.pitch)
    current_pitch = current_angles.pitch+ 1/2 * (current_rate_pitch + prev_rate_pitch) * time_interval

    current_rate_roll = filter_rates(gyro[0][0], thresholds.roll)
    prev_rate_roll = filter_rates(last_gyro[0][0], thresholds.roll)
    current_roll =current_angles.roll + 1/2 * (current_rate_roll + prev_rate_roll) * time_interval

    return gyro_measure(curr_time, current_roll, current_pitch, current_yaw)

def matrix_multiply(X, Y):

    result = [[0,0,0],
         [0,0,0],
         [0,0,0]]

    for i in range(len(X)):
   # iterate through columns of Y
        for j in range(len(Y[0])):
            # iterate through rows of Y
            for k in range(len(Y)):
                result[i][j] += X[i][k] * Y[k][j]
    return result

def rotate_accelerations(current_angles, a):

    rx = np.array([[1, 0, 0], [0, math.cos(current_angles.roll), -math.sin(current_angles.roll)], [0, math.sin(current_angles.roll), math.cos(current_angles.roll)]])
    ry = np.array([[math.cos(current_angles.pitch), 0, math.sin(current_angles.pitch)], [0, 1, 0], [-math.sin(current_angles.pitch), 0, math.cos(current_angles.pitch)]])
    rz = np.array([[math.cos(current_angles.yaw), -math.sin(current_angles.yaw), 0], [math.sin(current_angles.yaw), math.cos(current_angles.yaw), 0], [0, 0, 1]])

    a_new = matrix_multiply (rx, ry)
    a_new = matrix_multiply ( a_new, rz)
    a_new = matrix_multiply (a_new, a)

    r_acc = [[a_new[0][0]],[a_new[1][0]], [a_new[2][0]]]
    return r_acc

def filtered_acclerations(acc, thresholds):
    accX = filter_rates(acc[0][0], thresholds.X)
    accY = filter_rates(acc[1][0], thresholds.Y)
    accZ = filter_rates(acc[2][0],  thresholds.Z) 

    
    return [[accX],[accY], [accZ]]

def get_max(a,b):
    if a < b:
        return a
    
    return b

def track_velocity(acc, last_acc, current_velocity, last_time, time):
    time_interval = time-last_time

    X = acc
    Y = last_acc
    sum = [[X[i][j] + Y[i][j]  for j in range(len(X[0]))] for i in range(len(X))]
    change = 1/2 * np.multiply(sum,time_interval)

    X = change
    Y = current_velocity
    sum = [[X[i][j] + Y[i][j] for j in range(len(X[0]))] for i in range(len(X))]

    return sum

def sum_matrices(X, Y):
    return [[X[i][j] + Y[i][j] for j in range(len(X[0]))] for i in range(len(X))]

def get_next_rates_of_change(acc_i, gyro_i, gyro_list, acc_list):

    if acc_list[acc_i+1].time > gyro_list[gyro_i+1].time:
        curr_time = acc_list[acc_i+1].time

        add = 1
        sumX = acc_list[acc_i+add].accX
        sumZ = acc_list[acc_i+add].accZ
        sumY = acc_list[acc_i+add].accY

        while acc_list[acc_i+add].time == acc_list[acc_i+add +1].time:
            add +=1
            sumX += acc_list[acc_i+add+1].accX
            sumY += acc_list[acc_i+add+1].accY
            sumZ += acc_list[acc_i+add+1].accZ
        
        accX = sumX / add
        accY = sumY / add
        accZ = sumZ / add

        if gyro_i < 0:
            gyroX = 0
            gyroY =0
            gyroZ = 0
        else:
            gyroX = gyro_list[gyro_i].roll
            gyroY = gyro_list[gyro_i].pitch
            gyroZ = gyro_list[gyro_i].yaw

        acc_i += add

    elif acc_list[acc_i+1].time < gyro_list[gyro_i+1].time:
        curr_time = gyro_list[gyro_i+1].time

        add = 1
        sumX = gyro_list[gyro_i+add].pitch
        sumZ = gyro_list[gyro_i+add].yaw
        sumY = gyro_list[gyro_i+add].roll

       
        while gyro_list[gyro_i+add].time == gyro_list[gyro_i+add +1].time:
            add +=1
            sumX += gyro_list[gyro_i+add+1].roll
            sumY += gyro_list[gyro_i+add+1].pitch
            sumZ += gyro_list[gyro_i+add+1].yaw
        
        gyroX = sumX / add
        gyroY = sumY / add
        gyroZ = sumZ / add

        if acc_i < 0:
            accX = 0
            accY =0
            accZ = 0
        else:
            accX = acc_list[acc_i].accX
            accY = acc_list[acc_i].accY
            accZ= acc_list[acc_i].accZ

        gyro_i += add

    else:
        

        add = 1
        sumX = gyro_list[gyro_i+add].pitch
        sumZ = gyro_list[gyro_i+add].yaw
        sumY = gyro_list[gyro_i+add].roll

        while gyro_list[gyro_i+add].time == gyro_list[gyro_i+add +1].time:
            add = add + 1
            sumX += gyro_list[gyro_i+add+1].roll
            sumY += gyro_list[gyro_i+add+1].pitch
            sumZ += gyro_list[gyro_i+add+1].yaw
        
        gyroX = sumX / add
        gyroY = sumY / add
        gyroZ = sumZ / add
        gyro_i = gyro_i + add

        add = 1
        sumX = acc_list[acc_i+add].accX
        sumZ = acc_list[acc_i+add].accZ
        sumY = acc_list[acc_i+add].accY

        while acc_list[acc_i+add].time == acc_list[acc_i+add +1].time:
            add = add + 1
            sumX += acc_list[acc_i+add+1].accX
            sumY += acc_list[acc_i+add+1].accY
            sumZ += acc_list[acc_i+add+1].accZ
        
        accX = sumX / add
        accY = sumY / add
        accZ = sumZ / add
        acc_i = acc_i + add

        curr_time = gyro_list[gyro_i].time


    return [acc_i, gyro_i, accX, accY, accZ, gyroX, gyroY, gyroZ, curr_time]

def plot_data(list_a_x, list_a_y, list_a_z, list_g_x, list_g_y, list_g_z, list_vel_x, list_vel_y, list_vel_z, list_time, list_v_time, list_p_time):
    fig, (ax1, ax3,ax2) = plt.subplots(3)
    fig.suptitle('Vertically stacked subplots')
    ax1.plot(list_time, list_a_x, label = 'x')
    ax1.plot(list_time, list_a_y, label = 'y')
    ax1.plot(list_time, list_a_z, label = 'z')
    ax1.legend()
    plt.xlabel('Time (s)')
    ax1.set_ylabel('Acceleration (m/s/s)')
   

    ax2.plot(list_p_time, list_g_x, label = 'x')
    ax2.plot(list_p_time, list_g_y, label = 'y')
    ax2.plot(list_p_time, list_g_z, label = 'z')
    ax2.legend()
    
    ax2.set_ylabel('Position (m)')
 


    ax3.plot(list_v_time, list_vel_x, label = 'x')
    ax3.plot(list_v_time, list_vel_y, label = 'y')
    ax3.plot(list_v_time, list_vel_z, label = 'z')
    ax3.legend()
    
    ax3.set_ylabel('Velocity (m/s)')
    plt.show()

def plotSpectrum(y,Fs, name):
    """
    Plots a Single-Sided Amplitude Spectrum of y(t)
    """
    n = len(y) # length of the signal
    k = np.arange(n)
    T = n/Fs
    frq = k/T # two sides frequency range
    frq = frq[range(math.ceil(n/2))] # one side frequency range

    Y = fft.fft(y)/n # fft computing and normalization
    Y = Y[range(math.ceil(n/2))]
    
    plt.plot(frq,abs(Y),'r') # plotting the spectrum
    plt.xlabel('Freq (Hz)')
    plt.ylabel('|Y(freq)|')
    plt.title(name)
    plt.show()


def get_pos(list_vel_x, list_vel_y, list_vel_z, list_time):
    x =[]
    y = []
    z = []
    x_pos = 0
    y_pos = 0
    z_pos = 0
    last_time = 0
    for i in range(0, len(list_time)):
        
        x_pos = x_pos + list_vel_x[i]*(list_time[i]-last_time)
        y_pos = y_pos + list_vel_y[i]*(list_time[i]-last_time)
        z_pos = z_pos + list_vel_z[i]*(list_time[i]-last_time)

        x.append(x_pos)
        y.append(y_pos)
        z.append(z_pos)
        last_time = list_time[i]

    tt = list_time
    plt.plot(tt,x, label = 'vx')
    plt.plot(tt,y, label='vy')
    plt.plot(tt,z, label='vz')
    plt.legend()
    plt.show()

def track(gyro_list, acc_list, limit_high_pass, limit_low_pass, threshold_cane, cane_acc_x, cane_acc_y, cane_acc_z, timestamps):

    list_vel_x = []
    list_vel_y = []
    list_vel_z = []
    list_a_x = []
    list_a_y = []
    list_a_z = []
    list_time = []
    list_auf_x = []
    list_auf_y = []
    list_auf_z = []
    list_tuf = []

    accelerations_rolling_filter = [[],[],[],[]]
    velocities_rolling_filter = [[],[],[],[]]
    velocities = [[],[],[],[]]
    accelerations = [[],[],[],[]]
    positions = [[],[],[],[]]
    curr_time = 0
    acc_i =-1
    gyro_i = -1
    current_angles = gyro_measure(curr_time, initial_roll, initial_pitch, initial_yaw)
    last_acc = [[0],[0],[0]]
    acc = [[0],[0],[0]]
    acc_uf = acc
    thresholds = get_thresholds(gyro_list, acc_list)
    gyro = [[0],[0],[0]]
    max_time = get_max(gyro_list[-3].time, acc_list[-3].time)
    current_velocity = [[0],[0],[0]]
    fs = len(acc_list)/max_time
    i = 0
    
    while curr_time < max_time:
        [acc_i, gyro_i, accX, accY, accZ, gyroX, gyroY, gyroZ, t]=get_next_rates_of_change(acc_i, gyro_i, gyro_list, acc_list)
        last_time = curr_time
        curr_time = t
        

        last_gyro = gyro
        acc = [[accX], [accY], [accZ]]
        gyro = [[gyroX], [gyroY], [gyroZ]]
        current_angles = track_angle(curr_time-last_time, thresholds, last_gyro, gyro, current_angles, curr_time)
        
        acc_uf = rotate_accelerations(current_angles, acc)
        acc = filtered_acclerations(acc, thresholds)
        acc = rotate_accelerations(current_angles, acc)
        list_auf_x.append(acc_uf[0][0])
        list_auf_y.append(acc_uf[1][0])
        list_auf_z.append(acc_uf[2][0])
        list_tuf.append(curr_time)
        #acc_sum = sum_matrices(acc_sum, acc)
        

        i = i + 1
        accelerations_rolling_filter[0].append(acc[0][0])
        accelerations_rolling_filter[1].append(acc[1][0])
        accelerations_rolling_filter[2].append(acc[2][0])
        accelerations_rolling_filter[3].append(curr_time)

        current_velocity = track_velocity(acc, last_acc,current_velocity, last_time,curr_time)
        vel = velocity(curr_time, current_velocity[0][0], current_velocity[1][0], current_velocity[2][0])
        last_acc = acc

        velocities_rolling_filter[0].append(vel.x)
        velocities_rolling_filter[1].append(vel.y)
        velocities_rolling_filter[2].append(vel.z)
        velocities_rolling_filter[3].append(curr_time)

        
        trigger_analysis = True
        if curr_time == max_time:
            trigger_analysis = cane.scan_time_range_instant(timestamps, cane_acc_x, cane_acc_y, cane_acc_z, threshold_cane,curr_time)
        else:
            trigger_analysis = False

        #will be condition for acceleration - To do
        if trigger_analysis:
            times = velocities_rolling_filter[3]

            for k in range(0,3):
                #a = interp1d(times, accelerations_rolling_filter[k], kind='nearest')
                sig = accelerations_rolling_filter[k]
                sos = signal.butter(10,limit_low_pass, 'low', fs=500, output='sos')
                af = sig#signal.sosfilt(sos, sig)
                #finterp = InterpolatedUnivariateSpline(times,a_filtered, k = 1)
                #v = [finterp.integral(0, t) for t in times]
                sig = velocities_rolling_filter[k]
                sos = signal.butter(10, [limit_high_pass, limit_low_pass], 'bandpass', fs=fs, output='sos')
                vf = signal.sosfilt(sos, sig)

                finterp = InterpolatedUnivariateSpline(times,vf, k = 1)
                p = [finterp.integral(0, t) for t in times]
                sig = p
                sos = signal.butter(10, limit_high_pass, 'high', fs=fs, output='sos')
                pf = signal.sosfilt(sos, sig)

                for pos in pf:
                    positions[k].append(pos)

                for vel in vf:
                    velocities[k].append(vel)

                for accel in af:
                    accelerations[k].append(accel)

                if k == 0:
                    for t in times:
                        velocities[3].append(t)
                        positions[3].append(t)
                        accelerations[3].append(t)

            #SCANNING PARAMETERS FOR AKANKSHA
            ax = [velocities[0][-1]][0]
            ay = [velocities[1][-1]][0]
            az = [velocities[2][-1]][0]
            curr_vel = [ax, ay, az]
            print("Obstacle detection in progress...")
            print(obstacle.detect_obstacle(curr_vel))
            #height
            #step_length

            #resetting
            for j in range(0,4):
                accelerations_rolling_filter[j] = []
                velocities_rolling_filter[j] = []
            i = 0
            acc_sum = [[0],[0],[0]]
            last_acc = acc
        else:
            f = 0


        
    list_vel_x = velocities[0]
    list_vel_y = velocities[1]
    list_vel_z = velocities[2]
    list_v_time = velocities[3]

    list_a_x = accelerations[0]
    list_a_y = accelerations[1]
    list_a_z = accelerations[2]
    list_time = accelerations[3]


    list_p_x = positions[0]
    list_p_y = positions[1]
    list_p_z = positions[2]
    list_p_time = positions[3]
  
    plot_data(list_a_x, list_a_y, list_a_z, list_p_x, list_p_y, list_p_z, list_vel_x, list_vel_y, list_vel_z, list_time, list_v_time, list_p_time)
    
    
def main():
    sys.stdout.write("hello world\n")
    [list_gyro_data, list_acc_data] = read_data()

    #calibration steps for scanning times
    acc_dataFile = pd.read_csv('acc_cane.csv')
    timestamps = acc_dataFile['time'] 
    cane_acc_x = acc_dataFile['gFx'] 
    cane_acc_y = acc_dataFile['gFy'] 
    cane_acc_z = acc_dataFile['gFz'] 
    upright_positions = cane.calibrate_cane_upright(timestamps, cane_acc_x, cane_acc_y, cane_acc_z)
    thresholds_cane = cane.activity_start_range(cane_acc_x, cane_acc_y, cane_acc_z, upright_positions)

    
    track(list_gyro_data, list_acc_data,0.0102,3, thresholds_cane, cane_acc_z, cane_acc_y, cane_acc_z, timestamps)


main()

