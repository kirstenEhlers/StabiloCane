import sys
import csv
import numpy as np
import math

#Settings
num_skip = 300
num_still = 1000
initial_yaw = 0
initial_roll = -0.4157
initial_pitch = -0.1146
method_trapezoid = False
threshold_factor = 1.5

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
    def __init__(self, q1, q3, average):
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
    [list_gyro_data, time_subtract] = read_gyroscope_data("./Data Files/gyro_input.csv")
    sys.stdout.write("Read gyroscope data\n")
    list_acc_data = read_acc_data("./Data Files/lacm_input.csv", time_subtract)
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
    
def filter_rates(acc, acc_prev, acc_next, threshold):

    acc_out = 0
    if ((acc < threshold.min or acc > threshold.max) and (acc_prev > threshold.max or acc_prev < threshold.min)) and (acc_next > threshold.max or acc_next < threshold.min):
        acc_out = acc

    return acc_out

def track_angle(gyro_list, thresholds, i, current_angles):
    time_interval = gyro_list[i].time - gyro_list[i-1].time

    current_rate_yaw = filter_rates(gyro_list[i].yaw, gyro_list[i-1].yaw, gyro_list[i+1].yaw, thresholds.yaw)
    prev_rate_yaw = filter_rates(gyro_list[i-1].yaw, gyro_list[i-2].yaw, gyro_list[i].yaw, thresholds.yaw)
    current_yaw = current_angles.yaw + 1/2 * (current_rate_yaw + prev_rate_yaw) * time_interval

    current_rate_pitch = filter_rates(gyro_list[i].pitch, gyro_list[i-1].pitch, gyro_list[i+1].pitch, thresholds.pitch)
    prev_rate_pitch = filter_rates(gyro_list[i-1].pitch, gyro_list[i-2].pitch, gyro_list[i].pitch, thresholds.pitch)
    current_pitch = current_angles.pitch+ 1/2 * (current_rate_pitch + prev_rate_pitch) * time_interval

    current_rate_roll = filter_rates(gyro_list[i].roll, gyro_list[i-1].roll, gyro_list[i+1].roll, thresholds.roll)
    prev_rate_roll = filter_rates(gyro_list[i-1].roll, gyro_list[i-2].roll, gyro_list[i].roll, thresholds.roll)
    current_roll =current_angles.roll + 1/2 * (current_rate_roll + prev_rate_roll) * time_interval

    return gyro_measure(gyro_list[i].time, current_roll, current_pitch, current_yaw)

def rotate_accelerations(current_angles, a):

    rx = np.array([[1, 0, 0], [0, math.cos(current_angles.roll), -math.sin(current_angles.roll)], [0, math.sin(current_angles.roll), math.cos(current_angles.roll)]])
    ry = np.array([[math.cos(current_angles.pitch), 0, math.sin(current_angles.pitch)], [0, 1, 0], [-math.sin(current_angles.pitch), 0, math.cos(current_angles.pitch)]])
    rz = np.array([[math.cos(current_angles.yaw), -math.sin(current_angles.yaw), 0], [math.sin(current_angles.yaw), math.cos(current_angles.yaw), 0], [0, 0, 1]])


    return rx *ry * rz * a

def filtered_acclerations(acc, acc_prev, acc_next, thresholds):
    accX = filter_rates(acc.accX, acc_prev.accX, acc_next.accX, thresholds.X)
    accY = filter_rates(acc.accY, acc_prev.accY, acc_next.accY, thresholds.Y)
    accZ = filter_rates(acc.accZ, acc_prev.accZ, acc_next.accZ, thresholds.Z)


    sys.stdout.write(f"ACC {accX} {accY} {accZ}\n")
    return [[accX],[accY], [accZ]]

def get_max(a,b):
    if a > b:
        return a
    
    return b

def track_velocity(acc, last_acc, current_velocity, last_time, time):
    time_interval = time-last_time
    change = 1/2 * (acc + last_acc) * time_interval

    return current_velocity + change

def track(gyro_list, acc_list):

    curr_time = 0
    current_angles = gyro_measure(curr_time, initial_roll, initial_pitch, initial_yaw)
    acc = [[0],[0],[0]]
    thresholds = get_thresholds(gyro_list, acc_list)
    
    max_time = get_max(gyro_list[-1].time, acc_list[-1].time)

    i = 0
    j = 0
    last_time = 0
    current_velocity = [[0],[0],[0]]
    time = 0

    while i < len(gyro_list)-2 and j < len(acc_list)-2:
        last_acc = acc
        last_time = time
        current_angles = track_angle(gyro_list, thresholds, i, current_angles)
        acc = filtered_acclerations(acc_list[j], acc_list[j-1], acc_list[j+1], thresholds)
        acc = rotate_accelerations(current_angles, acc)
        time = get_max(acc_list[j].time, gyro_list[i].time)
        current_velocity = track_velocity(acc, last_acc,current_velocity, last_time,time)
        vel = velocity(time, current_velocity[0][0], current_velocity[1][0], current_velocity[2][0])
        vel.print()

        if(gyro_list[i+1].time < acc_list[j+1].time):
            i += 1
        else:
            j += 1

    
def main():
    sys.stdout.write("hello world\n")
    [list_gyro_data, list_acc_data] = read_data()
    track(list_gyro_data, list_acc_data)

main()
