from numpy import timedelta64

import math
from scipy import signal
import numpy as np
import matplotlib.pyplot as plt
import csv

class acc_measure:
    def __init__(self,timestamp, aX, aY, aZ):
        self.time = timestamp
        self.accX = aX
        self.accY= aY
        self.accZ = aZ

def read_acc_data(filename, time_subtract):
    with open(filename) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        list_acc_data =[]
        
        for row in csv_reader:

            if line_count > 0:
                if (True) >=0:
                    
                    data_point = acc_measure(0, float(row[1]), float(row[2]), float(row[3]))
                    list_acc_data.append(data_point)
                    line_count += 1
                else:
                    line_count += 1
            else:
                line_count += 1
        print(f'Processed {line_count} lines.')

        return list_acc_data


list_acc_data = read_acc_data("./Data Files/SNR_chest.csv", 0)

total_time = 14.018#foot12.953#hip13.021
#front15.05
FS = math.floor(len(list_acc_data)/total_time)
x = []
y = []
z = []
t=[]
SNR = []
for a in list_acc_data:
    x.append(a.accX)
    y.append(a.accY)
    z.append(a.accZ)
    t.append(a.time)

for sig in [x,y,z]:
    sos = signal.butter(10,4, 'low', fs=FS, output='sos')
    af = signal.sosfilt(sos, sig)
    sos = signal.butter(10,4, 'high', fs=FS, output='sos')
    an = signal.sosfilt(sos, sig)
    plt.plot(sig)
    plt.plot(af)
    plt.show()

    sum = 0
    #Calculate power
    for v in af:
        sum = sum + (v * v)

    power_f = sum/total_time

    sum = 0
    #Calculate power
    for v in an:
        sum = sum + (v * v)

    power_n = sum/total_time

    SNR.append(power_f/power_n)

for d in SNR:
    print(d)



