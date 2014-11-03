## we detect peaks (ball going by) using the methodologies described here:
## http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2631518/

##for graph
from pyqtgraph.Qt import QtGui, QtCore
import numpy as np
import pyqtgraph as pg
from pyqtgraph.dockarea import *

import sys, os, serial, datetime

from threading import Thread, Lock

mutex = Lock()

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    from math import factorial   
    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')
          
class capture(QtCore.QThread):
    def capture_line(self):
        global arduino_data, arduino_matrix
        line = self.ser.readline()
        mutex.acquire()
        arduino_line = line.rstrip().split(" ")
        if (len(arduino_line) == 4):
            try:
                arduino_data = map(int, arduino_line)
                arduino_matrix = np.roll(arduino_matrix, -1, axis=0)
                arduino_matrix[-1] = arduino_data
            except ValueError:
                print "Bad value read from serial..."
                arduino_data = [0,0,0,0]
            else:
                arduino_data = [0,0,0,0]
        mutex.release()
        
    def run(self):
        global arduino_data, arduino_matrix
        print "Starting arduino capture"
        self.ser = serial.Serial("/dev/ttyACM1", 115200, timeout=None)
        print "Init sequence..."
        for x in range(0,1000):
            self.capture_line()
        mutex.acquire()
        print "Signal 0 average: %d" % np.average(arduino_matrix[:,0][:1000])
        arduino_matrix[:,0] = np.average(arduino_matrix[:,0][-1000:])
        arduino_matrix[:,1] = np.average(arduino_matrix[:,1][-1000:])
        arduino_matrix[:,2] = np.average(arduino_matrix[:,2][-1000:])
        arduino_matrix[:,3] = np.average(arduino_matrix[:,3][-1000:])
        mutex.release()
        print "Monitoring sequence..."
        while (1):
            self.capture_line()

def reset():
    global arduino_matrix
    mutex.acquire()
    arduino_matrix = np.zeros(shape=(100001,4))
    mutex.release()
    
##register app and config window
app = QtGui.QApplication([])
win = QtGui.QMainWindow()
area = DockArea()
win.setCentralWidget(area)
win.resize(1000,900)
win.setWindowTitle('Shotdoc Monitor')
pg.setConfigOptions(antialias=True)

d1 = Dock("Control")
d2 = Dock("Sensors", size=(300,500))
d6 = Dock("Shot overlay", size=(300,200))
d7 = Dock("Shot Angle", size=(300,200))
area.addDock(d1)
area.addDock(d2)
area.addDock(d6)
area.addDock(d7,'right', d6)

w1 = pg.LayoutWidget()

resetBtn = QtGui.QPushButton('Reset Graphs')
resetBtn.clicked.connect(reset)

w1.addWidget(resetBtn)
d1.addWidget(w1)

p1 = pg.PlotWidget(title="Sensors")
p1.setRange(xRange=(0,10000),yRange=(500,1024))
curve1gas = p1.plot(pen='y')
curve2gas = p1.plot(pen='g')
curve3gas = p1.plot(pen='b')
curve4gas = p1.plot(pen='w')
d2.addWidget(p1)

p6 = pg.PlotWidget(title="Shot overlay")
d6.addWidget(p6)

p7 = pg.PlotWidget(title="Shot Angle")
d7.addWidget(p7)

win.show()

def signal_analysis(signal_array):
    if (len(signal_array) > 0):
        ##smoooooooth array
        savitzky = savitzky_golay(signal_array, window_size=31, order=4)
        ##stick the badboys on the baseline
        minimum = np.amin(savitzky)
        baseline = np.zeros_like(savitzky)
        baseline.fill(minimum)
        savitzky = savitzky# - baseline
        baseline = np.zeros_like(signal_array)
        baseline.fill(minimum)
        orig = signal_array - baseline
    else:
        orig = signal_array
        savitzky = np.zeros(shape=(1000))
    return (orig, savitzky)

def group_peaks(vals):
    peaks = []
    begin_peak = -1
    for idx, val in enumerate(vals):
        if val:
            if begin_peak >= 0:
                # we are already tracing a peak, so continue
                continue
            else:
                # begin trace
                begin_peak = idx
        else:
            if begin_peak >= 0 and idx-begin_peak>20:
                #end of a peak, add tuple to list if long enough
                peaks.append((begin_peak, idx))
                begin_peak = -1
            else:
                #chugging through no peak
                begin_peak = -1
    return peaks

def shot_find(sensor1, sensor2, sensor3, sensor4):
    #peak detect by finding aggressive slopes in cross-section
    slope_thresh = .5
    period = 5
    dotThresh = 4000
    if sensor1.any():
        dot = sensor1 * sensor2 * sensor3 * sensor4
    else:
        dot = np.zeros(shape=(arduino_matrix[-1].shape))
    #normalize
    dot[dot<dotThresh] = 0
    dot[dot>0] = 1
    #find peak indicies
    peak_groups = group_peaks(dot)

    #use peak groups to slice
    return peak_groups

def update():
    global arduino_data, arduino_matrix
    mutex.acquire()
    
    (original1, savitzky1) = signal_analysis(arduino_matrix[:,0])
    curve1gas.setData(savitzky1)
    ##lr1.setRegion([average + .1*average, average - .1*average])

    (original2, savitzky2) = signal_analysis(arduino_matrix[:,1])
    curve2gas.setData(savitzky2)

    (original3, savitzky3) = signal_analysis(arduino_matrix[:,2])
    curve3gas.setData(savitzky3)

    (original4, savitzky4) = signal_analysis(arduino_matrix[:,3])
    curve4gas.setData(savitzky4)
    
    slices = shot_find(savitzky1, savitzky2, savitzky3, savitzky4)
    #draw last shot only
    if slices:
        p6.clear()
        curve6_1 = p6.plot(pen='y')
        curve6_2 = p6.plot(pen='g')
        curve6_3 = p6.plot(pen='b')
        curve6_4 = p6.plot(pen='w')
        #draw curves
        curve6_1.setData(savitzky1[slices[-1][0]:slices[-1][1]])
        curve6_2.setData(savitzky2[slices[-1][0]:slices[-1][1]])
        curve6_3.setData(savitzky3[slices[-1][0]:slices[-1][1]])
        curve6_4.setData(savitzky4[slices[-1][0]:slices[-1][1]])
        #box peaks
        peak1 = np.argmax(savitzky1[slices[-1][0]:slices[-1][1]])
        peak2 = np.argmax(savitzky2[slices[-1][0]:slices[-1][1]])
        peak3 = np.argmax(savitzky3[slices[-1][0]:slices[-1][1]])
        peak4 = np.argmax(savitzky4[slices[-1][0]:slices[-1][1]])
        arrow1 = pg.ArrowItem(brush='y')
        arrow2 = pg.ArrowItem(brush='g')
        arrow3 = pg.ArrowItem(brush='b')
        arrow4 = pg.ArrowItem(brush='w')
        arrow1.setPos(peak1, savitzky1[slices[-1][0]:slices[-1][1]][peak1])
        arrow2.setPos(peak2, savitzky2[slices[-1][0]:slices[-1][1]][peak2])
        arrow3.setPos(peak3, savitzky3[slices[-1][0]:slices[-1][1]][peak3])
        arrow4.setPos(peak4, savitzky4[slices[-1][0]:slices[-1][1]][peak4])
        p6.addItem(arrow1)
        p6.addItem(arrow2)
        p6.addItem(arrow3)
        p6.addItem(arrow4)
        #now fill in angle
        p7.clear()
        #peaks are time values of when maxed
        #only count shots that are falling top to bottom
        if (min([peak1,peak2]) < min([peak3, peak4])):
            print "Shot angle detected!"
            top_sensor_val = sum([savitzky1[slices[-1][0]:slices[-1][1]][peak1], savitzky2[slices[-1][0]:slices[-1][1]][peak2]]) / 2
            bot_sensor_val = sum([savitzky3[slices[-1][0]:slices[-1][1]][peak3], savitzky4[slices[-1][0]:slices[-1][1]][peak4]]) / 2
            direction_curve = p7.plot(pen='w')
            direction_curve.setData([{'x': top_sensor_val, 'y': 1},
                                     {'x': bot_sensor_val, 'y': 0}])
        else:
            print "Not a valid shot"
        
    mutex.release()

thread = capture()

global arduino_data
#pointer for all plots
mutex.acquire()
arduino_data = []
mutex.release()
arduino_matrix = np.empty(shape=(10000,4))
arduino_matrix[:] = 555

timer = QtCore.QTimer()
timer.timeout.connect(update)
timer.start(10)

##start main loop
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        thread.start()
        QtGui.QApplication.instance().exec_()

