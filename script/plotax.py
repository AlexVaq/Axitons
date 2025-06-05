#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Update a simple plot as rapidly as possible to measure speed.
"""

from pyqtgraph.Qt import QtGui, QtCore, QtWidgets
import numpy as np
import pyqtgraph as pg
from pyqtgraph.ptime import time
# JAVI ADDED
import sys, os, h5py, re

################################################################################
# data = np.random.normal(size=(50,10000))
data = []
zdat = []

print('modes: theta [default], vtheta, saxion, vsaxion, saxion, eA, sP, real, imag')
mode = 'theta'
if len(sys.argv) == 2:
    if (sys.argv[-1] == 'eA'):
        mode = 'eA'
        print('Axion energy')
    if (sys.argv[-1] == 'vA'):
        mode = 'vA'
        print('Axion velocity')

prefileMeas = sorted([x for x in [y for y in os.listdir("./")] if re.search("axiton.[0-9]{5}$", x)])

fileMeas = []
for maes in prefileMeas:
	try:
		with h5py.File(maes, 'r') as f:
			if ('Field' in f) :
				fileMeas.append(maes)
	except:
		print('Error opening file: %s'%maes)


fileHdf5 = h5py.File(fileMeas[0], "r")
ctheta     = fileHdf5['Field/data']
r  = np.linspace(0,fileHdf5["/Physical"].attrs.get("Physical Length"),len(ctheta))
dx = r[1]-r[0]
fileHdf5.close()
del ctheta

data = np.zeros((len(fileMeas),len(r)))
zdat = np.zeros((len(fileMeas)))
mdat = np.zeros((len(fileMeas)))
i = 0
for meas in fileMeas:
#			print(meas)
    fileHdf5 = h5py.File(meas, "r")
    zR = fileHdf5["/Physical"].attrs.get("z")
    R  = fileHdf5["/Physical"].attrs.get("R")
    mA = fileHdf5["/Physical"].attrs.get("Axion mass")


    if (mode == 'theta'):
        data[i] = fileHdf5['Field/data'][:]/R
    elif (mode == 'vA'):
        data[i] = fileHdf5['Dev/data'][:]
    elif mode == 'eA':
        chi = mA*mA*R*R*R*R
        if mA == 0:
            den = 1/(R*R*R*R)
        else :
            den = mA/(R*R*R)
        ctheta     = fileHdf5['Field/data'].value
        grad = ctheta*0.0
        grad[1:-1] = (ctheta[2:]-ctheta[:-2])/(2*dx)
        grad[0]    = (ctheta[1]-ctheta[0])/(dx)
        grad[-1]   = (ctheta[-1]-ctheta[-2])/(dx)
        data[i]    = np.log10( (chi*(1 - np.cos(ctheta/R)) + 0.5*(fileHdf5['Dev/data'].value -ctheta/R)**2 + 0.5*grad**2)/(R*R*R*R)/den )
        # missing gradients

    # data.append(aData)
    zdat[i] = zR
    mdat[i] = mA*R
    fileHdf5.close()
    i = i + 1
    # size = size + 1
    # print(meas,zR)

# data  = np.array(data)
# zdat  = np.array(zdat)
sizeN = len(data[0])
sizeT = len(data)

print('sizeN %d',sizeN)
################################################################################

class MyWindow(pg.GraphicsWindow):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.timer = QtCore.QTimer()
        self.plot = self.addPlot()
        self.plot.setWindowTitle('pyqtgraph Axiton evol')
        self.plot.setLabel('bottom', 'Radius', units='L1')
        self.curve = self.plot.plot()
        self.timer.timeout.connect(self.update)
        self.tStep = 50
        self.fps   = None
        self.timer.start(self.tStep)
        self.lastTime = time()
        self.ptr   = 0
        self.pause = False
        QtWidgets.qApp.installEventFilter(self)

    def eventFilter(self, source, event):
      if event.type() == QtCore.QEvent.KeyPress:
        key = event.key()

        if key == QtCore.Qt.Key_Space:
            if self.pause:
                self.pause = False
                self.timer.start(self.tStep)
            else:
                self.pause = True
                self.timer.stop()
        elif key == QtCore.Qt.Key_N:
                self.tStep = self.tStep + 30
                if self.tStep > 500:
                    self.tStep = 500
                self.timer.setInterval(self.tStep)
                print("Increased step to", self.tStep)
        elif key == QtCore.Qt.Key_M:
                self.tStep = self.tStep - 30
                if self.tStep < 30:
                    self.tStep = 20
                self.timer.setInterval(self.tStep)
                print("Reduced step to", self.tStep)
      return super().eventFilter(source, event)

    def update(self):
        self.curve.setData(r, data[self.ptr%sizeT])
        # self.curve.setData(data[self.ptr%sizeT,0],data[self.ptr%sizeT,1])
        self.ptr += 1
        now = time()
        dt = now - self.lastTime
        self.lastTime = now
        if self.fps is None:
            self.fps = 1.0/dt
        else:
            s = np.clip(dt*3., 0, 1)
            self.fps = self.fps * (1-s) + (1.0/dt) * s
        self.plot.setTitle('z = %.3f     (%d) m_c = %.3f' % (zdat[self.ptr%sizeT],self.ptr%sizeT, mdat[self.ptr%sizeT]))
        QtGui.QApplication.processEvents()


app = pg.mkQApp()
win = MyWindow()

## Start Qt event loop unless running in interactive mode.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        #QtGui.QApplication.instance().exec_()
        app.exec_()
