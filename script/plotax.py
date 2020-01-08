#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Update a simple plot as rapidly as possible to measure speed.
"""

from pyqtgraph.Qt import QtGui, QtCore
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

for meas in fileMeas:
#			print(meas)
    fileHdf5 = h5py.File(meas, "r")
    zR = fileHdf5["/Physical"].attrs.get("z")

    if (mode == 'theta'):
        aData = fileHdf5['Field/data']
        aData = aData/zR
    elif (mode == 'vA'):
        aData  = fileHdf5['Dev/data'].value
    elif mode == 'eA':
        aData = (fileHdf5['Field/data'].value)**2 + (fileHdf5['Dev/data'].value)**2
        # missing gradients

    data.append(aData)
    zdat.append(zR)
    fileHdf5.close()
    # size = size + 1
    # print(meas,zR)

data  = np.array(data)
zdat  = np.array(zdat)
sizeN = len(data[0])
sizeT = len(data)

print('sizeN %d',sizeN)
################################################################################

app = QtGui.QApplication([])

p = pg.plot()
p.setWindowTitle('pyqtgraph Axiton evol')
p.setRange(QtCore.QRectF(0, -5, sizeN, 10))
p.setLabel('bottom', 'Radius', units='L1')
curve = p.plot()


ptr = 0
lastTime = time()
fps = None
def update():
    global curve, data, ptr, p, lastTime, fps
    # curve.setData(data[ptr%10])
    curve.setData(data[ptr%sizeT])
    ptr += 1
    now = time()
    dt = now - lastTime
    lastTime = now
    if fps is None:
        fps = 1.0/dt
    else:
        s = np.clip(dt*3., 0, 1)
        fps = fps * (1-s) + (1.0/dt) * s
    # p.setTitle('z = %.3f %d (%0.2f fps)' % (zdat[ptr%sizeT],ptr%sizeT, fps))
    p.setTitle('z = %.3f     (%d)' % (zdat[ptr%sizeT],ptr%sizeT))
    app.processEvents()  ## force complete redraw for every plot
timer = QtCore.QTimer()
timer.timeout.connect(update)
timer.start(50)



## Start Qt event loop unless running in interactive mode.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
