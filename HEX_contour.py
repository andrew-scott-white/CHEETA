# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 15:58:54 2021

@author: white
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp

# specify location of Nernst Data
path = r'C:\Users\white\Documents\Research\GTL\Sheets\HEXDragCalculator.xlsm'
sheet = 'ContourData'

# Nernst Data
read = pd.read_excel(path, sheet_name=sheet,header=0)   
table = np.asarray(read) 

D = np.asarray(table[1:11, 1:7])
Q = np.asarray(table[0:1,1:7])
A = np.asarray(table[1:11,0:1])

f1, f2 = np.meshgrid(Q,A)


mp.pyplot.contourf(f1,f2, D, levels=[0, 0.05, .1, .15, .2, .3, .4, .5, .6, .7, .8, .9, 1], cmap='RdGy')
plt.colorbar()
mp.pyplot.ylim(10,3000)


mp.pyplot.title('Contour of HEX Drag as a Fraction of Total Aircraft Drag')
mp.pyplot.xlabel('Thermal Load (kW)')
mp.pyplot.ylabel('HEX Face Area (sq in)')


#x, y = np.meshgrid(np.arange(10),np.arange(10))
#z = np.sqrt(x**2 + y**2)
#cs = plt.contourf(x,y,z,levels=[.1, .15, .2, .25, .3, .4, .5, .7, .9, 1.1])#

#proxy = [plt.Rectangle((0,0),1,1,fc = pc.get_facecolor()[0]) for pc in cs.collections]

#plt.legend(proxy, [.1, .15, .2, .25, .3, .4, .5, .7, .9, 1.1])
#plt.legend(loc='lower left')
#plt.show(plot)