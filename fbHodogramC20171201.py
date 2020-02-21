# -*- coding: utf-8 -*-
'''
Created on Wen Dec 26 23:45:30 2019

@author: Mengmeng Wang
'''

'''
Minimum Variance Analysis
Magnetic field
Version 0
case 20171201_
'''

'''
!!!Note
'''

import cdflib
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import os
###
import time

'''
function MVA:===========================================================================================================
'''


def mva(bx, by, bz):
    mm = np.ones((3, 3))
    mm[0][0] = np.mean(bx * bx) - np.mean(bx) * np.mean(bx)
    mm[0][1] = np.mean(bx * by) - np.mean(bx) * np.mean(by)
    mm[0][2] = np.mean(bx * bz) - np.mean(bx) * np.mean(bz)
    mm[1][0] = np.mean(by * bx) - np.mean(by) * np.mean(bx)
    mm[1][1] = np.mean(by * by) - np.mean(by) * np.mean(by)
    mm[1][2] = np.mean(by * bz) - np.mean(by) * np.mean(bz)
    mm[2][0] = np.mean(bz * bx) - np.mean(bz) * np.mean(bx)
    mm[2][1] = np.mean(bz * by) - np.mean(bz) * np.mean(by)
    mm[2][2] = np.mean(bz * bz) - np.mean(bz) * np.mean(bz)
    [D, V] = np.linalg.eig(mm)
    lambda1 = np.min(D)
    Imin = np.where(D == lambda1)
    lambda3 = np.max(D)
    Imax = np.where(D == lambda3)
    kLoop = np.arange(0, 3)
    for k in kLoop:
        if D[k] != lambda1 and D[k] != lambda3:
            lambda2 = D[k]
            Imid = k
    Imin = Imin[0]
    Imin = Imin[0]
    Imax = Imax[0]
    Imax = Imax[0]
    N = V[:, Imin]
    M = V[:, Imid]
    L = V[:, Imax]
    DNew = [D[Imin], D[Imid], D[Imax]]
    VNew = [N, M, L]
    return (lambda1, lambda2, lambda3, N, M, L, DNew, VNew)


'''
!!!Note:
arrange eigenvalues from small to large
lambda1, lambda2, lambda3
corresponding to eigenvectors
N, M, L
corresponding to variation of magnetic field 
minimum, intermediate, maximum
!!!Note
'''
'''
function MVA:===========================================================================================================
'''

pyStartTime = time.time()

timeInput1 = [2017, 12, 1, 14, 40, 10, 000]
timeInput2 = [2017, 12, 1, 14, 41, 40, 000]

fileNameFgm = cdflib.CDF('D:/data/mms/mms1/fgm/brst/l2/2017/12/01/mms1_fgm_brst_l2_20171201143933_v5.114.0.cdf')
infoFgm = fileNameFgm.cdf_info()
epochFgm = fileNameFgm.varget('Epoch')
bGSEVec = fileNameFgm.varget('mms1_fgm_b_gse_brst_l2')

epochStartType2000 = cdflib.cdfepoch.compute_tt2000(timeInput1)
epochStopType2000 = cdflib.cdfepoch.compute_tt2000(timeInput2)
timeTuple1 = np.where(epochFgm > epochStartType2000)
timeTuple2 = np.where(epochFgm < epochStopType2000)
timeArr1 = timeTuple1[0]
timeArr2 = timeTuple2[0]
timeArrToCal = np.intersect1d(timeArr1, timeArr2)

bxGSEVec = bGSEVec[:, 0]
byGSEVec = bGSEVec[:, 1]
bzGSEVec = bGSEVec[:, 2]
btGSEVec = bGSEVec[:, 3]

bxGSEToCal = bxGSEVec[timeArrToCal]
byGSEToCal = byGSEVec[timeArrToCal]
bzGSEToCal = bzGSEVec[timeArrToCal]
btGSEToCal = btGSEVec[timeArrToCal]
epochToCal = epochFgm[timeArrToCal]

[lambda1, lambda2, lambda3, N, M, L, eig, arrRot] = mva(bxGSEToCal, byGSEToCal, bzGSEToCal)

bGSEToRot = np.column_stack((bxGSEToCal, byGSEToCal, bzGSEToCal))
bLMN = np.dot(bGSEToRot, arrRot)

fig = plt.figure(figsize=(10, 8))
fig.subplots_adjust(top=0.94)
fig.subplots_adjust(bottom=0.07)
fig.subplots_adjust(left=0.13)
fig.subplots_adjust(right=0.93)
gs = GridSpec(4, 2, figure=fig)
ax0 = fig.add_subplot(gs[0, :])
ax1 = fig.add_subplot(gs[1, :])
ax2 = fig.add_subplot(gs[2, :])
ax3 = fig.add_subplot(gs[-1, 0])
ax4 = fig.add_subplot(gs[-1, -1])
# plt.subplots_adjust(wspace = 0, hspace = 0)


xlim = [min(epochToCal), max(epochToCal)]
colors = ['blue', 'green', 'red']
labelsGSE = ['bx', 'by', 'bz']
labelsLMN = ['bn', 'bm', 'bl']
kLoop = np.arange(0, 3)

for k in kLoop:
    ax0.plot(epochToCal, bGSEToRot[:, k], linewidth=0.8, color=colors[k], label=labelsGSE[k])
    ax0.legend(loc=2, fontsize=8, markerscale=1.0, handlelength=2, handleheight=0.5)
    #    ax1.plot(epochToCal, bLMN[:,k], linewidth = 0.8, color = colors[k], label = labelsLMN[k])

ax1.plot(epochToCal, bLMN[:, 0], linewidth=0.8, color=colors[0], label=labelsLMN[0])
ax1.plot(epochToCal, bLMN[:, 1], linewidth=0.8, color=colors[1], label=labelsLMN[1])
ax1.plot(epochToCal, (-1) * bLMN[:, 2], linewidth=0.8, color=colors[2], label=labelsLMN[2])
ax1.legend(loc=2, fontsize=8, markerscale=1.0, handlelength=2, handleheight=0.5)

xtick1 = [2017, 12, 1, 14, 40, 15, 000]
xtick2 = [2017, 12, 1, 14, 40, 30, 000]
xtick3 = [2017, 12, 1, 14, 40, 45, 000]
xtick4 = [2017, 12, 1, 14, 41, 0, 000]
xtick5 = [2017, 12, 1, 14, 41, 15, 000]
xtick6 = [2017, 12, 1, 14, 41, 30, 000]

xtick1 = cdflib.cdfepoch.compute_tt2000(xtick1)
xtick2 = cdflib.cdfepoch.compute_tt2000(xtick2)
xtick3 = cdflib.cdfepoch.compute_tt2000(xtick3)
xtick4 = cdflib.cdfepoch.compute_tt2000(xtick4)
xtick5 = cdflib.cdfepoch.compute_tt2000(xtick5)
xtick6 = cdflib.cdfepoch.compute_tt2000(xtick6)

epochTick = [xtick1, xtick2, xtick3, xtick4, xtick5, xtick6]
strTick = ['14:40:15', '14:40:30', '14:40:45', '14:41:00', '14:41:15', '14:41:30']

ax0.set_xlim(xlim)
ax0.set_ylim([-35, 35])
ax0.set_xticks(epochTick)
ax0.set_xticklabels(strTick, fontsize=15)
ax1.set_xlim(xlim)
ax1.set_ylim([-35, 35])
ax1.set_xticks(epochTick)
ax1.set_xticklabels(strTick, fontsize=15)

timeSatrtHodo = [2017, 12, 1, 14, 41, 0, 000]
timeStopHodo = [2017, 12, 1, 14, 41, 1, 000]
epochHodoStartType2000 = cdflib.cdfepoch.compute_tt2000(timeSatrtHodo)
epochHodoStopType2000 = cdflib.cdfepoch.compute_tt2000(timeStopHodo)
timeHodoTuple1 = np.where(epochToCal > epochHodoStartType2000)
timeHodoTuple2 = np.where(epochToCal < epochHodoStopType2000)
timeHodoArr1 = timeHodoTuple1[0]
timeHodoArr2 = timeHodoTuple2[0]
timeHodoArrToCal = np.intersect1d(timeHodoArr1, timeHodoArr2)

bN = bLMN[:, 0]
bM = bLMN[:, 1]
bL = bLMN[:, 2]
bNHodo = bL[timeHodoArrToCal]
bMHodo = bM[timeHodoArrToCal]
bLHodo = bN[timeHodoArrToCal]

xlimHodo = [min(timeHodoArrToCal), max(timeHodoArrToCal)]
ax2.plot(timeHodoArrToCal, bNHodo, linewidth=0.8, color=colors[0], label=labelsLMN[0])
ax2.plot(timeHodoArrToCal, bMHodo, linewidth=0.8, color=colors[1], label=labelsLMN[1])
ax2.plot(timeHodoArrToCal, bLHodo, linewidth=0.8, color=colors[2], label=labelsLMN[2])
ax2.legend(loc=2, fontsize=8, markerscale=1.0, handlelength=2, handleheight=0.5)
ax2.set_xlim(xlimHodo)
ax2.set_xticks([])


ax2.set_ylim([-7, 10])


ax3.plot(bLHodo, bMHodo, color='steelblue')
ax3.scatter(bLHodo[0], bMHodo[0], marker='o', color='steelblue')
ax3.scatter(bLHodo[bLHodo.size - 1], bMHodo[bMHodo.size - 1], marker='>', color='steelblue')
ax3.set_xlabel('$B_{L}$')
ax3.set_ylabel('$B_{M}$')
ax4.plot(bLHodo, bNHodo, color='darkorange')
ax4.scatter(bLHodo[0], bNHodo[0], marker='o', linewidth=5, color='darkorange')
ax4.scatter(bLHodo[bLHodo.size - 1], bNHodo[bNHodo.size - 1], marker='>', linewidth=5, color='darkorange')
ax4.set_xlabel('$B_{L}$')
ax4.set_ylabel('$B_{N}$')

#figurePath = 'D:\\wmmFigure\\FBs\\20171201\\twavpol\\'
figurePath = 'D:\\wmmFigure\\FBs\\20171201\\'
if not os.path.exists(figurePath):
    os.makedirs(figurePath)

plt.savefig(figurePath + '\\hodogram.png', dpi=200)
plt.show()

pyEndTime = time.time()
print(
    '==================================================================================================================')
print('Took %s seconds to run.' % (pyEndTime - pyStartTime))