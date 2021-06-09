import enum
import numpy as np
import sys
from functools import reduce


def readelectron(tag):
    raw = np.loadtxt('electrons.{}'.format(tag))[:, 4:]
    return [[[ row[i*2 ] + 1j* row[i*2+1] for i in range(len(row)//2)] for row in raw ]]

def readoffdiag():
    raw =  np.loadtxt('offdiag')
    t = np.zeros((len(raw), 20, 20))
    c = 0
    for i, row in enumerate(raw):
        for j in range(10):
            for k in range(j, 10):
                t[i][j][k] = t[i][k][j] = t[i][j+10][k+10] = t[i][k+10][j+10] = row[c]
                c += 1
    return t


def calpert(wfs, t):
    return sum([reduce(np.multiply, np.conj(wf), t[i], wf) for i, wf in enumerate(wfs)])

if __name__ =='__main__':
    for tag in range(1, 13):
        wfs, t = readelectron(tag), readoffdiag()
        print(calpert(wfs, t))