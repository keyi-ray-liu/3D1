import numpy as np
import matplotlib.pyplot as plt
import time
import os
from math import sqrt


# Throughout, NN means nearest neighbor, nNN means next to NN

# init() set up the original coefficients. The initial NN of dopant NN are assumed to be zero. And we start with a small central coeff

def set_para(OS):
    para = {
    'OS'        : OS,
    'orb_count' : 10,
    'NN_count'  : 1,
    'scale'     : 20, 
    'testcase'  : 100,
    'limit'     : 1.5, 
    'gradual'   : 1, 
    'core'      : '4221', 
    'numtype'   : 2,
    }
    
    para['dir'] = 'res/2sublat_' + para['core'] + '_gradual_' + str(para['gradual']) + '_NN_' + str(para['NN_count']) + '_iter_' + str(para['testcase']) + '/'
    if OS == 'Mac':
        location = '/Users/rayliu/Desktop/Code/ortho/'
        os.system('mkdir' + para['dir'])
    else:
        location = 'C:/Users/Ray/Desktop/Code/code1/'
        os.system('powershell.exe mkdir '+ para['dir'])
    ovlp = []
    if para['NN_count'] == 4:
        ovlp.append(location + para['core'] + '/alloverlap1.dat')
        ovlp.append(location + para['core'] + '/alloverlap2.dat')
    elif para['NN_count'] == 1:
        ovlp.append(location + para['core'] + '/simpleoverlap1.dat')
        ovlp.append(location + para['core'] + '/simpleoverlap2.dat')
    para['files'] = ovlp
    
    return para


def init(para):
    # old = np.zeros((para['orb_count'],NN_para['orb_count']))
    old = np.random.rand(para['numtype'], para['orb_count'], para['orb_count'] + para['NN_count'] * para['orb_count']) / para['scale']
    new = np.random.rand(para['numtype'], para['orb_count'], para['orb_count'] + para['NN_count'] * para['orb_count']) / para['scale']
    norm = np.ones((para['numtype'], para['orb_count']))
    for cen in range(10):
        old[0][cen][cen] = 1
        old[1][cen][cen] = 1
    return old, new, norm


def normalize(cen, coeff, ovlp, types, para):
    normalization = 0
    normflag = 1
    typeL = types
    typeR = types
    for siteL in range(1 + para['NN_count']):
        for orbL in range(para['orb_count']):
            for siteR in range(1 + para['NN_count']):
                for orbR in range(para['orb_count']):
                    normalization += innerprod(cen, siteL, siteR, orbL, orbR, ovlp, coeff, normflag, typeL, typeR, para)
    return sqrt(abs(normalization))


# This loads the bare ovlp integral of Sthe 5 NN to all 25 NN and 2nd NN of interest
# Organized as follows: ovlp[0][:][:] is dopant site, and 1-4 are NN. Similarly, ovlp[:][0-4] are ovlp with dopant and NN, and [5-8] are NN of 1st NN, so on.

def readovlp(para):
    # print(np.shape(raw))
    ovlp = np.zeros((2, 1 + para['NN_count'], 1 + para['NN_count'] + para['NN_count'] * para['NN_count'], para['orb_count'], para['orb_count']))
    for ifile in range(para['numtype']):
        raw = np.loadtxt(para['files'][ifile])
        raw = raw[:, 6:]
        for site in range(1 + para['NN_count']):
            for allsite in range(1 + para['NN_count'] + para['NN_count'] * para['NN_count']):
                for lorb in range(para['orb_count']):
                    for rorb in range(para['orb_count']):
                        ovlp[ifile][site][allsite][lorb][rorb] = raw[(1 + para['NN_count'] + para['NN_count'] * para['NN_count']) * site + allsite][ para['orb_count'] * lorb + rorb]
    return ovlp


# setmatrix() sets up the calculation of the all the necessary iNNer products and send them to main iteration loop

def setmatrix(cen, ovlp, coeff, types, para):
    M = np.zeros((para['orb_count'] + para['NN_count'] * para['orb_count'], para['orb_count'] + para['NN_count'] * para['orb_count']))
    normflag = 0
    # update the [0][0] entry
    typeL = types
    typeR = types
    for orbR in range(para['orb_count']):
        for newsite in range(1 + para['NN_count']):
            for neworb in range(para['orb_count']):
                M[orbR][newsite* para['orb_count'] + neworb] = innerprod(cen, newsite, 0, neworb, orbR, ovlp, coeff, normflag, typeL, typeR, para)
    for siteR in range(1, 1 + para['NN_count']):
        typeL = types
        typeR = int(not types)
        for orbR in range(para['orb_count']):
            for newsite in range(1 + para['NN_count']):
                for neworb in range(para['orb_count']):
                    M[siteR * para['orb_count'] + orbR][newsite* para['orb_count'] + neworb] = innerprod(cen, newsite, siteR, neworb, orbR, ovlp, coeff, normflag, typeL, typeR, para)
    return M


def cal_ovlp(cen, coeff, ovlp, types, para):
    res = np.zeros(para['orb_count'] + para['NN_count'] * para['orb_count'])
    normflag = 1
    # update the [0][0] entry
    for orbR in range(para['orb_count']):
        typeL = types
        typeR = types
        res[orbR] = 0
        for newsite in range(1 + para['NN_count']):
            for neworb in range(para['orb_count']):
                res[orbR] += innerprod(cen, newsite, 0, neworb, orbR, ovlp, coeff, normflag, typeL, typeR, para)
    for siteR in range(1, 1 + para['NN_count']):
        typeL = types
        typeR = int(not types)
        for orbR in range(para['orb_count']):
            res[siteR*para['orb_count'] + orbR] = 0
            for newsite in range(1 + para['NN_count']):
                for neworb in range(para['orb_count']):
                    res[siteR*para['orb_count'] + orbR] += innerprod(cen, newsite, siteR, neworb, orbR, ovlp, coeff, normflag, typeL, typeR, para)
    return res


# innerprod calculates the inner product between a singular orbital on the left and the wavefunction on the right. 
# The indices siteL runs from 0-4, siteR runs from 0-4 (cen, NN1, NN2, NN3, NN4)
def innerprod(cen, siteL, siteR, orbL, orbR, ovlp, coeff, normflag, typeL, typeR, para):
    res = 0
    if normflag == 1:
        coeffL = coeff
        coeffR = coeff
    else:
        coeffL = np.ones((para['numtype'], para['orb_count'], para['orb_count'] + para['NN_count'] * para['orb_count']))
        coeffR = coeff
    for allcen in range(para['orb_count'] + para['NN_count'] * para['orb_count']):
        res += coeffL[typeL][cen][siteL * para['orb_count'] + orbL] * coeffR[typeR][orbR][allcen] * ovlp[typeL][siteL][findovlp(siteR, allcen, para)][orbL][findorb(siteR,allcen, para)]
    return res

# utility function that finds the correct overlap index
def findovlp(siteR, allcen, para):
    if siteR == 0:
        index = allcen//para['orb_count']
    elif allcen<para['orb_count']:
        index = siteR
    else:
        index = 1 + para['NN_count'] + (siteR -1) * para['NN_count'] + allcen//para['orb_count'] -1
    return index

# utility that finds the correct orbital index 
def findorb(siteR, allcen, para):
    return allcen%para['orb_count']


# gen_rhs() provides the orthogonality condition, i.e. sets up the (0-ovlp) of the RHS of the matrix equation
def gen_rhs(cen, para):  # here cen indicates the orbital on the central site, which should range from 0 to 9
    rhs = np.zeros(para['orb_count'] + para['NN_count'] * para['orb_count'])
    rhs[cen] = 1
    return rhs


# test function to see if there's repeating pattern in the

def resplot(coeffarray, cen, para):
    ref = range(para['testcase'])
    print(ref)
    fig, ax = plt.subplots(2,1, sharex='all', sharey='all')
    cen1 = coeffarray[0, :, cen, 0]
    cen2 = coeffarray[1, :, cen, 0]
    all1 = coeffarray[0, :, cen, :]
    all2 = coeffarray[1, :, cen, :]
    ax[0].plot(ref, cen1, label='cen', linewidth=10)
    ax[1].plot(ref, cen2, label='cen', linewidth=10)
    ax[0].plot(ref, all1)
    ax[1].plot(ref, all2)
    ax[1].set_xlabel('No. of iterations')
    ax[1].set_ylabel('Orbital Coefficients')
    plt.ylim(-para['limit'], para['limit'])
    # plt.ylim(-0.5,0.0)
    figname = para['dir'] +  'orbital_' + str(cen) 
    plt.savefig(figname)
    plt.show()


def ovlpplot(ovlparray, cen, para):
    ref = range(para['testcase'])
    fig, ax = plt.subplots(2,1, sharex='all', sharey='all')
    ax[0].plot(ref, ovlparray[0, :, cen, :])
    ax[1].plot(ref, ovlparray[1, :, cen, :])
    ax[1].set_xlabel('No. of iterations')
    ax[1].set_ylabel('Wavefunction Overlap')
    plt.ylim(-1.1, 1.1)
    # plt.ylim(-0.5,0.0)
    figname = para['dir'] +  'orbital_' + str(cen) 
    plt.savefig(figname)
    plt.show()


def saveresult(new, norm, seed, para):
    temp = np.zeros((para['numtype'], para['orb_count'], para['orb_count'] * para['NN_count'] + para['orb_count'] + 1))
    txtname = para['dir'] +  'wf'
    seedname = para['dir'] +  'seed'
    for orb in range(para['orb_count']):
        for types in range(para['numtype']):
            temp[types][orb] = np.append(new[types][orb], norm[types][orb])
    np.save( txtname, temp)
    np.save( seedname, seed)

def solve(cen, ovlp, old, types, para):
    LHS = setmatrix(cen, ovlp, old, types, para)
    RHS = gen_rhs(cen, para)
    new = np.zeros((para['numtype'], para['orb_count'], para['orb_count'] + para['NN_count'] * para['orb_count'] ))
    if para['gradual'] == 0:
        new[types][cen] = np.linalg.solve(LHS, RHS)
    else:
        new[types][cen] = 0.8 * old[types][cen] + 0.2 * np.linalg.solve(LHS, RHS)
    norm = normalize(cen, new, ovlp, types, para)
    new = new / norm   
    return new[types][cen], norm

# main() iteratively solve the matrix equation until the old and new coefficients converge.
def main():
    para = set_para('Win')
    ovlp = readovlp(para)
    # print(np.shape(ovlp), ovlp)
    old, new, norm = init(para)
    seed = old
    coeffarray = np.zeros((para['numtype'], para['testcase'], para['orb_count'], para['orb_count'] + para['NN_count'] * para['orb_count']))
    ovlparray = np.zeros((para['numtype'], para['testcase'], para['orb_count'], para['orb_count'] + para['NN_count'] * para['orb_count']))
    start_time = time.time()
    for itest in range(para['testcase']):
        for types in range(para['numtype']):
            for cen in range(para['orb_count']):
                new[types][cen], norm[types][cen] = solve(cen, ovlp, old,  types, para)
                print(norm[types][cen], itest)
            old[types] = new[types]
            for cen in range(para['orb_count']):
                ovlparray[types][itest][cen] = cal_ovlp(cen, old, ovlp, types, para)
            coeffarray[types][itest] = new[types]
    # print(new[5])
    elapsed_time = time.time() - start_time
    print(elapsed_time)
    # set which central function to plot
    cen = 0
    resplot(coeffarray,  cen, para)
    ovlpplot(ovlparray,  cen, para)
    saveresult(new, norm, seed, para)


# np.savetxt('testout', newarray, depara['limit']er='')

main()
