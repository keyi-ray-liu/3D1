import numpy as np
import matplotlib.pyplot as plt
import time
from math import sqrt


# Throughout, NN means nearest neighbor, nNN means next to NN

# init() set up the original coefficients. The initial NN of dopant NN are assumed to be zero. And we start with a small central coeff

def set_para():
    global orb_count, NN_count, scale, testcase, limit, file1, file2, gradual, para, numtype
    orb_count, NN_count, scale, testcase, limit, gradual, para, numtype = 10, 1, 20, 2, 1.5, 1, "4221", 2
    if NN_count == 4:
        file1 = '/Users/rayliu/Desktop/Code/ortho/' + para + '/alloverlap1.dat'
        file2 = '/Users/rayliu/Desktop/Code/ortho/' + para + '/alloverlap2.dat'
    elif NN_count == 1:
        file1 = '/Users/rayliu/Desktop/Code/ortho/' + para + '/simpleoverlap1.dat'
        file2 = '/Users/rayliu/Desktop/Code/ortho/' + para + '/simpleoverlap2.dat'


def init():
    # old = np.zeros((orb_count,NN_orb_count))
    old = np.random.rand(numtype, orb_count, orb_count + NN_count * orb_count) / scale
    new = np.random.rand(numtype, orb_count, orb_count + NN_count * orb_count) / scale
    norm = np.ones((numtype, orb_count))
    for cen in range(10):
        old[0][cen][cen] = 1
        old[1][cen][cen] = 1
    return old, new, norm


def normalize(cen, coeff, ovlp, types):
    normalization = 0
    normflag = 1
    typeL = types
    typeR = types
    for siteL in range(1 + NN_count):
        for orbL in range(orb_count):
            for siteR in range(1 + NN_count):
                for orbR in range(orb_count):
                    normalization += innerprod(cen, siteL, siteR, orbL, orbR, ovlp, coeff, normflag, typeL, typeR)
    return sqrt(abs(normalization))


# This loads the bare ovlp integral of Sthe 5 NN to all 25 NN and 2nd NN of interest
# Organized as follows: ovlp[0][:][:] is dopant site, and 1-4 are NN. Similarly, ovlp[:][0-4] are ovlp with dopant and NN, and [5-8] are NN of 1st NN, so on.

def readovlp(file):
    raw = np.loadtxt(file)
    raw = raw[:, 6:]
    # print(np.shape(raw))
    ovlp = np.zeros((1 + NN_count, 1 + NN_count + NN_count * NN_count, orb_count, orb_count))
    for site in range(1 + NN_count):
        for allsite in range(1 + NN_count + NN_count * NN_count):
            for lorb in range(orb_count):
                for rorb in range(orb_count):
                    ovlp[site][allsite][lorb][rorb] = raw[(1 + NN_count + NN_count * NN_count) * site + allsite][
                        orb_count * lorb + rorb]
    return ovlp


# setmatrix() sets up the calculation of the all the necessary iNNer products and send them to main iteration loop

def setmatrix(cen, ovlp, coeff, types):
    M = np.zeros((orb_count + NN_count * orb_count, orb_count + NN_count * orb_count))
    normflag = 0
    # update the [0][0] entry
    typeL = types
    typeR = types
    for orbR in range(orb_count):
        for newsite in range(1 + NN_count):
            for neworb in range(orb_count):
                M[orbR][newsite* orb_count + neworb] = innerprod(cen, newsite, 0, neworb, orbR, ovlp, coeff, normflag, typeL, typeR)
    for siteR in range(1, 1 + NN_count):
        typeL = types
        typeR = int(not types)
        for orbR in range(orb_count):
            for newsite in range(1 + NN_count):
                for neworb in range(orb_count):
                    M[siteR * orb_count + orbR][newsite* orb_count + neworb] = innerprod(cen, newsite, siteR, neworb, orbR, ovlp, coeff, normflag, typeL, typeR)
    return M


def cal_ovlp(cen, coeff, ovlp, types):
    res = np.zeros(orb_count + NN_count * orb_count)
    normflag = 1
    # update the [0][0] entry
    for orbR in range(orb_count):
        typeL = types
        typeR = types
        res[orbR] = 0
        for newsite in range(1 + NN_count):
            for neworb in range(orb_count):
                res[orbR] += innerprod(cen, newsite, 0, neworb, orbR, ovlp, coeff, normflag, typeL, typeR)
    for siteR in range(1, 1 + NN_count):
        typeL = types
        typeR = int(not types)
        for orbR in range(orb_count):
            res[siteR*orb_count + orbR] = 0
            for newsite in range(1 + NN_count):
                for neworb in range(orb_count):
                    res[siteR*orb_count + orbR] += innerprod(cen, newsite, siteR, neworb, orbR, ovlp, coeff, normflag, typeL, typeR)
    return res


# innerprod calculates the inner product between a singular orbital on the left and the wavefunction on the right. 
# The indices siteL runs from 0-4, siteR runs from 0-4 (cen, NN1, NN2, NN3, NN4)
def innerprod(cen, siteL, siteR, orbL, orbR, ovlp, coeff, normflag, typeL, typeR):
    res = 0
    if normflag == 1:
        coeffL = coeff
        coeffR = coeff
    else:
        coeffL = np.ones((numtype, orb_count, orb_count + NN_count * orb_count))
        coeffR = coeff
    for allcen in range(orb_count + NN_count * orb_count):
        res += coeffL[typeL][cen][siteL * orb_count + orbL] * coeffR[typeR][orbR][allcen] * ovlp[typeL][siteL][findovlp(siteR, allcen)][orbL][findorb(siteR,allcen)]
    return res

# utility function that finds the correct overlap index
def findovlp(siteR, allcen):
    if siteR == 0:
        index = allcen//orb_count
    elif allcen<orb_count:
        index = siteR
    else:
        index = 1 + NN_count + (siteR -1) * NN_count + allcen//orb_count -1
    return index

# utility that finds the correct orbital index 
def findorb(siteR, allcen):
    return allcen%orb_count


# gen_rhs() provides the orthogonality condition, i.e. sets up the (0-ovlp) of the RHS of the matrix equation
def gen_rhs(cen):  # here cen indicates the orbital on the central site, which should range from 0 to 9
    rhs = np.zeros(orb_count + NN_count * orb_count)
    rhs[cen] = 1
    return rhs


# test function to see if there's repeating pattern in the

def resplot(coeffarray, testcase, cen):
    ref = range(testcase)
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
    plt.ylim(-limit, limit)
    # plt.ylim(-0.5,0.0)
    figname = './res/2para_' + para + '_orbital_' + str(cen) + '_gradual_' + str(gradual) + '_NN_' + str(
        NN_count) + '_iter_' + str(testcase)
    plt.savefig(figname)
    plt.show()


def ovlpplot(ovlparray, testcase, cen):
    ref = range(testcase)
    fig, ax = plt.subplots(2,1, sharex='all', sharey='all')
    ax[0].plot(ref, ovlparray[0, :, cen, :])
    ax[1].plot(ref, ovlparray[1, :, cen, :])
    ax[1].set_xlabel('No. of iterations')
    ax[1].set_ylabel('Wavefunction Overlap')
    plt.ylim(-1.1, 1.1)
    # plt.ylim(-0.5,0.0)
    figname = './res/2para_' + para + '_orbital_' + str(cen) + '_gradual_' + str(gradual) + '_NN_' + str(
        NN_count) + '_iter_' + str(testcase) + 'ovlp'
    plt.savefig(figname)
    plt.show()


def saveresult(new, norm, seed):
    temp = np.zeros((numtype, orb_count, orb_count * NN_count + orb_count + 1))
    txtname = '2para_' + para + '_gradual_' + str(gradual) + '_NN_' + str(NN_count) + '_iter_' + str(testcase) + '.dat'
    for orb in range(orb_count):
        for types in range(numtype):
            temp[types][orb] = np.append(new[types][orb], norm[types][orb])
    np.save('./res/' + txtname, temp)
    np.save('./res/seed' + txtname, seed)

def solve(cen, ovlp, old, gradual, types):
    LHS = setmatrix(cen, ovlp, old, types)
    RHS = gen_rhs(cen)
    new = np.zeros((numtype, orb_count, orb_count + NN_count * orb_count ))
    if gradual == 0:
        new[types][cen] = np.linalg.solve(LHS, RHS)
    else:
        new[types][cen] = 0.8 * old[types][cen] + 0.2 * np.linalg.solve(LHS, RHS)
    norm = normalize(cen, new, ovlp, types)
    new = new / norm   
    return new[types][cen], norm

# main() iteratively solve the matrix equation until the old and new coefficients converge.
def main():
    set_para()
    ovlp = np.zeros((numtype, 1 + NN_count, 1 + NN_count + NN_count*NN_count , orb_count, orb_count))
    ovlp[0], ovlp[1] = readovlp(file1), readovlp(file2)
    # print(np.shape(ovlp), ovlp)
    old, new, norm = init()
    seed = old
    coeffarray = np.zeros((2, testcase, orb_count, orb_count + NN_count * orb_count))
    ovlparray = np.zeros((2, testcase, orb_count, orb_count + NN_count * orb_count))
    start_time = time.time()
    for itest in range(testcase):
        for types in range(numtype):
            for cen in range(orb_count):
                new[types][cen], norm[types][cen] = solve(cen, ovlp, old, gradual, types)
                print(norm[types][cen], itest)
            old[types] = new[types]
            for cen in range(orb_count):
                ovlparray[types][itest][cen] = cal_ovlp(cen, old, ovlp, types)
            coeffarray[types][itest] = new[types]
    # print(new[5])
    elapsed_time = time.time() - start_time
    print(elapsed_time)
    # set which central function to plot
    cen = 0
    resplot(coeffarray, testcase, cen)
    ovlpplot(ovlparray, testcase, cen)
    saveresult(new, norm, seed)


# np.savetxt('testout', newarray, delimiter='')

main()
