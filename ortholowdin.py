import numpy as np
import matplotlib.pyplot as plt
import time
import os
import copy
from matplotlib.backends.backend_pdf import PdfPages
from math import sqrt


# Throughout, NN means nearest neighbor, nNN means next to NN

# init() set up the original coefficients. The initial NN of dopant NN are assumed to be zero. And we start with a small central coeff

def set_para(OS, NN, resdir):
    para = {
    'OS'        : OS,
    'orb_count' : 10,
    'NN_count'  : NN,
    'scale'     : 20, 
    'testcase'  : 300,
    'limit'     : 1.5, 
    'gradual'   : 0.05, 
    'core'      : '4221', 
    'numtype'   : 2,
    'sym_op'    : 0,
    'lower'     : 80,
    'upper'     : 150
    }

    if para['sym_op'] == 1:
        para['dir'] = 'res/highpre2sublat' +  '_sym_' +  str(para['sym_op']) + '_core_' + para['core'] + '_gradual_' + str(para['gradual']) + '_NN_1_iter_' + str(para['testcase']) + '/'
    else:
        para['dir'] = 'res/highpre2sublat' +  '_sym_' +  str(para['sym_op']) + '_core_' + para['core'] + '_gradual_' + str(para['gradual']) + '_NN_' + str(para['NN_count']) + '_iter_' + str(para['testcase']) + '/'
    
    if OS == 'Mac':
        location = '/Users/rayliu/Desktop/Code/ortho/'
        if resdir:
            os.system('mkdir' + para['dir'])
    else:
        location = 'C:/Users/Ray/Desktop/Code/code1/'
        if resdir:
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



def init_res_array(para):
    coeffarray = np.zeros((para['numtype'], para['testcase'], para['orb_count'], para['orb_count'] + para['NN_count'] * para['orb_count']))
    ovlparray = np.zeros((para['numtype'], para['testcase'], para['orb_count'], para['orb_count'] + para['NN_count'] * para['orb_count']))
    return coeffarray, ovlparray

def init_sym(para):
    sym = np.zeros((para['numtype'], para['orb_count'], para['orb_count'] + para['NN_count'] * para['orb_count']))
    symarray = pmarray()
    symovlparray = np.zeros((para['numtype'], para['testcase'], para['orb_count'], para['orb_count'] + para['NN_count'] * para['orb_count']))
    return sym, symarray, symovlparray

def init_data(para):
    wf = np.load(para['dir'] + 'wf.npy')
    ovlp = np.load(para['dir'] + 'ovlp.npy')
    sym = np.load(para['dir'] + 'symovlp.npy')
    return wf, ovlp, sym

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

def solve(cen, ovlp, old, types, para):
    LHS = setmatrix(cen, ovlp, old, types, para)
    RHS = gen_rhs(cen, para)
    new = np.zeros((para['numtype'], para['orb_count'], para['orb_count'] + para['NN_count'] * para['orb_count'] ))
    new[types][cen] = (1 - para['gradual']) * old[types][cen] + para['gradual'] * np.linalg.solve(LHS, RHS)
    norm = normalize(cen, new, ovlp, types, para)
    new = new / norm   
    return new[types][cen], norm

# symmetry takes the coefficients of the 1-NN calculations and propogates it to 4-NN according to symmetry. 
def symmetry(old, types, symarray, para, sympara):
    res = np.zeros((sympara['orb_count'], sympara['orb_count'] + sympara['NN_count'] * sympara['orb_count']))
    for cen in range(sympara['orb_count']):
        res[cen][:para['orb_count'] + para['NN_count'] * para['orb_count']] = old[types][cen][:]
        for restNN in range(sympara['NN_count'] - para['NN_count']):
            for orb in range(para['orb_count']):
                res[cen][para['orb_count'] * (2 + restNN) + orb] = old[types][cen][orb] * symarray[types][restNN * para['orb_count'] + orb]
    return res


def resplot(coeffarray, para):
    ref = range(para['testcase'])
    figname = para['dir'] +  'allorbital_'  + '.pdf'
    with PdfPages(figname) as pdf:
        for cen in range(para['orb_count']):
            plt.figure(figsize=(3, 3))
            fig, ax = plt.subplots(2,1, sharex='all', sharey='all')
            fig.suptitle('Central Orbital {}' .format(cen+1))
            cen1 = coeffarray[0, :, cen, cen]
            cen2 = coeffarray[1, :, cen, cen]
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
            pdf.savefig()
            plt.close()
    

def ovlpplot(ovlparray, para):
    ref = range(para['upper'] - para['lower'])
    figname = para['dir'] +  'allorbital_' + '_NN_' + str(para['NN_count']) +'_overlap' + '_zoom_' + str(para['lower']) + str(para['upper']) +  '.pdf' 
    with PdfPages(figname) as pdf:
        for cen in range(para['orb_count']):
            plt.figure(figsize=(3, 3))
            fig, ax = plt.subplots(2,1, sharex='all', sharey='all')
            fig.suptitle("Central Orbital Overlap {}" .format(cen+1))
            ax[0].plot(ref, ovlparray[0, para['lower']:para['upper'], cen, :])
            ax[1].plot(ref, ovlparray[1, para['lower']:para['upper'], cen, :])
            ax[1].set_xlabel('No. of iterations')
            ax[1].set_ylabel('Wavefunction Overlap')
            plt.ylim(-0.25, 0.25)
    # plt.ylim(-0.5,0.0)
            pdf.savefig()
            plt.close()


def saveresult(newarray, ovlparray, symovlparray, norm, seed, para):
    wfname = para['dir'] +  'wf'
    ovlpname = para['dir'] + 'ovlp'
    symname = para['dir'] + 'symovlp'
    seedname = para['dir'] +  'seed'
    normname = para['dir'] +  'norm'
    np.save( wfname, newarray)
    np.save( ovlpname, ovlparray)
    np.save( symname, symovlparray)
    np.save( seedname, seed)
    np.save( normname, norm)


def printresult(coeffarray, ovlparray, para):
    for itr, ovlp in enumerate(np.sum(abs(ovlparray[0, para['lower']:para['upper'], 0, 1:]), axis=1)/49):
        print(itr, ovlp)
    #np.savetxt(para['dir'] + 'coeff1.dat', np.reshape(coeffarray[0, 90, :, :], 500)[None], delimiter='  ')
    #np.savetxt(para['dir'] + 'coeff2.dat', np.reshape(coeffarray[1, 90, :, :], 500)[None], delimiter='  ')
    np.savetxt(para['dir'] + 'coeff1.dat', coeffarray[0, 90, :, :], delimiter='  ')
    np.savetxt(para['dir'] + 'coeff2.dat', coeffarray[0, 90, :, :], delimiter='  ')


# main() iteratively solve the matrix equation until the old and new coefficients converge.
def main(PostProcess, Repeat):

    
    para= set_para('Win', 4, 1)
    sympara = set_para('Win', 4, 0)
    ovlp = readovlp(para)
    symovlp = readovlp(sympara)

    if not PostProcess:
    # print(np.shape(ovlp), ovlp)
        seed, old, new, norm = init_coeff(para)

        if Repeat:
            old = np.load(para['dir'] + 'seed.npy')

        coeffarray, ovlparray = init_res_array(para)
    

        if sympara['sym_op'] == 1:
            sym, symarray, symovlparray = init_sym(sympara)
        else:
            symovlparray = np.zeros(1)

        start_time = time.time()
        for itest in range(para['testcase']):
            for types in range(para['numtype']):
                for cen in range(para['orb_count']):
                    new[types][cen], norm[types][cen] = solve(cen, ovlp, old,  types, para)
                    print(norm[types][cen], itest)
                old[types] = new[types]

                if sympara['sym_op'] == 1:
                    sym[types] = symmetry(old, types, symarray, para, sympara)

                for cen in range(para['orb_count']):
                    ovlparray[types][itest][cen] = cal_ovlp(cen, old, ovlp, types, para)

                    if sympara['sym_op'] == 1:
                        symovlparray[types][itest][cen] = cal_ovlp(cen, sym, symovlp, types, sympara)

                coeffarray[types][itest] = new[types]
    # print(new[5])
        elapsed_time = time.time() - start_time
        print(elapsed_time)
    # set which central function to plot
        saveresult(coeffarray, ovlparray, symovlparray, norm, seed, para)

    else:
        coeffarray, ovlparray, symovlparray = init_data(para)
        #resplot(coeffarray,  para)
        #ovlpplot(ovlparray,  para)
        printresult(coeffarray, ovlparray, para)
        if sympara['sym_op'] == 1:
            ovlpplot(symovlparray, sympara)
    
# np.savetxt('testout', newarray, depara['limit']er='')

main(1, 0)
