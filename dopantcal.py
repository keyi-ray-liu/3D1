import numpy as np
from math import sqrt
from math import factorial as fact
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits import mplot3d
from matplotlib.widgets import Slider
from matplotlib.widgets import TextBox
from matplotlib.backends.backend_pdf import PdfPages
from multiprocessing import Pool
from functools import partial
from scipy import integrate
import time
import sys
import os

def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path) #if onefile else relative_path

def orbital(n, l, Z, orb, coor, x, y ,z ):

    def normalization():
	    return sqrt((2 / n)** 3 * fact(n - l -1) / (2 * n * fact(n + l)))

    pi = np.pi
    r = [x, y, z]
    vec = np.subtract(r, coor)* 5.4307
    x, y, z = vec
    rr = [vec.dot(vec) if vec.any() else 1][0]

    if orb == 0:
        sph = 0.5 * sqrt(1/pi)
    elif orb == 1:
        sph = sqrt(0.75 / pi) * x / rr
    elif orb == 2:
        sph = sqrt(0.75 / pi) * y / rr
    elif orb == 3:
        sph = sqrt(0.75 / pi) * z / rr
    elif orb == 4:
        sph = 0.5 * sqrt(15 / pi) * x * y / rr ** 2
    elif orb == 5:
        sph = 0.5 * sqrt(15 / pi) * y * z / rr ** 2
    elif orb == 6:
        sph = 0.5 * sqrt(15 / pi) * x * z / rr ** 2
    elif orb == 7:
        sph = 0.25 * sqrt(15 / pi) * (x**2 - y**2) / rr ** 2
    elif orb == 8:
        sph = 0.25 * sqrt(5 / pi) * (2 * z**2 - x**2 - y**2) / rr ** 2
    elif orb == 9:
        sph = 0.5 * sqrt(1/pi)
    else:
        print('wrong orb!')
        

    alpha = Z * vec.dot(vec) / n
    return normalization() * Z ** 1.5 * laguerre(n - l - 1, 2 * l + 1, 2 * alpha) * np.exp(-alpha) * (2 * alpha) ** l * sph

def laguerre(n, k, x):
	res = 0
	for m in range(n + 1):
		res += (-1) ** m * fact(n + k) / (fact(n - m) * fact(k + m) * fact(m)) * x ** m
	return res

def wf(site, orb, numsite, numorb, n, l, Z, coors, x, y, z, coeff):
    temp = 0
    for subsite in range(numsite):
        for suborb in range(numorb):
            temp += coeff[site * numorb + orb][subsite * numorb + suborb] * orbital(n[suborb], l[suborb], Z[suborb], suborb, coors[subsite], x, y ,z)
    return temp

def setorb():
    n = [3] * 9 + [4]
    l = [0] + [1] * 3 + [2] *5 + [0]
    Z = [4.5] + [1.4] * 3 + [4.4] * 5 + [1.5] * 1
    return n, l, Z

def process(coeff, coors, den):
    numsite, numorb = 17, 10
    rup, rdown, iup, idown = 0, 0, 0, 0
    n, l, Z = setorb()
    for site in range(numsite):
        for orb in range(numorb):
            temp = wf(site, orb, numsite, numorb, n, l, Z, coors, 0.0, 0.0, 0.0, coeff)
            rup += temp * den[site][orb*2]
            iup += temp * den[site][orb*2 +1]
            rdown += temp * den[site][orb*2 + 20]
            idown += temp * den[site][orb*2 + 21]
            #print(temp, rup, iup, rdown, idown)
    res = rup ** 2 + rdown ** 2 + iup**2 + idown** 2
    print(res)

def readcoeffs():
    numsite = 17
    coeff = np.loadtxt(resource_path('coeff'))
    coors = np.loadtxt(resource_path('listpos'))
    coors *= 5.43
    den = np.loadtxt(resource_path('density'))
    den = np.reshape(den, (len(den)//numsite, numsite, len(den[0])))
    den = den[:, :, 4:]
    print(den.shape)
    return coeff, coors, den


def writewf(numsite , numorb, plotBoth, coors, coeff, sitedict, orbdict):
    n, l, Z = setorb()
    i = 0
    for site in range(numsite):
        x, y, z= coors[site]
        print('Printing statistics on site {}: \n\n'.format(sitedict[site]))
        for subsite in range(numsite):
            if subsite == site:
                print('-------------------')
            for suborb in range(numorb):
                i+=1
                ortho = wf(subsite, suborb, numsite, numorb, n, l, Z, coors, x, y, z, coeff)
                ori = orbital(n[suborb], l[suborb], Z[suborb], suborb, coors[subsite], x, y, z)
                if subsite==site :
                    print("{} **** Contribution from wavefunction of site {}, orb {} ------> site {} is {:.8f}, original is {:.8f} ".format(i, sitedict[subsite], orbdict[suborb], sitedict[site], ortho, ori))
                else:
                    print("{}   Contribution from wavefunction of site {}, orb {} ------> site {} is {:.8f}, original is {:.8f}".format(i, sitedict[subsite], orbdict[suborb], sitedict[site], ortho, ori))
            if subsite == site:
                print('-------------------')
        print('\n\n')


def plotwf(numsite, numorb, plotBoth, coors, coeff):
    n, l, Z = setorb()
    rs = np.linspace(-1.5, 1.5, 60)

    ctrlAx = [plt.axes([0.25, 0.1 ,0.3, 0.15]), plt.axes([0.25, 0.3 ,0.3, 0.15]), plt.axes([0.25, 0.5 ,0.3, 0.15])]
    txtz = TextBox(ctrlAx[0], 'z=', initial=0.0)
    txtorb = TextBox(ctrlAx[1], 'Orb (1-10):', initial=1)
    txtsite = TextBox(ctrlAx[2], 'Site (1-17):', initial=1)
    fig, ax2 = plt.subplots(subplot_kw={"projection":"3d"})
    X = np.array([[x for y in rs] for x in rs]).flatten()
    Y = np.array([[y for y in rs] for x in rs]).flatten()
    
    def submit(val):
        z = float(txtz.text)
        orb = int(txtorb.text) - 1
        site = int(txtsite.text) - 1
        ax2.clear()
        if plotBoth:
            wfs = np.array([[ wf(site, orb, numsite, numorb, n, l, Z, coors, x, y, z, coeff) for y in rs] for x in rs]).flatten()
            ax2.scatter(X, Y, wfs, marker='o' )
        ori = np.array([[ orbital(n[orb], l[orb], Z[orb], orb, coors[site], x, y, z) for y in rs] for x in rs]).flatten()
        ax2.scatter(X, Y, ori, marker='.' )
        
        #plt.legend()
        plt.title('site{}, orb{}, z={}, num of site={}, num of orb={}'.format(site, orb, z, numsite, numorb))
        fig.canvas.draw_idle()
    
    txtz.on_submit(submit)
    txtorb.on_submit(submit)
    txtsite.on_submit(submit)
    plt.show()

def plotprewf(numsite, numorb, plotBoth, coors, coeff, wfs, oris, grid):
    n, l, Z = setorb()
    rs = np.linspace(-2.8, 2.8, grid)

    zmax = wfs.max() + 0.1
    zmin = wfs.min() - 0.1

    auxAx = [[plt.axes([0.25, 0.1 ,0.3, 0.09]), plt.axes([0.25, 0.2 ,0.3, 0.09]), plt.axes([0.25, 0.3 ,0.3, 0.09])], \
        [plt.axes([0.25, 0.6 ,0.3, 0.09]), plt.axes([0.25, 0.7 ,0.3, 0.09]), plt.axes([0.25, 0.8 ,0.3, 0.09])]]
    sliderzL = Slider(auxAx[0][0], 'WF1 z: ( lat. const.)', -0.5, 0.5, valstep=0.25)
    txtorbL = TextBox(auxAx[0][1], 'WF1 Orb (1-10):', initial=1)
    txtsiteL = TextBox(auxAx[0][2], 'WF1 Site (1-17):', initial=1)

    sliderzR = Slider(auxAx[1][0], 'WF2 z: ( lat. const.)', -0.5, 0.5, valstep=0.25)
    txtorbR = TextBox(auxAx[1][1], 'WF2 Orb (1-10):', initial=1)
    txtsiteR = TextBox(auxAx[1][2], 'WF2 Site (1-17):', initial=1)

    fig1, ax1 = plt.subplots(subplot_kw={"projection":"3d"})
    fig2, ax2 = plt.subplots(subplot_kw={"projection":"3d"})

    X = np.array([[x for y in rs] for x in rs]).flatten()
    Y = np.array([[y for y in rs] for x in rs]).flatten()
    
    def submitL(val):
        z = sliderzL.val * 4
        orb = int(txtorbL.text) - 1
        site = int(txtsiteL.text) - 1
        ax1.clear()
        if plotBoth:
            wf = wfs [(site * 10 + orb) * 5 + int(z + 2.0)]
            ax1.scatter(X, Y, wf, marker='o', s=0.5 )
        ori = oris [(site * 10 + orb) * 5 + int(z + 2.0)]
        ax1.scatter(X, Y, ori, marker='.', s=0.5 )
        ax1.set_zlim(zmin, zmax)
        #plt.legend()
        ax1.set_title('site{}, orb{}, z={}'.format(site+1, orb+1, z))
        fig1.canvas.draw_idle()

    def submitR(val):
        z = sliderzR.val * 4
        orb = int(txtorbR.text) - 1
        site = int(txtsiteR.text) - 1
        ax2.clear()
        if plotBoth:
            wf = wfs [(site * 10 + orb) * 5 + int(z + 2.0)]
            ax2.scatter(X, Y, wf, marker='o', s=0.5 )
        ori = oris [(site * 10 + orb) * 5 + int(z + 2.0)]
        ax2.scatter(X, Y, ori, marker='.', s=0.5 )
        ax2.set_zlim(zmin, zmax)
        #plt.legend()
        ax2.set_title('site{}, orb{}, z={}'.format(site+1, orb+1, z))
        fig2.canvas.draw_idle()
    
    sliderzL.on_changed(submitL)
    txtorbL.on_submit(submitL)
    txtsiteL.on_submit(submitL)

    sliderzR.on_changed(submitR)
    txtorbR.on_submit(submitR)
    txtsiteR.on_submit(submitR)

    plt.show()

def pdfplot(numsite, numorb, coors, coeff, wfs, oris):
    n, l, Z = setorb()
    rs = np.linspace(-4.0, 4.0, 100)
    X = np.array([[x for y in rs] for x in rs]).flatten()
    Y = np.array([[y for y in rs] for x in rs]).flatten()
    zmax = wfs.max() + 0.1
    zmin = wfs.min() - 0.1
    start = time.time()

    with PdfPages('allplot.pdf') as pdf:
        for site in range(numsite):
            for orb in range(numorb):
                for z in range(5):
                    fig, ax2 = plt.subplots(subplot_kw={"projection":"3d"})
                    wf = wfs [(site * 10 + orb) * 5 + z]
                    ori = oris [(site * 10 + orb) * 5 + z]
                    ax2.scatter(X, Y, wf, marker='o', s=0.1 )   
                    ax2.scatter(X, Y, ori, marker='.', s=0.1 )
                    ax2.set_zlim(zmin, zmax)
                    ax2.set_title('site {}, orb {}, z={}'.format(site, orb, float(z-2)))
                    print('site = {}, orb= {}, z={} plot complete'.format(site, orb, z))
                    pdf.savefig()
                    plt.close()

    end = time.time()
    print('time elapsed: {}'.format(end - start))


def prepwf(plotsite, plotorb, coors, coeff, grid):
    pool = Pool()

        
    res = pool.map(partial(prepsinglewf, plotsite=plotsite, plotorb=plotorb, coors=coors, coeff=coeff, grid=grid), range(plotsite))
    wfs = np.concatenate([wf[0] for wf in res])
    oris = np.concatenate([wf[1] for wf in res])

    np.savetxt('wfs', wfs)
    np.savetxt('oris', oris)

    
def prepsinglewf(site, plotsite, plotorb, coors, coeff, grid):
    n, l, Z = setorb()
    numsite, numorb = 17, 10
    rs = np.linspace(-2.8, 2.8, grid)
    zs = np.linspace(-2.0, 2.0, 5) * 0.25 * 5.4307
    wfs = [0] * (plotorb * len(zs))
    oris = [0] * (plotorb * len(zs))
    for orb in range(plotorb):
        for i, z in enumerate(zs):
            wfs [orb * len(zs) + i] = np.array([[ wf(site, orb, numsite, numorb, n, l, Z, coors, x, y, z, coeff) for y in rs] for x in rs]).flatten()
            oris [orb * len(zs) + i] = np.array([[ orbital(n[orb], l[orb], Z[orb], orb, coors[site], x, y, z) for y in rs] for x in rs]).flatten()
            print('site = {}, orb= {}, z={} calc complete'.format(site, orb, z))
    return wfs, oris

def integraltest(numsite, numorb, coors, coeff):
    n, l, Z = setorb()
    integrallist = setlist()
    for site1, site2, orb1, orb2 in integrallist:
        print('start')
        def func(x, y, z):
            return wf(site1, orb1, numsite, numorb, n, l, Z, coors, x, y ,z , coeff) * wf(site2, orb2, numsite, numorb, n, l, Z, coors, x, y, z, coeff)

        print(integrate.tplquad(func, -5, 5, lambda x: -5, lambda x:5, lambda x, y: -5, lambda x, y: 5))

def setlist():
    res  = [[0, 0, 0, 0]]
    res += [[0, 0, 0, 1]]
    res += [[0, 0, 1, 0]]
    res += [[0, 0, 1, 1]]
    return res

def setdict():
    orbdict = {0: 's', 1: 'px', 2:'py', 3:'pz', 4:'dxy', 5:'dyz', 6:'dxz', 7:'dx2y2', 8:'dz2', 9:'s*'}
    sitedict = {0: 'dop', 1: 'NN1', 2:'NN2', 3:'NN3', 4:'NN4'}
    for i in range(12):
        sitedict[i+5] = 'nNN' + str(i+1)
    return sitedict, orbdict

def sumden(den):
    return den[0][0] ** 2+ den[0][1]** 2 + den[0][20]**2 + den[0][21]**2
    
if __name__ == '__main__':
    coeff, coors, dens = readcoeffs()
    sitedict, orbdict = setdict()
    test = 0
    if test==1:
        np.savetxt('pden', dens[0])
        process(coeff, coors, dens[0])
    elif test==2:
        for den in dens:
            print(sumden(den))
            #process(coeff, coors, den)
    site, orb, plotBoth, pre, calc, grid = np.loadtxt(resource_path('paras'), dtype=int)
    if not pre:
        #plotwf(site , orb, plotBoth, coors, coeff)
        #writewf(site , orb, plotBoth, coors, coeff, sitedict, orbdict)
        integraltest(site, orb, coors, coeff)
    elif not calc:
        wfs = np.loadtxt(resource_path('wfs'))
        oris = np.loadtxt(resource_path('oris'))
        plotprewf(site , orb, plotBoth, coors, coeff, wfs, oris, grid)
        #pdfplot(site, orb, coors, coeff, wfs, oris)
    else:
        start = time.time()
        prepwf(17, 10, coors, coeff, grid)
        end = time.time()
        print(end - start)
