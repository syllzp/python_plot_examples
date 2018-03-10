#!/usr/bin/env python

import os
import re
import numpy as np
from optparse import OptionParser

import matplotlib as mpl
from matplotlib.ticker import AutoMinorLocator
mpl.use('agg')
import matplotlib.pyplot as plt

import matplotlib.colors as mcolors
from matplotlib.patches import Polygon

mpl.rcParams['axes.unicode_minus'] = False

############################################################
__version__ = "1.0"
############################################################

def nelect_from_OUTCAR(infile='OUTCAR'):
    assert os.path.isfile(infile), '%s cannot be found!' % infile
    nele = [float(line.split()[2]) for line in open(infile) if 'NELECT' in line]
    return nele

def vacumm_from_OUTCAR(infile='OUTCAR'):
    assert os.path.isfile(infile), '%s cannot be found!' % infile
    vacumm = np.asarray([line.split()[-2:] for line in open(infile) if 'vacuum' in line],dtype=float)
    return vacumm

def Efermi_from_OUTCAR(infile='OUTCAR'):
    assert os.path.isfile(infile), '%s cannot be found!' % infile
    Ef = [float(line.split()[2]) for line in open(infile) if 'E-fermi' in line]
    return Ef

def WeightFromPro(infile='PROCAR'):
    """
    Contribution of selected atoms to the each KS orbital
    """

    assert os.path.isfile(infile), '%s cannot be found!' % infile
    FileContents = [line for line in open(infile) if line.strip()]

    # when the band number is too large, there will be no space between ";" and
    # the actual band number. A bug found by Homlee Guo.
    # Here, #kpts, #bands and #ions are all integers
    nkpts, nbands, nions = [int(xx) for xx in re.sub('[^0-9]', ' ', FileContents[1]).split()]
    nspin = len([line for line in FileContents if 'bands' in line])

    #Weights = np.asarray([line.split()[1:-1] for line in FileContents
    #                     if not re.search('[a-zA-Z]', line)], dtype=float)

    kpt_weight = np.asarray([line.split()[-1] for line in FileContents if 'weight' in line], dtype=float)
    
    energies = np.asarray([line.split()[-4] for line in FileContents if 'occ.' in line], dtype=float)

    #nlmax = Weights.shape[-1]
    #nspin = Weights.shape[0] / (nkpts * nbands * nions)
    #nspin /= 4 if lsorbit else 1

    #if lsorbit:
    #    Weights.resize(nspin, nkpts, nbands, 4, nions, nlmax)
    #    Weights = Weights[:,:,:,0,:,:]
    #else:
    #   Weights.resize(nspin, nkpts, nbands, nions, nlmax)

    kpt_weight.resize(nspin, nkpts)
    energies.resize(nspin, nkpts, nbands)
    
    return energies, kpt_weight

############################################################
def calc_Efermi_at_0K(infile='PROCAR'):
    ens, kptw = WeightFromPro(infile)
    nele = nelect_from_OUTCAR()
    aen=ens.copy(); aen.resize(aen.size);aen=np.sort(aen);
    aen[:] = aen[::-1]
    index=np.where(aen<0)
    i=0
    en = aen[index[0][i]] 
    occ=(np.multiply((ens<=en).sum(axis=2),kptw)).sum().sum()
    while(abs(occ-nele)>0.001):
        if(occ-nele<0):
            i=i+1
        else:
            i=i+1
        en = aen[index[0][i]] 
        occ=(np.multiply((ens<=en).sum(axis=2),kptw)).sum().sum()

    return en


############################################################
def command_line_arg():
    usage = "usage: %prog [options] arg1 arg2"  
    par = OptionParser(usage=usage, version= __version__)

    par.add_option("-d", '--dir', 
            action='append',  dest='dir',
            default=['.'],
            help='location of the PROCAR')

    return  par.parse_args( )

############################################################
def cal_workfunction(dir='.'):
    print(dir)
    Ef = calc_Efermi_at_0K(dir+'/PROCAR') 
    Ef_from_OUTCAR = Efermi_from_OUTCAR(dir+'/OUTCAR')
    vac =vacumm_from_OUTCAR(dir+'/OUTCAR')
    wf=vac-Ef
    print('fermi energy from OUTCAR is {Ef}'.format(Ef=Ef_from_OUTCAR))
    print('fermi energy at 0K is {Ef}'.format(Ef=Ef))
    print('work function is {0} and {1}'.format(wf[0][0],wf[0][1]))
############################################################
if __name__ == '__main__':
    from time import time
    opts, args = command_line_arg()

    #print( opts.dir)
    #print( len(opts.dir))
    if len(opts.dir)==1:
        cal_workfunction('.')
    else:
        for i in opts.dir[1:]:
            if os.path.exists(i):
                cal_workfunction(i)
            else:
                print('Directory {0} doesn\'t exist'.format(i))
