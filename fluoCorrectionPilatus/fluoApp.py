from fluoCorrectionPilatus import *
import argparse
from glob import glob

def parseArgs():
    filedct = {'filename':'image file to perform fluorescence correction on',
               'directory': 'directory to perform fluorescence correction loop in'}
    parser = argparse.ArgumentParser()
    parser.add_argument('file', type=str, help = 'filename or directory to run in (default current directory)', nargs='?', default='.')
    parser.add_argument('-p','--poni', type = str, help = 'poni file for measurement')
    parser.add_argument('-k','--k0', default=5e5, type = int, help= 'starting scaling factor for fluorescence to use')
    parser.add_argument('-i','--index',default=4800,type = int, help = 'xrd bin to try to flatten in optimisation (default 4800)')
    parser.add_argument('-r','--recurse',action='store_true',help='run recursively (only averaged files) (flag argument)')
    parser.add_argument('-so','--saveOriginal',action='store_true',help='save the non-integrated fluo subtracted image (flag argument)')
    parser.add_argument('-no', '--nooptimisation',action='store_true', help= 'turn off optimisation (just take starting k value)')
    args = parser.parse_args()
    filename = args.file
    poni = args.poni
    k0 = args.k0
    index = args.index
    recurse = args.recurse
    so = args.saveOriginal
    no = args.nooptimisation
    return filename, poni, k0, index, recurse, so, no


def runOptimise(file,poni, k0,index, saveOriginal=False, nooptimise = False):
    #file,poni, k0,index = parseArgs() #average image file
    print(file)
    if nooptimise:
        fluoSub(file, poni, k0, saveOriginal=saveOriginal)
        return
    optimiseFluoBins(file, poni, k0, 5000, index, saveOriginal=saveOriginal)

def runOptimiseDir(direc, poniFile, k0, index, saveOriginal=False, nooptimise = False):
    files = glob(f'{direc}/*.cbf')
    for file in files:
        runOptimise(file, poniFile,k0, index, saveOriginal=saveOriginal, nooptimise=nooptimise)

def runOptimiseRecurse(direc, poniFile, k0, index, saveOriginal=False):
    #direc, poniFile, k0, index = parseArgs()
    for root, dirs,files in os.walk(direc):
        if not 'average' in root or 'xye' in root:
            continue
        cbfs = glob(f'{root}/*.cbf')
        for cbf in cbfs:
            print(cbf)
            optimiseFluoBins(cbf, poniFile, k0, 5000, index, saveOriginal=saveOriginal)

def run():
    file, poniFile, k0, index, recurse, so, no = parseArgs()
    if os.path.isdir(file):
        if recurse:
            runOptimiseRecurse(file, poniFile, k0, index,so)
        else:
            runOptimiseDir(file,poniFile, k0, index,so, no)
    else:
        runOptimise(file, poniFile, k0, index,so, no)

def cakesubargpars():
    parser = argparse.ArgumentParser()
    parser.add_argument('cakefile',help='cake file to optimise fluo subtraction')
    parser.add_argument('-k','--k0', help='starting fluorescence subtraction constant',default=10**5,type=float)
    parser.add_argument('-p','--polarization_factor',help='polarization factor',default=0.85, type=float)
    a= parser.parse_args()
    k0 = a.k0
    cakefile = a.cakefile
    pfactor = a.polarization_factor
    return cakefile, k0, pfactor


def cakesub():
    cakefile, k0, pfactor = cakesubargpars()
    return fluoSub_integrated(cakefile, k0,pfactor)