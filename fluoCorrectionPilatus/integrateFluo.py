import maskGeneratorBM31 as mg
from .fluoCorrectionFunctions import optimiseFluoBins
import argparse
import os

def integrateFluo(direc, poni, mask, gainfile, scale, k0=10**5):
    '''
    Docstring for integrateFluo
    
    :param direc: directory to run in
    use the mask generator to integrate a stack of images in a directory, then apply a fluorescence correction to the integrated image
    '''
    mg.runAll(direc,poni,mask,gainfile,scale)
    avfiles = mg.getAvFiles(direc)
    optimiseFluoBins(avfiles[1], poni, k0, 5000, 4800)

def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('directory',type = str,default=os.path.curdir, help = 'directory to run in')
    parser.add_argument('-p', '--poni', type = str, help ='poni file')
    parser.add_argument('-m','--mask', type = str, help = 'mask file')
    parser.add_argument('-g','--gainfile', type = str, help = 'gain file to use')
    parser.add_argument('-s','--integrationScale',type = int, default=10**9, help = 'scale to apply to integrated data, default 10**9')
    parser.add_argument('-k','--k0', type = int, default=10**5, help = 'starting fluorescence scaling factor, default 10**5')
    args = parser.parse_args()
    direc = args.directory
    poni = args.poni
    mask = args.mask
    gainfile= args.gainfile
    scale = args.integrationScale
    k0 = args.k0
    return direc, poni, mask,gainfile,scale,k0


def main():
    direc, poni,mask,gainfile, scale, k0 = parseArgs()
    integrateFluo(direc,poni,mask,gainfile,scale, k0)



