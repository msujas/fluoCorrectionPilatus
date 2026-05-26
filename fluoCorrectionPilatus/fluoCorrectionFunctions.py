import numpy as np
import pyFAI.geometry
import cryio
import fabio
import os
from scipy.optimize import least_squares
from functools import partial

def integrate2d(data, mask, ponifile, filename = None, pfactor=  0.85):
    poni = pyFAI.load(ponifile)
    return poni.integrate2d(data= data, mask = mask, filename=filename,polarization_factor = pfactor,unit = "2th_deg",correctSolidAngle = True, 
                            method = 'bbox',npt_rad = 5000, npt_azim = 360, error_model = 'poisson', safe = False) #needs data, mask, filename

def solidAngle(poni1,poni2, d, px, py,psize = 172e-6):
    xpos = px*psize
    ypos = py*psize
    angle1 = np.arctan((np.abs(ypos-poni1)+psize/2)/d) - np.arctan((np.abs(ypos-poni1)-psize/2)/d)
    angle2 = np.arctan((np.abs(xpos-poni2)+psize/2)/d) - np.arctan((np.abs(xpos-poni2)-psize/2)/d)
    return angle1*angle2


def readPoni(poniFile):
    f = open(poniFile)
    string = [line.replace('\n','') for line in f.readlines() if '#' not in line]
    dct={}
    for s in string:
        ssplit = s.split(':')
        dct[ssplit[0]] = ':'.join(ssplit[1:])
        value = dct[ssplit[0]]
        try:
            dct[ssplit[0]] = float(value)
        except ValueError:
            continue
    return dct

def detectorShape(poniFile):
    geo = pyFAI.geometry.Geometry()
    geo.load(poniFile)
    shape = geo.get_shape()
    det = np.empty(shape = (*shape,2))
    for y in range(len(det)):
        for x in range(len(det[0])):
            det[y,x] = [y,x]
    det = det.astype('uint16')
    return det


def solidAngleMap(poniFile):
    det = detectorShape(poniFile)
    poniDct = readPoni(poniFile)
    poni1 = poniDct['Poni1']
    poni2 = poniDct['Poni2']
    d = poniDct['Distance']
    return solidAngle(poni1,poni2,d,det[:,:,1],det[:,:,0])


def fluoCorrection(poniFile, fluoK=1):
    saMap = solidAngleMap(poniFile)
    return fluoK*saMap/np.max(saMap)

def getSAmap(ponifile):
    geo = pyFAI.geometry.Geometry()
    geo.load(ponifile)
    return geo.solidAngleArray()

def getmaps(ponifile, pfactor = 0.85):
    geo = pyFAI.geometry.Geometry()
    geo.load(ponifile)
    return geo.twoThetaArray(), geo.solidAngleArray(), geo.polarization(factor = pfactor)

def fluoCorrectionPyfai(poniFile,fluoK=1):
    return getSAmap(poniFile)*fluoK

def readFile(imageFile):
    ext = os.path.splitext(imageFile)[-1]
    if ext == '.cbf':
        imageArray = cryio.cbfimage.CbfImage(imageFile).array
    elif ext == '.edf' or ext == '.tif':
        imageArray = fabio.open(imageFile).data
    else:
        raise ValueError('image type needs to be .cbf, .edf, or .tif')
    return imageArray

def fluoSub(imageFile,poniFile, fluoK, saveOriginal = False, originalFormat = 'cbf', pfactor = 0.85):
    imageArray = readFile(imageFile)
    fluoArray = fluoCorrectionPyfai(poniFile, fluoK)
    poni = pyFAI.load(poniFile)
    fluoCorr = imageArray - fluoArray
    ext =os.path.splitext(imageFile)[-1]
    direc = os.path.dirname(os.path.realpath(imageFile))
    outfilebase = os.path.basename(imageFile).replace(ext,'fluoSub')
    outfile = f'{direc}/xye/{outfilebase}.xye'
    outfile_2d = f'{direc}/xye/{outfilebase}.edf'
    mask = np.where(imageArray < 0, 1, 0)
    os.makedirs(f'{direc}/xye/', exist_ok = True)
    x,y,e = poni.integrate1d(data = fluoCorr, filename = outfile,mask = mask,polarization_factor = pfactor,unit = '2th_deg',
                    correctSolidAngle = True, method = 'bbox',npt = 5000, error_model = 'poisson', safe = False)
    clearPyFAI_header(outfile)
    result = poni.integrate2d(data = fluoCorr, filename = outfile_2d,mask = mask,polarization_factor = pfactor,unit = "2th_deg",
                    correctSolidAngle = True, method = 'bbox',npt_rad = 5000, npt_azim = 360, error_model = 'poisson', safe = False)
    bubbleHeader(outfile_2d,*result[:3], y,e)
    print(fluoK)
    if saveOriginal:
        match originalFormat:
            case 'cbf':
                #im = fabio.edfimage.EdfImage()
                im = cryio.cbfimage.CbfImage()
                im.array = np.where(fluoCorr<0, -1, fluoCorr)
                im.save(f'{direc}/{outfilebase}.cbf')
            case 'edf':
                im = fabio.edfimage.EdfImage()
                im.data = np.where(fluoCorr<0, -1, fluoCorr)
                im.save(f'{direc}/{outfilebase}.edf')
            case _:
                raise ValueError('originalFormat must be cbf or edf')
    return result

def fluoSub_integrated_base(cakeArray: np.ndarray, polcake:np.ndarray, fluoK:float):
    '''
    The total intensity in the diffraction image can be given by
    It = Isc*P*SA + k*SA
    It - total intensity
    Isc - scattered intensity
    P - beam polarisation effect
    SA - pixel solid angle
    k - fluorescence constant
    So the integrated pattern, without fluorescence correction is given by:
    It/(P*SA) = Isc + k/P
    And so k/P must be subtracted to correct the fluorescence in the already integrated pattern
    '''
    return cakeArray - (fluoK/polcake)

def saveFluosub(fluoSubArray, cakeFile, header):
    im = fabio.edfimage.EdfImage()
    im.data = fluoSubArray
    im.header = header
    im.save(cakeFile.replace('.edf','fluoSub.edf'))

def getMapsintegrated(poniFile, avarrayfile):
    poni = pyFAI.load(poniFile)
    _tthmap, saMap, polmap = getmaps(poniFile)
    avarray = cryio.cbfimage.CbfImage(avarrayfile).array
    mask = np.where(avarray < 0, 1, 0)
    saresult = poni.integrate2d(data = saMap, mask = mask, unit = "2th_deg", method = 'bbox',npt_rad = 5000, npt_azim = 360, 
                              correctSolidAngle=False, error_model = 'poisson', safe = False)
    polresult = poni.integrate2d(data = polmap, mask = mask, unit = "2th_deg", method = 'bbox',npt_rad = 5000, npt_azim = 360, 
                              correctSolidAngle=False, error_model = 'poisson', safe = False)
    return saresult, polresult

def polcorrection(tth, chi, pfactor):
    "0.5 * (1.0 + cos(tth)**2 - factor * cos(2.0 * (chi + axis_offset)) * (1.0 - cos(tth)**2))"
    tthr = tth*np.pi/180
    chir = chi*np.pi/180
    return 0.5*(1.0 + np.cos(tthr)**2 - pfactor * np.cos(2.0 *chir) * (1.0 - np.cos(tthr)**2))

def fluoSub_integrated(cakeFile, fluoK, pfactor):
    cake = fabio.open(cakeFile)
    #cakeArray = cake.data
    header = cake.header 
    tthchi = np.fromstring(header['Bubble_cake'], sep=' ')
    tthlen = int(header['Dim_1'])
    chilen = int(header['Dim_2'])
    tthrange = np.linspace(tthchi[0], tthchi[1], tthlen)
    
    chirange = np.linspace(tthchi[2],tthchi[3],chilen)
    tthmesh, chimesh = np.meshgrid(tthrange,chirange)
    polcake = polcorrection(tthmesh,chimesh,pfactor)
    '''
    saresult, polresult = getMapsintegrated(poniFile,avarrayfile)
    saIntegrated = saresult[0].transpose()
    polIntegrated = polresult[0].transpose()
    polIntegrated = np.where(polIntegrated == 0, 1, polIntegrated)
    '''
    fluosubarray = optimise_fluoIntegrated(cakeFile, polcake, fluoK) #fluoSub_integrated_base(cakeArray, polcake, fluoK)
    fluosub1d = np.nanmean(np.where(fluosubarray<=0, np.nan, fluosubarray),axis=0)
    fluosub1d = np.where(np.isnan(fluosub1d), 0, fluosub1d )
    fluosubarray = np.where(fluosubarray < 0, 0, fluosubarray)
    xyestring = header['Bubble_pattern']
    xyearray = np.fromstring(xyestring, sep=' ')
    lenarray = int(xyearray.shape[0]/3)
    
    x,_y,e = xyearray.reshape((lenarray,3)).transpose()
    
    xyearray = np.array([x,fluosub1d, e]).transpose().flatten()
    xyestring = ' '.join(str(i) for i in xyearray)
    header['Bubble_pattern'] = xyestring
    #saveFluosub(fluosubarray, cakeFile,header)
    bubbleHeader(cakeFile.replace('.edf','fluoSub.edf'), fluosubarray, tthrange, chirange, fluosub1d, e,flip=False)
    fname1d = cakeFile.replace('.edf', 'fluoSubCake.xye')
    np.savetxt(fname1d, np.array([x,fluosub1d, e]).transpose(), fmt = '%.6f')
    return fluosubarray

def optimise_fluoFormula(k0,imagefile, ponifile, index = 4800):
    result = fluoSub(imagefile, ponifile, k0)
    array = result[0]
    arrayline = array[:,index]
    indexes = np.where(arrayline == 0)
    arrayline = np.delete(arrayline,indexes)
    linemean = np.mean(arrayline)
    return (arrayline - linemean)**2

def optimise_fluo(imagefile, ponifile,k0, index = 4800, iters = 20):
    result = least_squares(optimise_fluoFormula,[k0], args = (imagefile, ponifile, index), max_nfev=iters, bounds = (0,np.inf))
    kopt = result['x'][0]
    return fluoSub(imagefile,ponifile,kopt)

def fluointegrated_lsquare(k:float, cake:np.ndarray, polcake:np.ndarray, index = 4800):
    array = fluoSub_integrated_base(cake,polcake, k)
    arrayline = array[:,index]
    arrayline = np.delete(arrayline, np.where(arrayline <= 0))
    linemean = np.mean(arrayline)
    return (arrayline - linemean)**2

def optimise_fluoIntegrated(cakefile, polcake, k0, index = 4800, iters = 1000):
    cakedata = fabio.open(cakefile)
    cakearray = cakedata.data
    header = cakedata.header
    '''
    saresult, polresult = getMapsintegrated(ponifile,avarrayfile)
    saintegrated = saresult[0].transpose()
    polIntegrated = polresult[0].transpose()
    polIntegrated = np.where(polIntegrated == 0, 1, polIntegrated)
    '''
    result = least_squares(fluointegrated_lsquare, [k0],args = (cakearray, polcake, index), max_nfev=iters)
    kopt = result['x'][0]
    print(kopt)
    fluosub = fluoSub_integrated_base(cakearray, polcake, kopt)
    #saveFluosub(fluosub,cakefile,header)
    return fluosub

def rebin(array, nbins):
    binsize = int(len(array)/nbins)
    return np.array([np.mean(array[i*binsize:(i+1)*binsize]) for i in range(nbins)])
    
def fluobinPrep(avfile, ponifile):
    array = readFile(avfile)
    tthmap, saMap, polmap = getmaps(ponifile)
    return array, tthmap, saMap, polmap

def fluoSubBins(fluoK, array, tthmap, saMap, polmap, nbins, index):
    if index > nbins:
        raise ValueError('index must be less than the number of bins')
    fluosubarray = array - (saMap*fluoK)
    fluosubarray = fluosubarray/(saMap*polmap)
    binsize = np.max(tthmap)/nbins
    binarray = ((tthmap+binsize/2)*(nbins-1)//np.max(tthmap)).astype(int)
    arrayline = fluosubarray[np.where((binarray == index) & (array >= 0))]
    #arrayline = rebin(arrayline, 200)
    linemean = np.mean(arrayline)
    modifier = np.where(arrayline >= 0, 0, arrayline**2) #np.bitwise_and(arrayline.astype(np.int32),-2**32)*100/2**32 #punish negative values
    return (arrayline - linemean)**2 + modifier

def optimiseFluoBins(avfile, ponifile,k0, nbins, index, saveOriginal=False):
    #This is the default version
    if index > nbins:
        raise ValueError('index must be less than the number of bins')
    array, tthmap, saMap, polmap = fluobinPrep(avfile,ponifile)
    result = least_squares(fluoSubBins, k0, args = (array, tthmap, saMap, polmap,nbins, index), bounds = (0,np.inf))
    kopt = result['x'][0]
    print(kopt)
    return fluoSub(avfile, ponifile, kopt, saveOriginal=saveOriginal)

def bubbleHeader(file2d,array2d, tth, eta, y, e, flip=True):
    xye = np.array([tth,y,e]).transpose().flatten()
    xyestring = ' '.join([str(i) for i in xye])
    header = {
    'Bubble_cake_version' : 3,
    'Bubble_cake' : f'{tth[0]} {tth[-1]} {eta[0]} {eta[-1]}',
    'Bubble_normalized': 1 ,
    'Bubble_pattern': xyestring
    }
    spacing = -1 if flip else 1
    f = fabio.edfimage.EdfImage(data = array2d[::spacing,:], header = header)
    f.write(file2d)

def clearPyFAI_header(file):
    x,y,e = np.loadtxt(file,unpack = True, comments = '#')
    np.savetxt(file,np.array([x,y,e]).transpose(), '%.6f')

def getChiMap(ponifile):
    geo = pyFAI.geometry.Geometry()
    geo.load(ponifile)
    return geo.chiArray()

import matplotlib.pyplot as plt
def plotBin(avfile, ponifile, nbins, index, fluoK = 0):
    array = readFile(avfile)
    tthmap, samap, polmap = getmaps(ponifile)
    array = array - (samap*fluoK)
    array = array/(samap*polmap)
    chiarray = getChiMap(ponifile)
    binsize = np.max(tthmap)/nbins
    binarray = ((tthmap+binsize/2)*(nbins-1)//np.max(tthmap)).astype(int)
    where = np.where((binarray == index ) & (array >= 0))
    arrayline = array[where]
    chiline = chiarray[where]
    plt.plot(chiline, arrayline, 'o', markersize = 1)
    plt.xlabel('chi (rad)')
    plt.ylabel('intensity')
    plt.show()