import numpy as np
import pyFAI.geometry
import cryio
import fabio
import os

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

def fluoCorrectionPyfai(poniFile,fluoK=1):
    geo = pyFAI.geometry.Geometry()
    geo.load(poniFile)
    saMap = geo.solidAngleArray()
    return saMap*fluoK

def fluoSub(imageFile,poniFile, fluoK):
    imageArray = cryio.cbfimage.CbfImage(imageFile).array
    fluoArray = fluoCorrectionPyfai(poniFile, fluoK)
    poni = pyFAI.load(poniFile)
    fluoCorr = imageArray - fluoArray
    ext =os.path.splitext(imageFile)[-1]
    direc = os.path.dirname(os.path.realpath(imageFile))
    outfilebase = os.path.basename(imageFile).replace(ext,'fluoSub')
    outfile = f'{direc}/xye/{outfilebase}.xye'
    outfile_2d = f'{direc}/xye/{outfilebase}.edf'
    mask = np.where(imageArray < 0, 1, 0)
    
    poni.integrate1d(data = fluoCorr, filename = outfile,mask = mask,polarization_factor = 0.99,unit = '2th_deg',
                    correctSolidAngle = True, method = 'bbox',npt = 5000, error_model = 'poisson', safe = False)
    clearPyFAI_header(outfile)
    result = poni.integrate2d(data = fluoCorr, filename = outfile_2d,mask = mask,polarization_factor = 0.99,unit = "2th_deg",
                    correctSolidAngle = True, method = 'bbox',npt_rad = 5000, npt_azim = 360, error_model = 'poisson', safe = False)
    bubbleHeader(outfile_2d,*result[:3])

def bubbleHeader(file2d,array2d, tth, eta):
    header = {
    'Bubble_cake' : f'{tth[0]} {tth[-1]} {eta[0]} {eta[-1]}',
    'Bubble_normalized': 1 
    }
    f = fabio.edfimage.EdfImage(data = array2d.transpose(), header = header)
    f.write(file2d)

def clearPyFAI_header(file):
    x,y,e = np.loadtxt(file,unpack = True, comments = '#')
    np.savetxt(file,np.array([x,y,e]).transpose(), '%.6f')