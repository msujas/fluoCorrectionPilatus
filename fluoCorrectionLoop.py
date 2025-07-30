from fluoCorrectionPilatus import fluoSub, optimise_fluo
import os
from glob import glob

scale = 1e6
direc = r'C:\Users\kenneth1a\Documents\beamlineData\boJiang_capillaries'
poniFile = r'C:\Users\kenneth1a\Documents\beamlineData\boJiang_capillaries/Si_0_MD.poni' #poni file used
optimise = True
index = 4800
iters = 10
for root, dirs,files in os.walk(direc):
    if not 'average' in root or 'xye' in root:
        continue
    cbfs = glob(f'{root}/*.cbf')
    for cbf in cbfs:
        if optimise:
            optimise_fluo(cbf, poniFile,scale,index,iters)
        else:
            fluoSub(cbf, poniFile, scale)
