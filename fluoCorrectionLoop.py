from fluoCorrectionPilatus import fluoSub, optimise_fluo,optimiseFluoBins
import os
from glob import glob

scale = 5e5
direc = r''
poniFile = r'' #poni file used
optimise = True #takes longer, but maybe useful if different samples have different amounts of fluorescence
index = 4800
iters = 10
for root, dirs,files in os.walk(direc):
    if not 'average' in root or 'xye' in root:
        continue
    cbfs = glob(f'{root}/*.cbf')
    for cbf in cbfs:
        if optimise:
            #optimise_fluo(cbf, poniFile,scale,index,iters)
            optimiseFluoBins(cbf, poniFile, scale, 5000, index)
        else:
            fluoSub(cbf, poniFile, scale)
