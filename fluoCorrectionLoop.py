from fluoCorrectionPilatus import fluoSub, optimise_fluo,optimiseFluoBins, optimise_fluoIntegrated
import os
from glob import glob

scale = 5e5
direc = r'C:\Users\kenneth1a\Documents\beamlineData\boJiang_capillaries'
poniFile = r'C:\Users\kenneth1a\Documents\beamlineData\boJiang_capillaries/Si_0_MD.poni' #poni file used
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
            '''
            basefile = os.path.basename(cbf)
            cake = cbf.replace(basefile,f'xye/{basefile}').replace('.cbf','.edf')
            optimise_fluoIntegrated(cake, poniFile, scale, cbf)
            '''
        else:
            fluoSub(cbf, poniFile, scale)
