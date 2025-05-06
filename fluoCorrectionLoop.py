from fluoCorrectionPilatus import fluoSub
import os
from glob import glob

scale = 2e6
direc = r'F:\NOV-2024'
poniFile = r'F:\NOV-2024/Si000_15tilt-IZ.poni' #poni file used

for root, dirs,files in os.walk(direc):
    if not 'average' in root or 'xye' in root:
        continue
    cbfs = glob(f'{root}/*.cbf')
    for cbf in cbfs:
        fluoSub(cbf, poniFile, scale)
