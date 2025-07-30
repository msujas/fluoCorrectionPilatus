import numpy as np
from fluoCorrectionPilatus import fluoSub, optimise_fluo

poni = r'C:\Users\kenneth1a\Documents\beamlineData\boJiang_capillaries/Si_0_MD.poni'
file = r'C:\Users\kenneth1a\Documents\beamlineData\boJiang_capillaries\pdf_5\average/5_0001_average_gainCorrected.cbf'
k0 = 1e6
index = 4900
optimise_fluo(file,poni,k0, index)