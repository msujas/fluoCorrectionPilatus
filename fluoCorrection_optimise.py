import numpy as np
from fluoCorrectionPilatus import fluoSub, optimise_fluo

poni = r''
file = r''
k0 = 1e6
index = 4900
optimise_fluo(file,poni,k0, index)