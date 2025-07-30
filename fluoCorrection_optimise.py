import numpy as np
from fluoCorrectionPilatus import fluoSub, optimise_fluo

file = r''
poni = r''
k0 = 5e5
index = 4900
iters = 10
optimise_fluo(file,poni,k0, index, iters)