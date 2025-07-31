import numpy as np
from fluoCorrectionPilatus import optimise_fluo, optimiseFluoBins

file = r''
poni = r''
k0 = 5e5
index = 4800
iters = 10
#optimise_fluo(file,poni,k0, index, iters)
optimiseFluoBins(file, poni, k0, 5000, 4800)
