a script for manual subtraction of fluorescence signal for a Pilatus detector.

First install by doing 'pip install -e .' in the directory.

Use: in the fluoCorrection.py or fluoCorrection_optimise.py script input the average_gainCorrected cbf image (not integrated) and the poni file for the experiment, then choose a scale factor for the fluorescence, this subtracts a background which is the scale factor multiplied by the solid angle of each pixel (calculated with pyFAI) (maximum normalised), then integrates the image in 2D and 1D with pyFAI. Look at the fluoCorr 2d integrated image in Fiewer, and iteratively adjust the scale factor until it looks like it's been properly subtracted: judge by the azimuthal intensity distribution which should be similar to the background sample.

The fluoCorrectionLoop.py script will search for all averaged files in the directory recursively. Can be done with or without the optimiser.

The optimiser works by trying to flatten a bin of the 2d integrated pattern.