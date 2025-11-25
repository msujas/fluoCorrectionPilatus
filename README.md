a script for manual subtraction of fluorescence signal for a Pilatus detector.

First install by doing 'pip install -e .' in the directory.

Use: in the fluoCorrection.py or fluoCorrection_optimise.py script input the average_gainCorrected cbf image (not integrated) and the poni file for the experiment, then choose a scale factor for the fluorescence, this subtracts a background which is the scale factor multiplied by the solid angle of each pixel (calculated with pyFAI) (maximum normalised), then integrates the image in 2D and 1D with pyFAI. Look at the fluoCorr 2d integrated image in Fiewer, and iteratively adjust the scale factor until it looks like it's been properly subtracted: judge by the azimuthal intensity distribution which should be similar to the background sample.

There is also a command line program: fluoCorrection.
```
usage: fluoCorrection [-h] [-p PONI] [-k K0] [-i INDEX] [-r] [-so] [file]

positional arguments:
  file                 filename or directory to run in (default current directory)

options:
  -h, --help           show this help message and exit
  -p, --poni PONI      poni file for measurement
  -k, --k0 K0          starting scaling factor for fluorescence to use
  -i, --index INDEX    xrd bin to try to flatten in optimisation (default 4800)
  -r, --recurse        run recursively (only averaged files)
  -so, --saveOriginal  save the non-integrated fluo subtracted image
```
The fluoCorrectionLoop.py script will search for all averaged files in the directory recursively. Can be done with or without the optimiser.

The optimiser works by trying to flatten a bin of the 2d integrated pattern.

Non corrected:
<img width="1473" height="585" alt="image" src="https://github.com/user-attachments/assets/3420fbf7-f972-4ff9-b6f5-d4258f16aa8c" />

Corrected:
<img width="1474" height="584" alt="image" src="https://github.com/user-attachments/assets/a81253c3-af7d-47ae-ab54-a94c46430311" />
