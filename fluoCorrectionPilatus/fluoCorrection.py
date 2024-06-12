import fluoCorrectionFunctions as fluoCorrectionFunctions

imageFile = r'C:\Users\kenneth1a\Documents\beamlineData\øystein\pdf\pdf\50_CeBr3_G0p3mm_Bern\pdf\xrd\average/average_gainCorrected.cbf'
poniFile = r'C:\Users\kenneth1a\Documents\beamlineData\øystein\pdf/Si_0_15tilt.poni'
fluoCorrectionFunctions.fluoSub(imageFile,poniFile, fluoK=7.5e5)