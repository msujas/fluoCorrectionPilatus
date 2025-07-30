from fluoCorrectionPilatus import fluoSub



imageFile = r'' #average_gainCorrected file from maskGenerator
poniFile = r'' #poni file used
fluoScale = 2e6


fluoSub(imageFile,poniFile, fluoK=fluoScale)
