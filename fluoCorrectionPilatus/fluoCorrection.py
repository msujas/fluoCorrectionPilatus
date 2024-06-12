import fluoCorrectionFunctions as fluoCorrectionFunctions



imageFile = r'' #average_gainCorrected file from maskGenerator
poniFile = r'' #poni file used
fluoScale = 7.5e5


fluoCorrectionFunctions.fluoSub(imageFile,poniFile, fluoK=fluoScale)
