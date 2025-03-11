import fluoCorrectionFunctions



imageFile = r'' #average_gainCorrected file from maskGenerator
poniFile = r'' #poni file used
fluoScale = 5e2


fluoCorrectionFunctions.fluoSub(imageFile,poniFile, fluoK=fluoScale)
