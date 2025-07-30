from fluoCorrectionPilatus import fluoSub



imageFile = r'F:\NOV-2024\pdf_CSnO_DME_1C_0p5V\average/CSnO_DME_1C_0p5V_001_average_gainCorrected.cbf' #average_gainCorrected file from maskGenerator
poniFile = r'F:\NOV-2024/Si000_15tilt-IZ.poni' #poni file used
fluoScale = 2e6


fluoSub(imageFile,poniFile, fluoK=fluoScale)
