#This project is for radiation calculations in LIGNUM forest.
#It is placed in the same directory as LignumForest.pro since
#it is only a slight modification of it (growth parts: allocation,
#loop etc are not used in this project).
#The real changes are in the file GrowthLoopRadiationI.h. It
#contains 1) modifications for radiation, and 2) special
#functors EvaluateRadiationForCfTreeSegment_1s, EvaluateRadiationForCfTreeSegment_2
#for radiation.
#Other files are only (at least at the time of creation 13.9.2009) only for a
#separate project:
#lignum-radiation.cc - to include GrowthLoopRadiation.h
#GrowthLoop.h        - to include GrowthLoopRadiationI.h
# + pro files

#target : lig-radiation

TEMPLATE = subdirs
SUBDIRS = LignumForestSubdirs.pro LignumRadiationMain.pro 
