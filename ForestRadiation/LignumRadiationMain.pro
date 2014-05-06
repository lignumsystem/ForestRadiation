######################################################################
# Automatically generated by qmake (2.00a) Thu Aug 17 15:33:51 2006
######################################################################
CONFIG -= app_bundle
CONFIG += qt
QT += xml
TEMPLATE = app
TARGET = lig-radiation 
INCLUDEPATH += . include ../c++adt/include ../stl-lignum/include ../Firmament/include ../stl-voxelspace/include ../LEngine/include ../Pine ../XMLTree  ../Graphics ../Radiation/Shoot
DEPENDPATH += $$INCLUDEPATH
LIBS += -L../c++adt/lib -L../stl-lignum/lib -L../Firmament/lib -L../LEngine/lib -L../stl-voxelspace/lib   -lsky -lL -lvoxel -lLGM  -lcxxadt 
     
macx:LIBS +=  -L../Graphics -lVisual -F/usr/local/Trolltech/Qt-4.1.4/lib -framework GLUT -framework OpenGL
win32:CONFIG += console
#HEADERS += include/CalculateLight.h \
#           include/ScotsPine.h \
#           include/SomeFunctors.h \
#           include/TreeLocations.h \
#           include/StandDescriptor.h \
#	   GrowthLoopRadiation.h \
#           include/BorderForest.h \
#HEADERS +=  ../stl-voxelspace/include/VoxelSpace.h

SOURCES += generate-tree-locations.cc \
           lignum-radiation.cc \
           src/metabolism.cc \
           src/borderforest.cc 

HEADERS += \
    GrowthLoopRadiationI.h \
    GrowthLoopRadiation.h \
    include/CalculateLight.h \
    include/CalculateLightI.h