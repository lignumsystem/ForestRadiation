//Include Lignum implementation
#include <Lignum.h>
#include <MainProgramAsClass.h> 
//Include the implementation of the tree segment and bud
#include <ScotsPine.h>

#if defined (__APPLE__) || defined(__MACOSX__)
#include <VisualFunctor.h>
//Impelements VisualizeLGMTree
#include <GLSettings.h>
#include <OpenGLUnix.h>
#include <LGMVisualization.h>
#endif
//Includes all kinds of stuff, turtle graphics etc.
#include <lengine.h>

//and for pine, see also pine9bp.L in lsys.
namespace Pine{
#include <LSystem.h>

}

bool no_compartments;
string comp_tree;

int ran3_seed = -9648383; //-23797843;

int main(int argc, char** argv)
{
  MainProgramAsClass<ScotsPineTree,ScotsPineSegment,ScotsPineBud > main_program;

  ran3(&ran3_seed);

  main_program.setVerbose(true);
  //Check and parse command line, the  command line includes switches to control
  //the simulation

  main_program.checkCommandLine(argc, argv);

  main_program.parseCommandLine(argc,argv);
  if(main_program.getOnlySelf()) {
    main_program.calculateRadiationOnlySelf();
  }
  else {
    if(main_program.getTreesFromFile()) {
      main_program.getTreesAndPositions();
    }
    else if(main_program.getTreesFromFilePeriodical()) {
      main_program.getTreesAndPositionsPeriodicalBoundary();
    }
    else {
      main_program.setTreeLocations();
      main_program.createTrees();  //to locations set above
   
      //Write only the file about positions & trees in this case
      if(main_program.getWriteOnlyFile() && main_program.getManyTrees())
	exit(0);
    }
  }
  if(!main_program.getOnlyPositions()) {
    main_program.initializeVoxelSpace();

    main_program.setVoxelSpaceAndBorderForest();

    // main_program.createTargetTree();

    main_program.calculateRadiation();
  }
  else {
    main_program.calculateRadiationToPoint(argc,argv);
  }
  exit(0);
}
