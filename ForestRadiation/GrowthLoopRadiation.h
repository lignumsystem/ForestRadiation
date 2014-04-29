#ifndef GROWTH_LOOPRADIATION_H
#define GROWTH_LOOPRADIATION_H
#include <cmath>
#include <cstdio>
#include <fstream>
#include <utility>
#include <mathsym.h>
#include <Point.h>
#include <PositionVector.h>
#include <Bisection.h>
#include <VoxelSpace.h>
#include <TreeLocations.h>
#include <CalculateLight.h>
#include <VoxelSpace.h>
#include <XMLTree.h>
#include <StandDescriptor.h>
#include <VoxelSpace.h>
#include <BorderForest.h>
#include <star_mean.h>

template <class TREE, class TS, class BUD>//,class LSYSTEM>
  class GrowthLoop{
 public:
  GrowthLoop()
    :vs(NULL),verbose(false),iterations(0),
    num_parts(1.0), tree_distance(0.0),
    writevoxels(false), generate_locations(false),
    no_trees(0), wood_voxel(true), gap_radius(0.0), middle_stand(0.0,0.0),
    target_tree_rad(0.5) {}
    ~GrowthLoop();
    void usage()const;
    void checkCommandLine(int argc, char** argv)const;
    void parseCommandLine(int argc, char** argv);
    void setVerbose(bool wordy){verbose=wordy;}
    bool isVerbose()const{return verbose;}
    void initRan3(int init)const{
      int neg_init=-abs(init);
      ran3(&neg_init);
    }
    void createTrees();
    void initializeTrees();
    void initializeVoxelSpace();
    void initializeFunctions();
    void initializeGrowthLoop();
    void growthLoop();
    //    void increaseXi(TREE& t);
   void setTreeLocations();
    void photosynthesis(TREE& t);
    void respiration(TREE& t);
    void writeProductionBalance(TREE& t,const string& file)const;
    void setVoxelSpaceAndBorderForest();
    void calculateRadiation();
    void calculateRadiationOnlySelf();
    void getTreesAndPositions();
    StandDescriptor<TREE>& getStand() {return stand;}
    bool getOnlySelf() {return only_self;}
    bool getWriteOnlyFile() {return write_only_file;}
    bool getManyTrees() {return many_trees;}
    bool getTreesFromFile() {return trees_from_file;}
 private:
    vector<TREE*> vtree;//vector of trees
    vector<pair<double,double> > locations;
    vector<ofstream*> vdatafile;
    VoxelSpace *vs;//The voxel space
    StandDescriptor<TREE> stand;
    BorderForest border_forest;
    bool verbose;
    int iterations;
    int num_parts;
    string location_file;
    double tree_distance;//Minimum distance between trees in meters

    double treeAf; //Foliage area of the tree
    DCLData dcl;//Diameter and heigth at the crown base.
    string metafile;
    string voxelfile;
    string xmlfile; //XML file where the tree can be saved and restored from
    bool writevoxels;
    ParametricCurve K;//The 'K' function 
    bool generate_locations;
    int no_trees;
    bool noWoodVoxel;
    bool wood_voxel;
    string phprodfile;
    LGMdouble voxboxside;
    bool dump_self;
    string input_tree_file;
    LGMdouble k_border_conifer;
    LGMdouble gap_radius;
    pair<double,double> middle_stand;
    LGMdouble target_tree_rad;
    bool evaluate_LAI;
    int rad_method;
    bool calculate_star;
    int num_star_calculations;
    bool voxel_tree;
    bool mean_inclination;
    LGMdouble rel_tt_gap;
    bool print_box_cf_data;
    string box_cf_data_file;
    bool box_dir_effect;
    bool tree_info;
    LGMdouble constant_star; //If this > 0, k = constant_star else k = star_mean
    bool correct_star;       //If star_eq -> star correction is done
    bool one_time;        //If production balance is written to -phprodfile in trunc-mode or app-mode
    bool only_self;       //If radiation is analyzed only for -inputTree
    bool many_trees;      //If many shading trees
    string many_tree_file;
    vector<string> tree_files;
    bool write_only_file;
    bool trees_from_file;
    string trees_pos_file;
    bool radius_only;
    LGMdouble max_radius;
};


class Compartments{
 public:
 Compartments() : buds(0), segments(0), bps(0), axes(0) {}
 Compartments(const Compartments& copy_from) : buds(copy_from.buds),
    segments(copy_from.segments), bps(copy_from.bps), axes(copy_from.axes) {}
 
    Compartments& operator += (const Compartments& add_this) {
      buds += add_this.buds;
      segments += add_this.segments;
      bps += add_this.bps;
      return *this;
   }

  int buds;
  int segments;
  int bps;
  int axes;
};

template <class TS,class BUD>
  class AddUpCompartments {
 public:
  Compartments&
    operator()(Compartments& n,TreeCompartment<TS,BUD>* tc)const
    {
      if (dynamic_cast<Axis<TS,BUD>*>(tc)){
	n.axes++;
      }

      else if (dynamic_cast<BranchingPoint<TS,BUD>*>(tc)){
	n.bps++;
      }

      else if (dynamic_cast<TS*>(tc))
	n.segments++;

      else if (dynamic_cast<Bud<TS,BUD>*>(tc))
	n.buds++;

      else
	;
      return n;
    }
};


#endif

#include <GrowthLoopRadiationI.h>
