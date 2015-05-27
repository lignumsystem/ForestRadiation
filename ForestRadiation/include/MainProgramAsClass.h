#ifndef MAIN_PROGRAM_AS_CLASS_H
#define MAIN_PROGRAM_AS_CLASS_H
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
#include <SomeFunctors.h>
#include <CalculateLight.h>
#include <VoxelSpace.h>
#include <XMLTree.h>
#include <StandDescriptor.h>
#include <BorderForest.h>
#include <star_mean.h>

template <class TREE, class TS, class BUD>
  class MainProgramAsClass{
 public:
  MainProgramAsClass()
    :vs(NULL),verbose(false),
    num_parts(1.0), tree_distance(0.0),
    writevoxels(false), generate_locations(false),
    no_trees(0), wood_voxel(true), gap_radius(0.0), middle_stand(0.0,0.0),
    target_tree_rad(0.5) {}
    ~MainProgramAsClass();
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
    void initializeMainProgramAsClass();
    void growthLoop();
    //    void increaseXi(TREE& t);
   void setTreeLocations();
    void photosynthesis(TREE& t);
    void respiration(TREE& t);
    void writeProductionBalance(TREE& t,const string& file)const;
    void setVoxelSpaceAndBorderForest();
    void calculateRadiation();
    void calculateRadiationOnlySelf();
    void calculateRadiationToPoint();
    void getTreesAndPositions();
    void getTreesAndPositionsPeriodicalBoundary();
    void createTargetTree();
    StandDescriptor<TREE>& getStand() {return stand;}
    bool getOnlySelf() {return only_self;}
    bool getWriteOnlyFile() {return write_only_file;}
    bool getManyTrees() {return many_trees;}
    bool getTreesFromFile() {return trees_from_file;}
    bool getTreesFromFilePeriodical() {return trees_from_file_periodical;}
    bool getOnlyPositions() {return only_positions;}
    void getThis(int& i) {
      i++;
    }
 private:
    vector<TREE*> vtree;//vector of trees
    vector<pair<double,double> > locations;
    vector<ofstream*> vdatafile;
    VoxelSpace *vs;//The voxel space
    StandDescriptor<TREE> stand;
    BorderForest border_forest;
    bool verbose;
    int num_parts;
    string location_file;
    double tree_distance;//Minimum distance between trees in meters

    double treeAf; //Foliage area of the tree
    DCLData dcl;//Diameter and heigth at the crown base.
    string xmlfile; //XML file where the tree can be saved and restored from
    bool writevoxels;
    ParametricCurve K;//The 'K' function 
    bool generate_locations;
    int no_trees;
    bool wood_voxel;
    bool dir_star; // variable to determine whether to activate the calculations of the directional STAR calculations or not.
    string resultfile;
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
    bool one_time;        //If production balance is written to -resultfile in trunc-mode or app-mode
    bool only_self;       //If radiation is analyzed only for -inputTree
    bool many_trees;      //If many shading trees
    string many_tree_file;
    vector<string> tree_files;
    bool write_only_file;
    bool trees_from_file;
    string trees_pos_file;
    bool radius_only;
    LGMdouble max_radius;
    bool calculateDirectionalStar; //To calculate directional STAR values
    LGMdouble vs_x, vs_y, vs_z;    //Dimensions of the voxelspace
    bool zero_woody_radius;        //If woody radius of all shading trees is about 0
    string trees_pos_periodical_file;
    bool trees_from_file_periodical;
    bool only_positions;
    string only_positions_file;


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


vector<pair<LGMdouble,LGMdouble> > translateCoordinates(const LGMdouble& x, const LGMdouble& y,
					       const LGMdouble& tr_x, const LGMdouble& tr_y)
{
  vector<pair<LGMdouble,LGMdouble> > points(8);
  points[0] = pair<LGMdouble,LGMdouble>(x, y + tr_y);
  points[1] = pair<LGMdouble,LGMdouble>(x, y - tr_y);
  points[2] = pair<LGMdouble,LGMdouble>(x + tr_x, y - tr_y);
  points[3] = pair<LGMdouble,LGMdouble>(x + tr_x, y);
  points[4] = pair<LGMdouble,LGMdouble>(x + tr_x, y + tr_y);
  points[5] = pair<LGMdouble,LGMdouble>(x - tr_x, y - tr_y);
  points[6] = pair<LGMdouble,LGMdouble>(x - tr_x, y);
  points[7] = pair<LGMdouble,LGMdouble>(x - tr_x, y + tr_y);

  return points;
}


#endif

#include <MainProgramAsClassI.h>
