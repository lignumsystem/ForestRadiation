#ifndef MAIN_PROGRAM_AS_CLASSI_H
#define MAIN_PROGRAM_AS_CLASSI_H

class SetSf;  //This is defined at the end of this file (declared here before use)
class SetSf2;
class CalculateSTAR;

extern int ran3_seed;
extern bool no_compartments;
extern string comp_tree;

//This is MainProgramAsClassI.h for radiation calculations. See LignumRadiation.pro for information
// about LignumRadiation project.

template<class TREE, class TS, class BUD>
MainProgramAsClass<TREE,TS,BUD>::~MainProgramAsClass()
{
  for (unsigned int i = 0; i < vtree.size(); i++){
    delete vtree[i];
  }

  for (unsigned int i = 0; i < vdatafile.size(); i++){
    vdatafile[i]->close();
    delete vdatafile[i];
  }
  
}

template<class TREE, class TS, class BUD>
  void MainProgramAsClass<TREE,TS,BUD>::usage()const
{
  cout << "Usage:  ./lig-radiation [-numParts <parts>]  [-treeDist <dist>]" <<endl;
  cout << "[-xml <filename>] [-writeVoxels] [-noWoodVoxel]" <<endl;
  cout << "[-treeFile <filename>] [-generateLocations  <num>] [-treeLocations <file>]" << endl;
  cout << "[-resultfile <file>] [-Voxboxside <value>]" << endl;
  cout << "[-dumpSelf] [-inputTree <filename>] [-kBorderConifer <value>] [-GapRadius <value>]" << endl;
  cout << "[-targetTreeRad <value>] [-evaluateLAI] [-radMethod <num>] [-calculateSTAR <num>] [-calculateDirectionalStar] " << endl;
  cout << "[-voxelTree] [-boxDirEffect] [-treeInfo] [-segmentInfo <file>]" << endl;
  cout << "[-correctSTAR] [-constantSTAR <value>] [-appendMode] [-self] [-manyTrees <file>]" << endl;
  cout << "[-writeOnlyFile] [-getTreesPos <file>] [-radiusOnly <m>]" << endl;
  cout << "[-X <value>] [-Y <value>] [-Z <value>] [-evaluateLAI] [-zeroWoodyRadius]" << endl;
  cout << "-generateLocations <num>  In this case <num> trees will be generated to random locations. If this" << endl;
  cout << "          is not on, tree locations will be read from file Treelocations.txt. This file can be changed" << endl;
  cout << "          by -treeLocations <file>. If location file is not found program stops." << endl;
  cout << "-noWoodVoxel            Woody parts are NOT dumped to voxels (default is are dumped)" << endl;
  cout << "-calculateDirectionalStar If directional star needs to be calculated then use true (default = false)  "<<endl;
  cout << "-treeDist <dist>          Minimum distance between two trees (default = 0), works only with -generateLocations." << endl;  
  cout << "-numParts <parts>         Segments can be dumped to voxels by parts (i.e. they may belong to different voxels," << endl;
  cout << "-targetTree <num>         Any one of the trees can be identified as target tree (default = 0)" << endl;
  cout << "-resultfile <file>        File to store radiation calculations"  << endl;
  cout << "-Voxboxside <value>     Side length of voxel box - default is 0.2 m"  << endl;
  cout << "-dumpSelf              If the subject tree is dumped to vox-space (default no)" << endl;
  cout << "-inputTree <filename>  Input tree for radiation calculations" << endl;
  cout << "-kBorderConifer <value>       Extinction coefficient for conifer foliage in border forest (default = 0.14)" << endl;
  cout << "-targetTreeRad <value>        Distance of furthest segments from stem of the target tree." << endl;
  cout << "-evaluateLAI           If vertical LAI distn of generated forest is evaluated; if it is, nothing else is done" << endl;
  cout << "-radMethod <num>       The version of radiation calculations to be used." << endl;
  cout << "-calculateSTAR <num>   STAR values are clculated for each shoot in the first tree in in the tree vecor and"
          "                       program stops. <num> specifies no. runs in Monte Carlo evaluation of STAR for each"
          "                       shoot. Results are written in file STAR.dat. Note that if <num> is large,"
          "                       the program runs a long time."   << endl; 
  cout << "-voxelTree             In this case a tree is made which has one segment in each voxel."
          "                       NOTE Treelocations file or -generate MUST be such that there is only one tree."
          "                       Side length of voxelbox _MUST_ be > 0.15 (length of segments in voxels), otherwise program"
          "                       will crash."
       << endl;
  cout << "-PrintBoxCfData <file> Writes voxelboxcontents to file with VoxelSpace.PrintBoxCfData(). Then stops"
       << endl;
  cout << "-boxDirEffect          If effect of mean direction of segments in box considered (default = no)"
       << endl;
  cout << "-treeInfo              Writes H, Dbh, Dbase, Hcrown_base Dcrown_base Wf Af NAD of tree on console and stops." << endl;
  cout << "-correctSTAR           If the descrepancy with STATR from eq. and STAR (ca. STAR = -0.041 + 1.056*STAR_eq) is corr."
       << endl;
  cout << "-constantSTAR <value>  STAR has constant value <value> (may be corrected by -correctSTAR)."  << endl;
  cout << "-appendMode            If results are written to -resultfile in trunc-mode or app-mode." << endl;
  cout << "                       trunc mode is default, in app-mode information about voxels is written also." << endl;
  cout << "-self                  Calculates the radiation conditions only for -inputTree" << endl;
  cout << "-manyTrees <file>      Many shading trees, they are given in <file>. In this case -self has no effect" << endl; 
  cout << "-writeOnlyFile         Writes only positions & trees to runfile.dat and exits, requires -manyTrees" << endl;
  cout << "-getTreesPos <file>    Reads trees (xml files) and their positions from file" << endl;
  cout << "-radiusOnly <r>        Only trees at max distance r m are used in the calculation." << endl;
  cout << "-X, -Y, -Z <value>     X, y, z sidelengths of VoxelSpace (defaults are 10 m, 10 m, 3 m)" << endl;
  cout << "-zeroWoodyRadius       The woody radii in all shading trees set to ca. zero: SetValue(*ts,LGAR, 0.0001);" << endl; 
  cout  << endl;
}


template<class TREE, class TS, class BUD>
void MainProgramAsClass<TREE,TS,BUD>::checkCommandLine(int argc, char** argv)const
{
  //At least three  mandatory arguments required 
  if (argc < 1){
    usage();
    exit(0);
  }
  else if (verbose){
    cout << "Command line O.K." <<endl;
  } 
}

template<class TREE, class TS, class BUD>
  void MainProgramAsClass<TREE,TS,BUD>::parseCommandLine(int argc, char** argv)
{
  if (verbose){
    cout << "parseCommandLine begin" <<endl;
  }

   checkCommandLine(argc,argv);

  string clarg;

  //Read/set parameters of the voxelspace and the are trees are in
  voxboxside = 0.2;     //default value = 0.2 m
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-Voxboxside", clarg))
    voxboxside = atof(clarg.c_str());

  clarg.clear();
  //Side lengths of voxelspace, used sometimes (e.g. setTreeLocations())
  //to define the ares the trees are in
  //Default sixe of the voxelspace is 10 x 10 x 3 m3
  vs_x = 10.0;     
  if (ParseCommandLine(argc,argv,"-X", clarg))
    vs_x = atof(clarg.c_str());
 
  clarg.clear();
  vs_y = 10.0;     
  if (ParseCommandLine(argc,argv,"-Y", clarg))
    vs_y = atof(clarg.c_str());

  clarg.clear();
  vs_z = 3.0;     
  if (ParseCommandLine(argc,argv,"-Z", clarg))
    vs_z = atof(clarg.c_str());

  middle_stand.first = vs_x/2.0;
  middle_stand.second = vs_y/2.0;
  
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-treeDist", clarg))
    tree_distance = atof(clarg.c_str());
  
  //Number of segment parts used to assess Qabs
  num_parts = 1;
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-numParts", clarg))
    num_parts = atoi(clarg.c_str());
  
  wood_voxel = true;
  if (CheckCommandLine(argc,argv,"-noWoodVoxel"))
    wood_voxel = false;

  //Initialize ran3
  //ran3_seed is a global variable
  int s_ini;
  if (ParseCommandLine(argc,argv,"-seed", clarg)){
    if (clarg.length() > 0){
      s_ini = atoi(clarg.c_str());
      ran3_seed = -abs(s_ini);
    }
  }

  generate_locations = false;
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-generateLocations", clarg)) {
    generate_locations = true;
    no_trees = atoi(clarg.c_str());
  }
  else {
    clarg.clear();
    if (ParseCommandLine(argc,argv,"-treeLocations", clarg)) 
      location_file = clarg;
    else
      location_file = "Treelocations.txt";
  }

  tree_distance = 0.0;
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-treeDist", clarg))
    tree_distance = atof(clarg.c_str());

  resultfile = "result.dat";
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-resultfile", clarg))
    resultfile = clarg;

  dump_self = false;
  if (CheckCommandLine(argc,argv,"-dumpSelf"))
    dump_self = true;

  clarg.clear();
  if (ParseCommandLine(argc,argv,"-inputTree", clarg)) 
    input_tree_file = clarg;
  else {
    //    cout << "No input tree! " << endl;
    //    exit(-1);
    input_tree_file = "";
  }

  k_border_conifer = 0.14;
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-kBorderConifer", clarg))
    k_border_conifer = atof(clarg.c_str());

  gap_radius = 0.0;
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-GapRadius", clarg))
    gap_radius = atof(clarg.c_str());

  target_tree_rad = 0.5;
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-targetTreeRad", clarg))
    target_tree_rad = atof(clarg.c_str());
  
  evaluate_LAI = false; 
  if (CheckCommandLine(argc,argv,"-evaluateLAI")) {
    evaluate_LAI = true;
  }
  
  rad_method = 3;
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-radMethod", clarg))
    rad_method = atoi(clarg.c_str());

  clarg.clear();
  calculate_star = false;
  num_star_calculations = 1;
  if (ParseCommandLine(argc,argv,"-calculateSTAR",clarg)) {
    calculate_star = true;
    num_star_calculations = atoi(clarg.c_str());
  }

  voxel_tree = false;
  if (CheckCommandLine(argc,argv,"-voxelTree")) {
    voxel_tree = true;
  }

  print_box_cf_data = false;
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-PrintBoxCfData", clarg)) {
    print_box_cf_data = true;
    box_cf_data_file = clarg;
  }

  box_dir_effect = false;
  if (CheckCommandLine(argc,argv,"-boxDirEffect"))
    box_dir_effect = true;

  no_compartments = false;
  if (ParseCommandLine(argc,argv,"-noCompartments", comp_tree))
    no_compartments = true;

  tree_info = false;
  if (CheckCommandLine(argc,argv,"-treeInfo"))
    tree_info = true;

  clarg.clear();
  constant_star = -1.0;
  if (ParseCommandLine(argc,argv,"-constantSTAR",clarg))
    constant_star = atof(clarg.c_str());

  correct_star = false;
  if (CheckCommandLine(argc,argv,"-correctSTAR"))
    correct_star = true;

  one_time = true;
  if (CheckCommandLine(argc,argv,"-appendMode"))
    one_time = false;

  only_self = false;
  if (CheckCommandLine(argc,argv,"-self"))
    only_self = true;

  clarg.clear();
  many_trees = false;
  if(ParseCommandLine(argc,argv,"-manyTrees",clarg)) {
    many_tree_file = clarg;
    many_trees = true;
    only_self = false;
    ifstream mtf(many_tree_file.c_str());
    if(!mtf){
      cout << "Could not open tree tree information file " <<  many_tree_file << endl;
      exit(0);
    }
    string treef;
    while (mtf.eof() == false)
      {
	getline(mtf,treef);
	treef.erase(std::remove(treef.begin(), treef.end(), ' '),
		    treef.end());                //strip dangerous spaces from file name
	if(mtf.eof() == true)
	  break;
	tree_files.push_back(treef);
      }
    mtf.close();

  }  


  write_only_file = false;
  if (CheckCommandLine(argc,argv,"-writeOnlyFile"))
    write_only_file= true;

  trees_from_file = false;
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-getTreesPos", clarg)) {
    trees_from_file = true;
    trees_pos_file = clarg;
  }

  radius_only = false;
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-radiusOnly", clarg)) {
    radius_only = true;
    max_radius  = atof(clarg.c_str());
  }

  calculateDirectionalStar = false;
  if (CheckCommandLine(argc,argv,"-calculateDirectionalStar"))
    calculateDirectionalStar = true;

  zero_woody_radius = false;
  if(CheckCommandLine(argc,argv,"-zeroWoodyRadius")) {
    zero_woody_radius = true;
  }

  if (verbose){
    cout << "parseCommandLine end" <<endl;

  }   
}  //if(ParseCommand... )

//================================================================================
//Generate tree locations, or read them from a file
//Establish also stand corners with this information
//They are set in StandDescriptor and BorderForest
//================================================================================
template<class TREE, class TS,class BUD>
  void MainProgramAsClass<TREE, TS,BUD>::setTreeLocations()
{  
  if(generate_locations) {
    //In this case plot dimensions are taken from voxelspace x and y side lengths

    //Corners of the stand vs_x, vs_y has been set in ParseCommandLine()
    Point l(0.0, 0.0, 0.0);
    Point r(vs_x, vs_y, 0.0);
    stand.setLlCorner(l);
    stand.setUrCorner(r);
    stand.evaluateArea();
    border_forest.setCornerL(l);
    border_forest.setCornerR(r);

    //ForestGap is here only for consistency with use of GenerateLocations in Lig-Crobas
    //middle_stand has been set in ParseCommandLine() 
    ForestGap gap(pair<double,double>(middle_stand.first,middle_stand.second),gap_radius);

    //number of trees may decrease due to hard core
    int no_trees_0 = no_trees;
    GenerateLocations(no_trees,0.0,0.0,vs_x,vs_y,tree_distance,gap,locations);

    if (verbose){
      cout << "Number of trees" << locations.size() <<endl 
	   << " Density/ha wanted: " << (double)no_trees_0/(vs_x*vs_y/10000.0)
           << " Density/ha created: " << (double)no_trees/(vs_x*vs_y/10000.0) <<endl;
      cout << " Minimum tree distance: " << tree_distance <<endl; 
    }
  } //  if(generate_  ...)
  else {
    ifstream location_stream(location_file.c_str());
    if(!location_stream) {
      cout << "Could not open tree location file " << location_file << endl;
      exit(0);
    }

    string line;
    getline(location_stream,line);
    LGMdouble cx, cy;
    location_stream >> cx >> cy;
    getline(location_stream,line);
    Point l(0.0, 0.0, 0.0);
    Point r(cx, cy, 0.0);
    stand.setLlCorner(l);
    stand.setUrCorner(r);
    stand.evaluateArea();
    border_forest.setCornerL(l);
    border_forest.setCornerR(r);

    getline(location_stream,line);
    no_trees = 0;
    bool stop = false;

    while(!stop) {
      LGMdouble x, y;
      location_stream >> x >> y;
      getline(location_stream,line);
      if(!location_stream.eof()) {
	no_trees++;
	locations.insert(locations.end(), pair<double,double>(x,y));
      }
      else
	stop = true;
    }
    if(no_trees < 1) {
      cout << "Reading tree locations from " << location_file << " did not succeed." << endl;
      exit(0);
    }

    //Defines the middle of the stand rectangle
    middle_stand.first = cx/2.0;
    middle_stand.second = cy/2.0;
  }
}

//==============================================================================================
//Create the  trees.
//===============================================================================================
template<class TREE, class TS,class BUD>
  void MainProgramAsClass<TREE, TS,BUD>::createTrees()
{
  //In the case of many trees & generated (random) positions
  //store positions & trees in the positions to be able to
  //repeat the run with the same configuration
  if(many_trees && generate_locations) {
    ofstream of("runfile.dat", ofstream::trunc);
    of << " x  y  treefile" << endl;
    of.close();
  }

  XMLDomTreeReader<ScotsPineSegment,ScotsPineBud> createTrees_reader;
  for (int i = 0; i < no_trees; i++){
    pair<double,double> p = locations[i];

    TREE* t = new TREE(Point(p.first,p.second,0.0),PositionVector(0,0,1),
		       "sf.fun","fapical.fun","fgo.fun",
		       "fsapwdown.fun","faf.fun","fna.fun", "fwd.fun",
		       "flr.fun");

    if(no_compartments) {
      XMLDomTreeReader<ScotsPineSegment,ScotsPineBud> reader;
      reader.readXMLToTree(*t, comp_tree);

      if(zero_woody_radius) {
	ForEach(*t,ZeroWoodyRadius());
      }

      Compartments result;
      result = Accumulate(*t,result,AddUpCompartments<ScotsPineSegment,ScotsPineBud>());
      cout << "Buds " << result.buds << endl;
      cout << "Segments " << result.segments << endl;
      cout << "Branching points " << result.bps << endl;
      cout << "Axes " << result.axes << endl;

      exit(0);
    }

    //If a voxel space tree (segments are in voxelboxes) is made, above only a empty tree is generated,
    //and the the tree is definedconstructed below in initializeVoxelSpace()


    if(!voxel_tree) {
      if(!many_trees) {
	createTrees_reader.readXMLToTree(*t, input_tree_file);
	if(zero_woody_radius) {
	  ForEach(*t,ZeroWoodyRadius());
	}
      }
      else {
	LGMdouble num_trees = (double)tree_files.size();
	//trees are taken randomly => all trees are used in the same proportion
	unsigned int tree = (unsigned int) (ran3(&ran3_seed)* num_trees);
	createTrees_reader.readXMLToTree(*t, tree_files[tree]);
	if(zero_woody_radius) {
	  ForEach(*t,ZeroWoodyRadius());
	}

	//In the case of many trees & generated (random) positions
	//store positions & trees in the positions to be able to
	//repeat the run with the same configuration
	if(generate_locations) {
	  ofstream of("runfile.dat", ofstream::app);
	  of <<  p.first << " " << p.second << " " << tree_files[tree] << endl;
	  of.close();
	}
      }
   
      //If calculation of STAR is required, it is calculated for the first tree in the tree vector
      //and the program stops

      if(calculate_star) {
	cout << "STAR values for the first tree in the tree vector and then exit" << endl;
	ForEach(*t, CalculateSTAR("STAR.dat",num_star_calculations)); 
	exit(0);
      }

      if(tree_info) {
	LGMdouble treeAf = 0.0;
	treeAf = Accumulate(*t,treeAf,CollectFoliageArea<TS,BUD>());
	DCLData dcld;
	dcld.clear();
	AccumulateDown(*t,dcld,AddBranchWf(),DiameterCrownBase<TS,BUD>());
	LGMdouble d_cb = dcld.DCrownBase();
	LGMdouble h_cb =  dcld.HCrownBase();
	LGMdouble dbh = GetValue(*t,LGADbh);
	LGMdouble h = GetValue(*t,LGAH);
	LGMdouble Wf = 0.0;
	Accumulate(*t,Wf,CollectFoliageMass<TS,BUD>());
	LGMdouble d_base = GetValue(*t,LGADbase);
	LGMdouble dk = 2.0 + 1.25*100.0*dbh;
	LGMdouble cl = h - h_cb;
	LGMdouble Wf_Repola = -2.385 + 15.022*dk/(dk+4.0)-11.979*h/(h+1.0)
	  + 1.116*log(cl)+ 0.5*(0.034 + 0.095);
	Wf_Repola = exp(min(max(-20.0,Wf_Repola),20.0));
	dk = 2.0 + 1.25*100.0*d_base;
	LGMdouble Wf_Repolab = -2.385 + 15.022*dk/(dk+4.0)-11.979*h/(h+1.0)
	  + 1.116*log(cl)+ 0.5*(0.034 + 0.095);
	Wf_Repolab = exp(min(max(-20.0,Wf_Repolab),20.0));
        CrownVolume<TS,BUD> cv;
        LGMdouble cvol = cv(*t);
	LGMdouble NAD = 0.0;
	if(cvol > 0.0)
	  NAD = treeAf/cvol;

	cout << " Tree H (m) Dbh (cm), Dbase (cm) Hcrown_base (m) Dcrown_base (cm)  Wf (kg dm)  Af (m2) NAD(m2/m3) Wf_Repola"
	  " Wf_Repolab"  << endl;
	cout <<  input_tree_file << " "
	     << h << " " << 100.0*dbh << " " << 100.0*d_base << " " << h_cb << " "
	     << 100.0*d_cb << " " << 2.0*Wf << " " << treeAf << " " << NAD <<  " " << Wf_Repola
	     << " " << Wf_Repolab << endl;

	exit(0);
      }

      LGMdouble Af = 0.0;
      Af = Accumulate(*t,Af,CollectFoliageArea<TS,BUD>());
      cout << "Af " << Af << endl;

      Point mov = Point(p.first,p.second,0.0)-GetPoint(*t);
//      MoveTree<TS,BUD> move(Point(p.first,p.second,0.0)-GetPoint(*t),*t);
      MoveTree<TS,BUD> move(mov,*t);

      ForEach(*t, move);

/*       cout << "After " << GetPoint(*t); */

      Axis<TS,BUD>& ax = GetAxis(*t);
/*       TreeSegment<TS,BUD>* ts = GetFirstTreeSegment(ax); */
/*       cout << "my_point " << GetPoint(*ts); */

    if (verbose){
      cout << "Created a tree at: " << p.first << " " << p.second <<endl;
    }

    } //  if(!voxel_tree) ...

    vtree.push_back(t);   //this must here (also for voxel_tree) since in the case of voxel_tree the first
                          //tree in vtree is reshaped as voxeltree in initializeVoxelSpace()
  } //for(int i = ...

}    //::createTrees(


//====================================================================================
//
// initializeVoxelSpace() just creates the voxelspace vs and in the case of voxel_tree reshapes *vtree[0] as
// voxeltree according to the vs

template<class TREE, class TS,class BUD>
  void MainProgramAsClass<TREE, TS,BUD>::initializeVoxelSpace()
{

  // vs_x, vs_y, vs_z, voxboxside have been set in ParseCommandLine()
  vs = new VoxelSpace(Point(0,0,0),Point(vs_x,vs_y,vs_z),
		      voxboxside,voxboxside,voxboxside,
		      static_cast<int>(vs_x/voxboxside)+1,static_cast<int>(vs_y/voxboxside)+1,
		      static_cast<int>(vs_z/voxboxside)+1,GetFirmament(*vtree[0]));

  if(voxel_tree) { //Now the structure of the tree in voxel space can be specified
    cout << endl << "Radiation calculations now with voxel tree in voxel space" << endl;

    TREE* t = vtree[0];       //Note absolutely only one tree in case of voxel_tree
    Axis<TS,BUD>& axis = GetAxis(*t);


    int nx = vs->Xn;
    int ny = vs->Yn;
    int nz = vs->Zn;

    LGMdouble deltaZ = vs->getZSideLength();
 
    double fol_den = 150.0;
    double L = 0.15;
    double Rw = 0.000001;
    double Rf = 0.025;
    double sf = 24.0;
    double Wf = PI_VALUE*Rf*Rf*L*fol_den/sf;
    double Rh = 0.0;
    Point ms(middle_stand.first, middle_stand.second, 0.0);
    PositionVector up(0.0,0.0,1.0);

    cout << " Tree segment: L Rf Wf Af " << L << " " << Rf << " " << Wf << " " << Wf*sf << endl;

    MoveTree<TS,BUD> move(ms-GetPoint(*t),*t);
    ForEach(*t, move);

    Point bud_p  = ms + Point(0.0,0.0, vs->getUpperRightCorner().getZ());
    BUD* bud = new BUD(bud_p, up, 1.0, t); 
    InsertTreeCompartment(axis, bud);
    
    TS* ts;
    BranchingPoint<TS,BUD>* bp;

    for(int i = 0; i < nx; i++)
      for(int j = 0; j < ny; j++)
	for(int k = 0; k < nz; k++) {
	  Point p = vs->voxboxes[i][j][k].getCenterPoint();
	  Point loc = p - Point(0.0,0.0,deltaZ/2.0-0.0001); //a bit (0.0001) above the bottom
	  ts = new TS(loc, up, 1.0, L, Rw, Rh, t);
	  SetValue(*ts,LGAR,Rw);
	  SetValue(*ts,LGARh,Rh);
	  SetValue(*ts,LGAL,L);
	  SetValue(*ts,LGAsf,sf);
	  SetValue(*ts,LGAWf,Wf);
	  SetValue(*ts,LGARf,Rf);
	    
	  bp = new BranchingPoint<ScotsPineSegment,ScotsPineBud>(loc,up,1.0,t);
	  InsertTreeCompartmentSecondLast(axis, ts);
	  InsertTreeCompartmentSecondLast(axis, bp);
	}

    //Write the voxeltree to a file
    XMLDomTreeWriter<TS,BUD> writer;
    writer.writeTreeToXML(*t,"voxel-tree.xml");

    //Dump voxeltree here into voxelspac and NOT in setVoxelSpaceAndBorderForest()
    DumpCfTree(*vs, *t, num_parts, wood_voxel);
    //After dumping all trees to VoxelSpace it is necessary to evaluate
    //sum quantities (e.g. STAR_mean)
    vs->updateBoxValues();

    LGMdouble fa = 0.0;
    fa = Accumulate(*t,fa,CollectFoliageArea<TS,BUD>());
    cout << "Foliage area of the voxel tree " << fa << "  m2" << endl;

  }  //if(voxel_tree) ....

} // end of initializeVoxelSpace() 

template<class TREE, class TS,class BUD>
void MainProgramAsClass<TREE, TS,BUD>::initializeFunctions()
{  
  K.install("K.fun");
  if (verbose){
    cout << "K() O.K: " << K.ok() << endl; //" stems_ha() O.K: " << stems_ha.ok() 
  }
}


//===================================================================
// Prepare for radiation calculations:
// Set BoundingBox and VoxelSpace for trees in the stand.
// Dump also the foliage (except the first tree) into voxels and
// set BorderForest
//===================================================================
template<class TREE, class TS,class BUD>
  void MainProgramAsClass<TREE, TS,BUD>::setVoxelSpaceAndBorderForest()
{
  //Check first if -Radius only is set and remove from tree vector trees
  //that are further away than max_radius.

  if(radius_only) {
    typename vector<TREE*>::iterator I;
    for(I = vtree.begin(); I != vtree.end(); I++) {
      if(sqrt(pow(GetPoint(**I).getX()-middle_stand.first,2.0)+
	      pow(GetPoint(**I).getY()-middle_stand.second,2.0)) > max_radius) {
	delete *I;
	vtree.erase(I);
	no_trees--;
	I--;
      }
    }
    if(no_trees == 0) {
      cout << endl << "Radius " << max_radius << "  too small, no trees in the circle." << endl;
      exit(0);
    }
  }
    
  cout << endl << "no_trees " << no_trees << endl;

  BoundingBox bb;
  FindCfBoundingBox<TS,BUD> fb;      //Bounding box including the segments that have no foliage

  for (unsigned int k = 0; k < (unsigned int)no_trees; k++) {
    bb = Accumulate(*vtree[k], bb, fb);
  }

  Point ll = bb.getMin();
  Point ur = bb.getMax();
  cout << ll;
  cout << ur;

  if(!voxel_tree) {   //In the case of voxeltree it is dumped already in initializeVoxelSpace()
    vs->resize(ll, ur);

    vs->reset();          //this resets all data in voxelboxes before dumping

    for (unsigned int k = 0; k < (unsigned int)no_trees; k++) {
      DumpCfTree(*vs, *vtree[k], num_parts, true);
    }
    vs->updateBoxValues(); 

  }  // end of if(!voxel_tree) ...

  //Now the shading tree information is in the voxelspace
  //write out information about it
  cout << endl << "Information about voxelspace and dumped shading trees: " << endl;
  cout << "getArea() " << vs->getArea() 
       << " getLowerLeftCorner() " << vs->getLowerLeftCorner();
  cout << "getUpperRightCorner() " << vs->getUpperRightCorner();
  cout << "getNumberOfBoxes() " << vs->getNumberOfBoxes() << " getNumberOfFilledBoxes() "
       << vs->getNumberOfFilledBoxes() << " getNumberOfTreeSegments() "
       << vs->getNumberOfTreeSegments() << endl;
  cout << "getBoxVolume() " << vs->getBoxVolume() << " getXSideLength() " << vs->getXSideLength() << endl;
  cout << "getYSideLength() " << vs->getYSideLength() << " getZSideLength() " << vs->getZSideLength() << endl;
  cout << "getNoBoxX() " << vs->getNoBoxX() << " getNoBoxY() " << vs->getNoBoxY()
       << " getNoBoxZ() " << vs->getNoBoxZ() << endl;
  cout << "getNeedleArea() " << vs->getNeedleArea() << " getLeafArea() " <<  vs->getLeafArea() << endl;
  cout << endl;

  if(evaluate_LAI) {
    //Lasketaan LAI VoxelSpacen avulla
    //Tehdaan oma VoxelSpace sita varten
    //Target puu (alla: calculateRadiation()) tehdaan niin, etta segmentit ovat
    //tasoilla Hcb, ..., H, yhteensa 10 kpl. Kun tehdaan voxelboxeja 9 kpl
    //valilla [Hcb,H] niin vastaa target puun jakoa: 1. voxelboxtaso  varjostaa 2. anturia,
    //9. taso anturia 10.
    //Dumpataan kaikki puut voxelspaceen: target puu asetetaan sinne vasta myohemmin
    //(calculateRadiation()).
    //Ei ole niin valia, kuinka monta boxia x ja y suunnassa on (voisi olla 1).

    Point lai_ll = vs->getLowerLeftCorner();
    Point lai_ur = vs->getUpperRightCorner();

    VoxelSpace lai_vs(lai_ll, lai_ur, 10, 10, 9, GetFirmament(*vtree[0]));
    for (unsigned int k = 0; k < (unsigned int)no_trees; k++)
      DumpCfTree(lai_vs, *vtree[k], num_parts, true);

    vector<pair<LGMdouble,LGMdouble> > lai_h;
    LGMdouble Hmin, Hmax;
    int n_levels;
    lai_vs.evaluateVerticalNeedleAreaDensity(Hmax, Hmin, n_levels, lai_h);
    cout << endl << "This is output about LAI etc calculated in a LAI voxelspace: " << endl;
    cout << "Needle area in space: " << lai_vs.getFoliageArea() << endl;
    cout << "Hmin " << Hmin << "  Hmax " << Hmax << endl;
    cout << "left corner (" << lai_ll.getX() << " , " << lai_ll.getY() << " )     right corner  ( " 
	 << lai_ur.getX() << " , " << lai_ur.getY() << " )" << endl;   
    cout << "Box_center    lai cum_lai" << endl;
    LGMdouble cum_lai = 0.0;
    for(int i = 0; i < n_levels; i++) {
      cum_lai += lai_h[i].second;
      cout << lai_h[i].first << " " << lai_h[i].second << " " << cum_lai << endl;
    }
    cout << "LAI according to voxelspace: " << lai_vs.getFoliageArea()/ ((lai_ur.getX()-lai_ll.getX())*
									 (lai_ur.getY()-lai_ll.getY())) << endl;

    cout << endl << "This voxel space for LAI:" << endl << endl;
    cout << "getArea() " << lai_vs.getArea() 
	 << "getLowerLeftCorner() " << lai_vs.getLowerLeftCorner();
    cout << "getUpperRightCorner() " << lai_vs.getUpperRightCorner();
    cout << " getNumberOfBoxes() " << lai_vs.getNumberOfBoxes() << " getNumberOfFilledBoxes() "
	 << lai_vs.getNumberOfFilledBoxes() << " getNumberOfTreeSegments() "
	 << lai_vs.getNumberOfTreeSegments() << endl;
    cout << " getBoxVolume() " << lai_vs.getBoxVolume() << " getXSideLength() "
	 << lai_vs.getXSideLength() << endl;
    cout << "getYSideLength() " << lai_vs.getYSideLength() << " getZSideLength() "
	 << lai_vs.getZSideLength() << endl;
    cout << "getNoBoxX() " << lai_vs.getNoBoxX() << " getNoBoxY() " << lai_vs.getNoBoxY()
	 << "getNoBoxZ() " << lai_vs.getNoBoxX() << endl;
    cout << "getNeedleArea() " << lai_vs.getNeedleArea() << " getLeafArea() "
	 <<  lai_vs.getLeafArea() << endl;
    cout << endl;

    exit(0);
  } // end of if(evaluate_LAI) ..

  //The border forest top is set according to voxelspace
  // (i.e. bounding box)
  // so that the dimesions match. The dimensions of top height of the stand could
  // also be used but it is not exactly the same as of bounding box
  // since it considers also needles: bounding box is higher
  // than stand top height by length of needles. 
  border_forest.setH(bb.getMax().getZ());
  border_forest.setHcb(stand.getMinCrownLimit());
  border_forest.setLAI(stand.getLAI());

}  // end of setVoxelSpaceAndBorderForest() { ...

//  calculateDirectionalStar = false;
//  if (CheckCommandLine(argc,argv,"-calculateDirectionalStar"))
//    calculateDirectionalStar = true;


//===================================================================
// calculateRadiation 
//===================================================================
template<class TREE, class TS,class BUD>
  void MainProgramAsClass<TREE, TS,BUD>::calculateRadiation()
{

  //Find max height (Hmax) and minimum height of crown base (Hcbmin)
  //and construct then the target tree

  BoundingBox bb;
  //  FindCfBoundingBox<TS,BUD> fb(true);
  FindCfBoundingBox<TS,BUD> fb(false);

  for (unsigned int k = 0; k < (unsigned int)no_trees; k++) {
    bb = Accumulate(*vtree[k], bb, fb);
  }

  cout << "ur " << bb.getMax();
  cout << "ll " << bb.getMin();


  LGMdouble Hmax = bb.getMax().getZ()  - 0.001;  //Stay 0.001 inside voxelspace
  LGMdouble Hcb = bb.getMin().getZ()   + 0.001;

  vector<double> heights(10);
  vector<double> fii(13);
  fii[0] = 0.0;
  for(int i = 1; i < 6; i++)
    fii[i] = fii[i-1] + PI_VALUE/3.0;
  fii[6] = PI_VALUE/6.0;
  for(int i = 7; i < 12; i++)
    fii[i] = fii[i-1] + PI_VALUE/3.0;
  fii[12] = 0.0;

  vector<double> r(13);
  r[0] = r[1] = r[2] = r[3] = r[4] = r[5] = 1.0;
  r[6] = r[7] = r[8] = r[9] = r[10] = r[11] = 0.618;
  r[12] = 0.0;

  Point mp = Point(middle_stand.first,middle_stand.second,0.0);
  PositionVector up(0.0,0.0,1.0);

  TREE* t_t = new TREE(mp,up,"sf.fun","fapical.fun","fgo.fun",
		       "fsapwdown.fun","faf.fun","fna.fun", "fwd.fun",
		       "flr.fun");   //target_tree

  //GetFirmament(*t_t).resize(15,15,1200.0);

  Axis<ScotsPineSegment,ScotsPineBud>& t_t_a = GetAxis(*t_t);
    
  ScotsPineBud* bud = new ScotsPineBud(Point(mp.getX(),mp.getY(),Hmax), up, 1.0, t_t); 
  InsertTreeCompartment(t_t_a, bud);

  double fol_den = 100.0;
  double L = 0.01;
  double Rw = 1.0e-20;//0.000001;
  double Rf = 0.01;
  double sf = 14.0;
  double Wf = PI_VALUE*Rf*Rf*L*fol_den/sf;
  double Rh = 0.0;

  for(int i = 0; i < 10; i++) {
    LGMdouble z = Hcb + ((double)i)*(Hmax-Hcb)/9.0 - L/2.0;
    for(int j = 0; j < 13; j++) {
      LGMdouble x = target_tree_rad*r[j]*cos(fii[j]) + mp.getX();
      LGMdouble y = target_tree_rad*r[j]*sin(fii[j]) + mp.getY();
      ScotsPineSegment* ts = new ScotsPineSegment(Point(x,y,z),up,1.0,L,Rw,Rh,t_t);
      SetValue(*ts,LGAsf,sf);
      SetValue(*ts,LGAWf,Wf);
      SetValue(*ts,LGARf,Rf);
      SetValue(*ts,LGAR,Rw); 
      SetValue(*ts,LGARh,Rh);

      BranchingPoint<ScotsPineSegment,ScotsPineBud>* bp =
	new BranchingPoint<ScotsPineSegment,ScotsPineBud>(Point(x,y,z),up,1.0,t_t);
    
      InsertTreeCompartmentSecondLast(t_t_a, ts);
      InsertTreeCompartmentSecondLast(t_t_a, bp);
    }
  }

  //Finally set the target tree as the first tree in the vtree vector 
  vtree.insert(vtree.begin(),t_t);
  no_trees++;

  cout << endl << "Target tree at " << mp;
  cout << "Between heights " << Hmax << " and " << Hcb << endl;
  cout << "There are 10 levels with 13 ""sensors"" at each level" << endl;


  //Tama taytyy jatkossa paremmin mutta simulointien perusteella nayttaa, etta boxissa
  //k/STAR_mean = 0.97 + segment_length/Box_edge. Kaytetaan segment_length = 0.1 m. Josta
  //k = (0.97 + (0.1/3-juuri(box_vol)))*0.14
  //Taytyy parantaa jatkossa
  //  LGMdouble green_ext = (0.97 + 2.0* 0.1/pow(vs->getBoxVolume(), 0.33333))*0.14;


  LGMdouble a, b;
  ifstream ab_file("box-radiation.txt");
   char dummy[256];
   ab_file.getline(dummy,256);    //header
   ab_file >> a >> b;
   ab_file.close();

   cout << endl;
   cout << " a  b " << a << " " << b << endl;


   //There need to me more trees than only target tree
  if(no_trees <= 1) {
    cout << "Only target tree or nothing in calculateRadiation() - exiting." << endl;
    exit(0);
  }

   if(print_box_cf_data) {
     PrintBoxCfData(*vs, box_cf_data_file, false);
    exit(0);
   }


  //EvaluateRadiationForCfTreeSegment_1<ScotsPineSegment,ScotsPineBud> Rad(K,vs,&border_forest, green_ext);
  //Alla on pairwise kaikille puille
   K = ParametricCurve("K.fun");
   EvaluateRadiationForCfTreeSegment_2<ScotsPineSegment,ScotsPineBud, TREE> Rad2(K,vtree);
  //Tassa pairwise & voxelit
  //EvaluateRadiationForCfTreeSegment_1s<ScotsPineSegment,ScotsPineBud> Rad(K,vs,&border_forest, green_ext);

  //Tassa vain voxel
   bool virittely_dump = false;                    //HUOM virittely



   //Before radiation calculations set up the sky
 for (unsigned int k = 0; k < (unsigned int)no_trees; k++) {
   GetFirmament(*vtree[k]).resize(15,15,1200.0);
 }

 SetValue(*vtree[0],TreeQinMax,GetFirmament(*t_t).diffuseBallSensor());    //target tree is the first



   EvaluateRadiationForCfTreeSegment_3<ScotsPineSegment,ScotsPineBud>
     Rad3(K, vs, &border_forest, false, a, b, virittely_dump, k_border_conifer,
         box_dir_effect, wood_voxel, correct_star, constant_star,calculateDirectionalStar);


   SetStarMean<TS,BUD> setstar(ParametricCurve(0.14));
   ResetQinQabs<TS,BUD> RQQ;


   for (unsigned int k = 0; k < (unsigned int)no_trees; k++) {
     TREE* t = vtree[k];
     ForEach(*t,setstar);
     ForEach(*t,RQQ);
    
     if(rad_method == 2) {
       cout << "Radiation method _2" << endl;
       ForEach(*t,Rad2);
     }
     else if(rad_method == 3) {
       cout << "Radiation method _3" << endl;
       ForEach(*t,Rad3);
     }
     else {
       cout << "No radiation calculation method." << endl;
       exit(0);
     }

     //Tassa tulostetaan vtree vektorin 1. puun tiedot. Se on target puu

  if(one_time) {
   ResultsToFile rtf(resultfile);
    rtf.setQmax(GetFirmament(*t).diffuseBallSensor());

     ForEach(*t, rtf);
  }
  else {
    //int bs;                 // box STAR  0 = constant value,
                          //           1 = value of box
                          //           2 = value of box & mean direction effect
    //int corr;               // If STAR_eq -> STAR correction considered (0 = no, 1 = yes)
    //LGMdouble vox;          // Size of voxel (=length of side, assuming voxels are cubes)

    ResultsToFileInfo info;
    if(correct_star)
      info.corr = 1;
    else
      info.corr = 0;

    if(constant_star > 0.0) {
      info.bs = 0;
      info.corr = 0;
    }
    else {
      if(box_dir_effect)
	info.bs = 2;
      else
	info.bs = 1;
    }

    info.vox = voxboxside;
    info.location_file = location_file;
    info.tree_file = input_tree_file;

    ResultsToFile rtf(resultfile, info);
    rtf.setQmax(GetFirmament(*t).diffuseBallSensor());

    ForEach(*t, rtf);
  }
 
     exit(0);

   }
}  //end of calculateRadiation()  { ..





//===================================================================
// calculateRadiationOnlySelf()
//===================================================================

template<class TREE, class TS,class BUD>
  void MainProgramAsClass<TREE, TS,BUD>::calculateRadiationOnlySelf(){

  //1) Create the tree
  TREE* t = new TREE(Point(10.0,10.0,0.0),PositionVector(0,0,1),
		     "sf.fun","fapical.fun","fgo.fun",
		     "fsapwdown.fun","faf.fun","fna.fun", "fwd.fun",
		     "flr.fun");

  XMLDomTreeReader<ScotsPineSegment,ScotsPineBud> reader;
  reader.readXMLToTree(*t, input_tree_file);
  if(zero_woody_radius) {
    ForEach(*t,ZeroWoodyRadius());
  }

  //2) First and only tree in vector vtree

  vtree.insert(vtree.begin(),t);
  no_trees = 1;

  MoveTree<TS,BUD> move(Point(10.0,10.0,0.0)-GetPoint(*t),*t);
  ForEach(*t, move);

  //3) Radiation Calclulation

  if(rad_method == 2) {
    K = ParametricCurve("K.fun");
    EvaluateRadiationForCfTreeSegment_2<ScotsPineSegment,
      ScotsPineBud, TREE> Rad2(K,vtree, only_self);
    TREE* t = vtree[0];
    ForEach(*t,Rad2);
  }
  else {
    TREE* t = vtree[0];

  // vs_x, vs_y, vs_z, voxboxside have been set in ParseCommandLine()
    vs = new VoxelSpace(Point(0,0,0),Point(vs_x,vs_y,vs_z),
			voxboxside,voxboxside,voxboxside,
			static_cast<int>(vs_x/voxboxside),static_cast<int>(vs_y/voxboxside),static_cast<int>(vs_z/voxboxside),
			GetFirmament(*vtree[0]));

    BoundingBox bb;
    //  FindCfBoundingBox<TS,BUD> fb(!wood_voxel);
    FindCfBoundingBox<TS,BUD> fb(false);
    bb = Accumulate(*t, bb, fb);
    Point ll = bb.getMin();
    Point ur = bb.getMax();

    vs->resize(ll, ur);

    DumpCfTree(*vs, *t, num_parts, true);

    bool only_self_dump = true;   // is dump_self
    K = ParametricCurve("K.fun");
    double a = 0.0, b = 0.0;

    EvaluateRadiationForCfTreeSegment_3<ScotsPineSegment,ScotsPineBud>
      Rad3(K, vs, &border_forest, false/*border forest*/, a, b, only_self_dump, k_border_conifer,
       box_dir_effect, wood_voxel, correct_star, constant_star,calculateDirectionalStar);
    ForEach(*t,Rad3);
  }

  // 4) Output

  if(one_time) {
    ForEach(*t, ResultsToFile(resultfile));
  }
  else {
    //int bs;                 // box STAR  0 = constant value,
    //           1 = value of box
    //           2 = value of box & mean direction effect
    //int corr;               // If STAR_eq -> STAR correction considered (0 = no, 1 = yes)
    //LGMdouble vox;          // Size of voxel (=length of side, assuming voxels are cubes)

    ResultsToFileInfo info;
    if(correct_star)
      info.corr = 1;
    else
      info.corr = 0;

    if(constant_star > 0.0) {
      info.bs = 0;
      info.corr = 0;
    }
    else {
      if(box_dir_effect)
	info.bs = 2;
      else
	info.bs = 1;
    }

    info.vox = voxboxside;
    info.location_file = location_file;
    info.tree_file = input_tree_file;
    
    ResultsToFile rtf(resultfile, info);
    rtf.setQmax(GetFirmament(*t).diffuseBallSensor());
    ForEach(*t, rtf);
  }
 
}

/* #undef HIT_THE_FOLIAGE */
/* #undef NO_HIT */
/* #undef HIT_THE_WOOD */

class SetSf {
 public:
  SetSf() {multiplier = 1.0;}
  SetSf(LGMdouble mu): multiplier(mu) {}
  TreeCompartment<ScotsPineSegment,ScotsPineBud>*
    operator()(TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc) const
    {
      if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
	LGMdouble h = GetPoint(*ts).getZ();
	LGMdouble L = 9.3/(1.0+exp(5.0*(h-1.8)));
	//This is from Palmrooth & Hari in Palmrooth's PhD thesis
	LGMdouble sf = multiplier*(24.0 + (L/10.0)*(34.0 - 24.0));
	SetValue(*ts, LGAsf, sf);
      }
      return tc;
    }
 private:
  LGMdouble multiplier;
};
class SetSf2 {
 public:
  SetSf2() {multiplier = 1.0;}
  SetSf2(LGMdouble mu): multiplier(mu) {}
  TreeCompartment<ScotsPineSegment,ScotsPineBud>*
    operator()(TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc) const
    {
      if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
	LGMdouble h = GetPoint(*ts).getZ();
	LGMdouble L = 12.4/(1.0+exp(4.5*(h*-1.8)));            
	//This is from Palmrooth & Hari in Palmrooth's PhD thesis
	LGMdouble sf = multiplier*(24.0 + (L/12.4)*(34.0 - 24.0));
	SetValue(*ts, LGAsf, sf);
      }
      return tc;
    }
 private:
  LGMdouble multiplier;
};


//Calculates STAR values of segments and writes to a file. 
//Argument nruns in the constructor specifies the number of trials (beams of
//radiation shot through the segment randomly) for the accurate STAR calculation.
//If it is high, it takes a long time to calculate.

class CalculateSTAR {
 public:
  CalculateSTAR(const string& filename, int nruns): fname(filename), runs(nruns)
    {
      ofstream f(fname.c_str() , ofstream::trunc);
      f << "Height Dist_top Dist_stem age Rf L Af Wf Vf fol_den STAR STAR_eq" << endl;
      f.close();
    }

    TreeCompartment<ScotsPineSegment,ScotsPineBud>*
      operator()(TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc) const
      {
        if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
          LGMdouble Wf = GetValue(*ts,LGAWf);
          if(Wf > 0.0) {
            LGMdouble Rf = GetValue(*ts,LGARf);
            LGMdouble L = GetValue(*ts,LGAL);
            LGMdouble Age = GetValue(*ts, LGAage);
            LGMdouble Af = GetValue(*ts, LGAAf);
            Tree<ScotsPineSegment,ScotsPineBud>& t = GetTree(*ts);
            LGMdouble H = GetValue(t, LGAH);
            Point pp = GetPoint(*ts);
            LGMdouble ver_dist = H - pp.getZ();
            Axis<ScotsPineSegment,ScotsPineBud>& ax = GetAxis(t);
            TreeCompartment<ScotsPineSegment,ScotsPineBud>* fc = GetFirstTreeCompartment(ax);
            LGMdouble hor_dist = sqrt(pow(pp.getX()-GetPoint(*fc).getX(),2.0) +
                                      pow(pp.getY()-GetPoint(*fc).getY(),2.0));
	    LGMdouble Vf = GetValue(*ts,LGAVf);
	    LGMdouble fol_den = Af/Vf;
	    LGMdouble star = star_mean(runs,Rf,L,fol_den, ran3_seed);
	    LGMdouble s_sum = 0.0, c_sum = 0.0;
	    LGMdouble Sf = GetValue(*ts,LGAsf);
	    LGMdouble Wf = GetValue(*ts, LGAWf);
	    for(int k = 0; k < 11 ; k++) {
	      double theta = PI_VALUE * (double)k /(2.0 * 10.0);
	      s_sum += cos(theta)*S(theta,Sf,Wf,Rf,L);
	      c_sum += cos(theta);
	    }
	    LGMdouble star_eq = s_sum/c_sum;

        
            ofstream f(fname.c_str() , ofstream::app);
            f << pp.getZ() << " " << ver_dist << " " << hor_dist << " " << Age <<  " " << 100.0*Rf
	      << " " << 100.0*L << " " << 10000.0*Af << " " << 2000.0*Wf << " " 
	      << 1000000.0*Vf << " " << fol_den << " " << star << " " << star_eq << endl;
            f.close();
          }
        }
        return tc;
      }

 private:
    string fname;
    int runs;
};


template<class TREE, class TS, class BUD>
  void MainProgramAsClass<TREE,TS,BUD>::getTreesAndPositions() {  
  ifstream tpf(trees_pos_file.c_str());
  if(!tpf) {
    cout << "Could not open tree and position input file " << trees_pos_file << endl;
    exit(0);
  }

  string line;
  getline(tpf,line); //First line
  cout << line << endl;

  no_trees= 0;
  while (tpf.eof() == false)
    {
      getline(tpf,line);

      if(tpf.eof() == true)
        break;

      LGMdouble x, y;
      string xml_file;
      istringstream ist(line);
      ist >> x >> y >> xml_file;
      
      XMLDomTreeReader<ScotsPineSegment,ScotsPineBud> reader;

      TREE* t = new TREE(Point(0.0,0.0,0.0),PositionVector(0,0,1),
  		       "sf.fun","fapical.fun","fgo.fun",
  		       "fsapwdown.fun","faf.fun","fna.fun", "fwd.fun",
  		       "flr.fun");

      reader.readXMLToTree(*t, xml_file);
      if(zero_woody_radius) {
	ForEach(*t,ZeroWoodyRadius());
      }

      MoveTree<ScotsPineSegment,ScotsPineBud>
	move(Point(x,y,0.0)-GetPoint(*t),*t);
      ForEach(*t, move);

      vtree.push_back(t);
	
      no_trees++;

      if (verbose){
	cout << "Tree " << xml_file << " at (" << x << ", "
	     << y << ")" << endl;
      }
    }
}

#endif
