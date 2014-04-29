#ifndef GROWTHLOOPRADIATIONI_H
#define GROWTHLOOPRADIATIONI_H

class SetSf;  //This is defined at the end of this file (declared here before use)
class SetSf2;
class CalculateSTAR;

extern int ran3_seed;
extern bool no_compartments;
extern string comp_tree;

//This is GrowthLoopI.h for radiation calculations. See LignumRadiation.pro for information
// about LignumRadiation project.
// GrowthLoopRadiationI.h differs from GrowthLoopI.h in that
// 1) there are a few modifications for radiation calculations, marked by comments
//  virittelya ....
// and 2) functors  EvaluateRadiationForCfTreeSegment_1s EvaluateRadiationForCfTreeSegment_2
// are in this file instead of include/CalculateLight.h.
// They presumably are not used in growth simulations.

template<class TREE, class TS, class BUD>
GrowthLoop<TREE,TS,BUD>::~GrowthLoop()
{
  for (unsigned int i = 0; i < vtree.size(); i++){
    delete vtree[i];
  }

/*   for (unsigned int i = 0; i < vlsystem.size(); i++){ */
/*     delete vlsystem[i]; */
/*   } */

  for (unsigned int i = 0; i < vdatafile.size(); i++){
    vdatafile[i]->close();
    delete vdatafile[i];
  }
  
  //  stand_output->close();
  // delete stand_output;
  
}

template<class TREE, class TS, class BUD>
  void GrowthLoop<TREE,TS,BUD>::usage()const
{
  cout << "Usage:  ./lig-forest -iter <value>  -metafile <file>  -voxelspace <file>" <<endl;
  cout << "[-numParts <parts>]  [-treeDist <dist>] [-hw <hw_start>] [-viz]" <<endl;
  cout << "[-xml <filename>] [-writeVoxels]" <<endl;
  cout << "[-treeFile <filename>] [-generateLocations  <num>] [-woodVoxel] [-treeLocations <file>]" << endl;
  cout << "[-phprodfile <file>] [-Voxbox <value>]" << endl;
  cout << "[-dumpSelf] [-inputTree <filename>] [-kBorderConifer <value>] [-GapRadius <value>]" << endl;
  cout << "[-targetTreeRad <value>] [-evaluateLAI] [-radMethod <num>] [-calculateSTAR <num>] [-calculateDirectionalStar] " << endl;
  cout << "[-voxelTree] [-boxDirEffect] [-treeInfo] [-segmentInfo <file>]" << endl;
  cout << "[-correctSTAR] [-constantSTAR <value>] [-appendMode] [-self] [-manyTrees <file>]" << endl;
  cout << "[-writeOnlyFile] [-getTreesPos <file>] [-radiusOnly <m>]" << endl;
   cout << endl;
  cout << "-generateLocations <num>  In this case <num> trees will be generated to random locations. If this" << endl;
  cout << "          is not on, tree locations will be read from file Treelocations.txt. This file can be changed" << endl;
  cout << "          by -treeLocations <file>. If location file is not found program stops." << endl;
  cout << "-woodVoxel                If woody parts are dumped to voxels (default = true)" << endl;
  cout << "-calculateDirectionalStar If directional star needs to be calculated then use true (default = false)  "<<endl;
  cout << "-treeDist <dist>          Minimum distance between two trees (default = 0), works only with -generateLocations." << endl;  
  cout << "-numParts <parts>         Segments can be dumped to voxels by parts (i.e. they may belong to different voxels," << endl;
  cout << "-targetTree <num>         Any one of the trees can be identified as target tree (default = 0)" << endl;
  cout << "-phprodfile <file>        File to store radiation calculations"  << endl;
  cout << "-Voxbox <value>        Side length of voxel box - overrides the one given in VoxelSpace.txt"  << endl;
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
  cout << "-voxelTree             In this case a tree is made which has one segment in each voxel (exept in the volume of the"
          "                       target tree). NOTE Treelocations file or -generate MUST be such that the is only one tree."
       << endl;
  cout << "-PrintBoxCfData <file> Writes voxelboxcontents to file with VoxelSpace.PrintBoxCfData(). Then stops"
       << endl;
  cout << "-boxDirEffect          If effect of mean direction of segments in box considered (default = no)"
       << endl;
  cout << "-treeInfo              Writes H, Dbh, Dbase, Hcrown_base Dcrown_base Wf Af of tree on console and stops." << endl;
  cout << "-correctSTAR           If the descrepancy with STATR from eq. and STAR (ca. STAR = -0.041 + 1.056*STAR_eq) is corr."
       << endl;
  cout << "-constantSTAR <value>  STAR has constant value <value> (may be corrected by -correctSTAR)."  << endl;
  cout << "-appendMode            If production balance is written to -phprodfile in trunc-mode or app-mode." << endl;
  cout << "                       trunc mode is default, in app-mode information about voxels is written also." << endl;
  cout << "-self                  Calculates the radiation conditions only for -inputTree" << endl;
  cout << "-manyTrees <file>      Many shading trees, they are given in <file>. In this case -self has no effect" << endl; 
  cout << "-writeOnlyFile         Writes only positions & trees to runfile.dat and exits, requires -manyTrees" << endl;
  cout << "-getTreesPos <file>    Reads trees (xml files) and their positions from file" << endl;
  cout << "-radiusOnly <r>        Only trees at max distance r m are used in the calculation." << endl;
  cout  << endl;

}


template<class TREE, class TS, class BUD>
void GrowthLoop<TREE,TS,BUD>::checkCommandLine(int argc, char** argv)const
{
  //At least three  mandatory arguments required 
  if (argc < 4){
    cout << "Three mandatory arguments are required!" << endl << endl;
    usage();
    exit(0);
  }
  else if (CheckCommandLine(argc,argv,"-iter") == false){
    cout << "Mandatory -iter <num> option missing" << endl;
    exit(0);
  }
  else if (CheckCommandLine(argc,argv,"-metafile") == false){
    cout << "Mandatory -metafile <MetaFile.txt> option missing" << endl;
    exit(0);
  }
  else if (CheckCommandLine(argc,argv,"-voxelspace") == false){
    cout << "Mandatory -voxelspace <VoxelSpace.txt> option missing" << endl;
    exit(0);
  }
  else if (verbose){
    cout << "Command line O.K." <<endl;
  } 
}

template<class TREE, class TS, class BUD>
void GrowthLoop<TREE,TS,BUD>::parseCommandLine(int argc, char** argv)
{
  if (verbose){
    cout << "parseCommandLine begin" <<endl;
  }

  checkCommandLine(argc,argv);

  //Mandatory arguments
  string clarg;
  if (ParseCommandLine(argc,argv,"-iter", clarg)){
    iterations = atoi(clarg.c_str());
  }

  clarg.clear();
  if (ParseCommandLine(argc,argv,"-metafile", clarg)){
    metafile = clarg;
  }
  //Read here and set middle_stand (may be aread and set also in other
  //places
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-voxelspace", clarg)){
    voxelfile = clarg;
    ifstream vf(voxelfile.c_str());
    int vx,vy,vz; vx = vy = vz = 0;
    LGMdouble s1,s2,s3,b1,b2; s1 = s2 = s3 = 0.0; b1=b2=0.0;
    vf >> vx >> vy >> vz >> s1 >> s2 >> s3 >> b1 >> b2;
    vf.close();

    //Tama valonlaskentaa varten ja maarittelee standin keskipisteen
    middle_stand.first = vx/2.0;
    middle_stand.second = vy/2.0;
  }

  //End of mandatory arguments 

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

   phprodfile = "phprod.dat";
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-phprodfile", clarg))
      phprodfile = clarg;

  voxboxside = -1.0;
  clarg.clear();
  if (ParseCommandLine(argc,argv,"-Voxbox", clarg))
       voxboxside = atof(clarg.c_str());

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
	if(mtf.eof() == true)
	  break;
	tree_files.push_back(treef);
      }
    mtf.close();

  }  //if(ParseCommand... )


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
  //cout<<"THS IS THE VALE "<<bool(calculateDirectionalStar)<<endl;



  if (verbose){
    cout << "parseCommandLine end" <<endl;

  }   
}

//================================================================================
//Generate tree locations, or read them from a file
//Establish also stand corners with this information
//They are set in StandDescriptor and BorderForest
//================================================================================
template<class TREE, class TS,class BUD>
  void GrowthLoop<TREE, TS,BUD>::setTreeLocations()
{  
  if(generate_locations) {
    //In this case the take the plot dimensions from voxelspace file: Read the voxel space file
    ifstream vf(voxelfile.c_str());
    int vx,vy,vz; vx = vy = vz = 0;
    LGMdouble s1,s2,s3,b1,b2; s1 = s2 = s3 = 0.0; b1=b2=0.0;
    vf >> vx >> vy >> vz >> s1 >> s2 >> s3 >> b1 >> b2;
    
    //Corners of the stand
    Point l(0.0, 0.0, 0.0);
    Point r(vx, vy, 0.0);
    stand.setLlCorner(l);
    stand.setUrCorner(r);
    stand.evaluateArea();
    border_forest.setCornerL(l);
    border_forest.setCornerR(r);


    //Defines the middle of the stand rectangle
    middle_stand.first = vx/2.0;
    middle_stand.second = vy/2.0;

    //ForestGap is here only for consistency with use of GenerateLocations in Lig-Crobas
    ForestGap gap(pair<double,double>(middle_stand.first,middle_stand.second),gap_radius);

    //number of trees may decrease due to hard core
    int no_trees_0 = no_trees;
    GenerateLocations(no_trees,0.0,0.0,vx,vy,tree_distance,gap,locations);

    if (verbose){
      cout << "Number of trees" << locations.size() <<endl 
	   << " Density/ha wanted: " << (double)no_trees_0/(vx*vy/10000.0)
           << " Density/ha created: " << (double)no_trees/(vx*vy/10000.0) <<endl;
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
  void GrowthLoop<TREE, TS,BUD>::createTrees()
{
  //In the case of many trees & generated (random) positions
  //store positions & trees in the positions to be able to
  //repeat the run with the same configuration
  if(many_trees && generate_locations) {
    ofstream of("runfile.dat", ofstream::trunc);
    of << " x  y  treefile" << endl;
    of.close();
  }

  for (int i = 0; i < no_trees; i++){
    pair<double,double> p = locations[i];

    TREE* t = new TREE(Point(p.first,p.second,0.0),PositionVector(0,0,1),
		       "sf.fun","fapical.fun","fgo.fun",
		       "fsapwdown.fun","faf.fun","fna.fun", "fwd.fun",
		       "flr.fun");

    if(no_compartments) {
      XMLDomTreeReader<ScotsPineSegment,ScotsPineBud> reader;
      reader.readXMLToTree(*t, comp_tree);

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
      XMLDomTreeReader<ScotsPineSegment,ScotsPineBud> reader;
      if(!many_trees) {
	reader.readXMLToTree(*t, input_tree_file);
      }
      else {
	LGMdouble no_many_trees = (LGMdouble)tree_files.size();
	unsigned int tree = (unsigned int)((int)(ran3(&ran3_seed) * no_many_trees));
	reader.readXMLToTree(*t, tree_files[tree]);
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
	cout << h << endl;
	exit(0);
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


	cout << " Tree H (m) Dbh (cm), Dbase (cm) Hcrown_base (m) Dcrown_base (cm)  Wf (kg dm)  Af (m2) Wf_Repola"
	     " Wf_Repolab"  << endl;
	cout <<  input_tree_file << " "
	     << h << " " << 100.0*dbh << " " << 100.0*d_base << " " << h_cb << " "
	     << 100.0*d_cb << " " << 2.0*Wf << " " << treeAf << " " << Wf_Repola
	     << " " << Wf_Repolab << endl;

	exit(0);
      }

      LGMdouble Af = 0.0;
      Af = Accumulate(*t,Af,CollectFoliageArea<TS,BUD>());
      cout << "Af " << Af << endl;

      MoveTree<TS,BUD> move(Point(p.first,p.second,0.0)-GetPoint(*t),*t);
      ForEach(*t, move);

    } //  if(!voxel_tree) ...


    if (verbose){
      cout << "Created a tree at: " << p.first << " " << p.second <<endl;
    }
    vtree.push_back(t);
  } //for(int i = ...
}    //::createTrees(

template<class TREE, class TS,class BUD>
  void GrowthLoop<TREE, TS,BUD>::initializeVoxelSpace()
{
  ifstream vf(voxelfile.c_str());
  LGMdouble vx,vy,vz; vx = vy = vz = 0.0;
  LGMdouble s1,s2,s3,b1,b2; s1 = s2 = s3 = 0.0; b1=b2=0.0;
  vf >> vx >> vy >> vz >> s1 >> s2 >> s3 >> b1 >> b2;
  vf.close();
  //vx,vy,vz define the voxel space  dimensions, s1, s2, s3 define the
  //voxel box dimensions
  if (verbose){
    cout << "Voxel Space: " << vx << " " << vy << " " << vz << " " 
	 << s1 << " " << s2 << " " << s3 << " " << b1 << " " << b2 <<endl;
 
    //NOTE: this overrides VoxelSpace.txt
    if(voxboxside > 0) {
      s1 = s2 = s3  = voxboxside;
    }
  }
  vs = new VoxelSpace(Point(0,0,0),Point(vx,vy,vz),
		      s1,s2,s3,
		      static_cast<int>(vx/s1),static_cast<int>(vy/s2),static_cast<int>(vz/s3),
		      GetFirmament(*vtree[0]));

   //Now the structure of the tree in voxel space is specified

  if(voxel_tree) {

    cout << "Radiation calculations now with voxel tree in voxel space accordig to VoxelSpace.txt" << endl;

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
	  if(sqrt(pow(p.getX()-ms.getX(),2.0)+pow(p.getY()-ms.getY(),2.0)) > target_tree_rad) {
	    Point loc = p - Point(0.0,0.0,deltaZ/2.0);
	    
	    Point p1 = p - ms;
	    
	    //	    double ker = (0.5/sqrt(2.0))/sqrt(pow(p1.getX(),2.0)+pow(p1.getY(),2.0));
	    //	    PositionVector suunta(ker*p1.getX(),ker*p1.getY(),sqrt(1.0-0.25));

	    ts = new TS(loc, up, 1.0, L, Rw, Rh, t);
	    //	    ts = new TS(loc, suunta, 1.0, L, Rw, Rh, t);

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
	}

    //Write the voxeltree to a file
    XMLDomTreeWriter<TS,BUD> writer;
    writer.writeTreeToXML(*t,"voxel-tree.xml");

    //and information about voxelspace 
    //    DumpCfTree(*vs, *t, num_parts, wood_voxel);
    DumpCfTree(*vs, *t, num_parts, true);
    //After dumping all trees to VoxelSpace it is necessary to evaluate
    //sum quantities (e.g. STAR_mean)
    vs->updateBoxValues();


    cout << "Xn Yn Zn " << nx << " " << ny << " " << nz << endl;
    vs->writeVoxelSpaceContents();

    LGMdouble fa = 0.0;
    fa = Accumulate(*t,fa,CollectFoliageArea<TS,BUD>());
    cout << "Foliage area of the voxel tree " << fa << "  m2" << endl;

  }  //if(voxel_tree) ....

  //End of virittely valoa varten
}  //End of    ::createTrees() ...

template<class TREE, class TS,class BUD>
void GrowthLoop<TREE, TS,BUD>::initializeFunctions()
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
  void GrowthLoop<TREE, TS,BUD>::setVoxelSpaceAndBorderForest()
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
  //  FindCfBoundingBox<TS,BUD> fb(!wood_voxel);
  FindCfBoundingBox<TS,BUD> fb(false);

  for (unsigned int k = 0; k < (unsigned int)no_trees; k++) {
    bb = Accumulate(*vtree[k], bb, fb);
  }

  Point ll = bb.getMin();
  Point ur = bb.getMax();

  if(!voxel_tree) {                        //tehty jo voxel treelle
    vs->resize(ll, ur);

    //Note: not first tree, since it will be the first to be calculated
    //and its foliage won't be in the voxelspace 
    // unless dump_self is set

    vs->reset();

    for (unsigned int k = 0; k < (unsigned int)no_trees; k++) {
      //      DumpCfTree(*vs, *vtree[k], num_parts, wood_voxel);
      DumpCfTree(*vs, *vtree[k], num_parts, true);
    }

    //After dumping all trees to VoxelSpace it is necessary to evaluate
    //sum quantities (e.g. STAR_mean)
    vs->updateBoxValues();
 
    cout << endl << "getArea() " << vs->getArea() 
	 << " getLowerLeftCorner() " << vs->getLowerLeftCorner();
    cout << " getUpperRightCorner() " << vs->getUpperRightCorner();
    cout << " getNumberOfBoxes() " << vs->getNumberOfBoxes() << " getNumberOfFilledBoxes() "
	 << vs->getNumberOfFilledBoxes() << " getNumberOfTreeSegments() "
	 << vs->getNumberOfTreeSegments() << endl;
    cout << " getBoxVolume() " << vs->getBoxVolume() << " getXSideLength() " << vs->getXSideLength() << endl;
    cout << " getYSideLength() " << vs->getYSideLength() << " getZSideLength() " << vs->getZSideLength() << endl;
    cout << " getNoBoxX() " << vs->getNoBoxX() << " getNoBoxY() " << vs->getNoBoxY()
	 << " getNoBoxZ() " << vs->getNoBoxX() << endl;
    cout << " getNeedleArea() " << vs->getNeedleArea() << " getLeafArea() " <<  vs->getLeafArea() << endl;
    cout << endl;
  }


  //Lasketaan LAI VoxelSpacen avulla
  //Tehdaan oma VoxelSpace sita varten
  //Target puu (alla: calculateRadiation()) tehdaan niin, etta segmentit ovat
  //tasoilla Hcb, ..., H, yhteensa 10 kpl. Kun tehdaan voxelboxeja 9 kpl
  //valilla [Hcb,H] niin vastaa target puun jakoa: 1. voxelboxtaso  varjostaa 2. anturia,
  //9. taso anturia 10.
  //Dumpataan kaikki puut voxelspaceen: target puu asetetaan sinne vasta myohemmin
  //(calculateRadiation()).
  //Ei ole niin valia, kuinka monta boxia x ja y suunnassa on (voisi olla 1).
  
  cout << "ll " << ll;
  cout << "ur " << ur;
  VoxelSpace lai_vs(ll, ur, 10, 10, 9, GetFirmament(*vtree[0]));
  for (unsigned int k = 0; k < (unsigned int)no_trees; k++)
    //    DumpCfTree(lai_vs, *vtree[k], num_parts, wood_voxel);
    DumpCfTree(lai_vs, *vtree[k], num_parts, true);

  vector<pair<LGMdouble,LGMdouble> > lai_h;
  LGMdouble Hmin, Hmax;
  int n_levels;
  lai_vs.evaluateVerticalNeedleAreaDensity(Hmax, Hmin, n_levels, lai_h);

  cout << "Needle area in space: " << lai_vs.getFoliageArea() << endl << endl;
    
  cout << "Hmin " << Hmin << "  Hmax " << Hmax << endl;
  cout << "left corner (" << ll.getX() << " , " << ll.getY() << " )     right corner  ( " 
       << ur.getX() << " , " << ur.getY() << " )" << endl;   
  cout << "Box_center    lai cum_lai" << endl;
  LGMdouble cum_lai = 0.0;
  for(int i = 0; i < n_levels; i++) {
    cum_lai += lai_h[i].second;
    cout << lai_h[i].first << " " << lai_h[i].second << " " << cum_lai << endl;
  }
  cout << "This for LAI voxel space" << endl << endl;
    cout << endl << "getArea() " << lai_vs.getArea() 
	 << " getLowerLeftCorner() " << lai_vs.getLowerLeftCorner();
    cout << " getUpperRightCorner() " << lai_vs.getUpperRightCorner();
    cout << " getNumberOfBoxes() " << lai_vs.getNumberOfBoxes() << " getNumberOfFilledBoxes() "
	 << lai_vs.getNumberOfFilledBoxes() << " getNumberOfTreeSegments() "
	 << lai_vs.getNumberOfTreeSegments() << endl;
    cout << " getBoxVolume() " << lai_vs.getBoxVolume() << " getXSideLength() "
	 << lai_vs.getXSideLength() << endl;
    cout << " getYSideLength() " << lai_vs.getYSideLength() << " getZSideLength() "
	 << lai_vs.getZSideLength() << endl;
    cout << " getNoBoxX() " << lai_vs.getNoBoxX() << " getNoBoxY() " << lai_vs.getNoBoxY()
	 << " getNoBoxZ() " << lai_vs.getNoBoxX() << endl;
    cout << " getNeedleArea() " << lai_vs.getNeedleArea() << " getLeafArea() "
	 <<  lai_vs.getLeafArea() << endl;
    cout << endl;


  if(evaluate_LAI) {
    exit(0);
  }

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
// calculateRadiation UnDumps a tree (except the first tree, since
// it was not dumped in setVoxelspace()) and calculates radiation
// 1) by pairwise for itself and 2) through voxelspace and
// 3) borderforest outside itself
// and then dumps it back to voxelspace
//===================================================================
template<class TREE, class TS,class BUD>
  void GrowthLoop<TREE, TS,BUD>::calculateRadiation()
{

  //Find max height (Hmax) and minimum height of crown base (Hcbmin)
  //and construct then the target tree
     // cout<<"this is inside calculateRadiation";

  BoundingBox bb;
  //  FindCfBoundingBox<TS,BUD> fb(true);
  FindCfBoundingBox<TS,BUD> fb(false);

  for (unsigned int k = 0; k < (unsigned int)no_trees; k++) {
    bb = Accumulate(*vtree[k], bb, fb);
  }

  cout << "ur " << bb.getMax();
  cout << "ll " << bb.getMin();

  LGMdouble Hmax = bb.getMax().getZ();
  LGMdouble Hcb = bb.getMin().getZ();

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

  GetFirmament(*t_t).resize(15,15,1200.0);
  SetValue(*t_t,TreeQinMax,GetFirmament(*t_t).diffuseBallSensor());

 for (unsigned int k = 0; k < (unsigned int)no_trees; k++) {
   GetFirmament(*vtree[k]).resize(15,15,1200.0);
 }

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

  //Finally set the target tree as the fisrt tree in the vtree vector 
  vtree.insert(vtree.begin(),t_t);
  no_trees++;

  cout << "H Hcb " << Hmax << " " << Hcb << endl;
  cout << endl << "Target tree at " << mp;

  //Tama taytyy jatkossa paremmin mutta simulointien perusteella nayttaa, etta boxissa
  //k/STAR_mean = 0.97 + segment_length/Box_edge. Kaytetaan segment_length = 0.1 m. Josta
  //k = (0.97 + (0.1/3-juuri(box_vol)))*0.14
  //Taytyy parantaa jatkossa
  //  LGMdouble green_ext = (0.97 + 2.0* 0.1/pow(vs->getBoxVolume(), 0.33333))*0.14;

  //Virittelya valonlaskentaa varten
  //myos valonlaskennan parametreja, joita luetaan tiedostosta

  LGMdouble a, b;
  ifstream ab_file("box-radiation.txt");
   char dummy[256];
   ab_file.getline(dummy,256);    //header
   ab_file >> a >> b;
   ab_file.close();

   cout << endl;
   cout << " a  b " << a << " " << b << endl;

   //  BoundingBox bb;
   //  FindCfBoundingBox<TS,BUD> fb1(!wood_voxel);
  FindCfBoundingBox<TS,BUD> fb1(false);

  for (unsigned int k = 0; k < (unsigned int)no_trees; k++) {
    bb = Accumulate(*vtree[k], bb, fb1);
  }

  Point ll = bb.getMin();
  Point ur = bb.getMax();

  vs->resize(ll, ur);

    //Note: not first tree, since it will be the first to be calculated
    //and its foliage won't be in the voxelspace 
    // unless dump_self is set

    vs->reset();

     for (unsigned int k = 1; k < (unsigned int)no_trees; k++)
       //      DumpCfTree(*vs, *vtree[k], num_parts, wood_voxel);
      DumpCfTree(*vs, *vtree[k], num_parts, true);


    cout << endl << "getArea() " << vs->getArea() 
	 << " getLowerLeftCorner() " << vs->getLowerLeftCorner();
    cout << " getUpperRightCorner() " << vs->getUpperRightCorner();
    cout << " getNumberOfBoxes() " << vs->getNumberOfBoxes() << " getNumberOfFilledBoxes() "
	 << vs->getNumberOfFilledBoxes() << " getNumberOfTreeSegments() "
	 << vs->getNumberOfTreeSegments() << endl;
    cout << " getBoxVolume() " << vs->getBoxVolume() << " getXSideLength() " << vs->getXSideLength() << endl;
    cout << " getYSideLength() " << vs->getYSideLength() << " getZSideLength() " << vs->getZSideLength() << endl;
    cout << " getNoBoxX() " << vs->getNoBoxX() << " getNoBoxY() " << vs->getNoBoxY()
	 << " getNoBoxZ() " << vs->getNoBoxX() << endl;
    cout << " getNeedleArea() " << vs->getNeedleArea() << " getLeafArea() " <<  vs->getLeafArea() << endl;
    cout << endl;



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

     //Virittelya valonlaskentaa varten
     //Tassa tulostetaan vtree vektorin 1. puun tiedot. Se on target puu

  if(one_time) {
     ForEach(*t, SegmentProductionBalance(phprodfile));
  }
  else {
    //int bs;                 // box STAR  0 = constant value,
                          //           1 = value of box
                          //           2 = value of box & mean direction effect
    //int corr;               // If STAR_eq -> STAR correction considered (0 = no, 1 = yes)
    //LGMdouble vox;          // Size of voxel (=length of side, assuming voxels are cubes)

    SegmentProductionBalanceInfo info;
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

    ForEach(*t, SegmentProductionBalance(phprodfile, info));
  }
 
     exit(0);

     //end of virittelya

   }
}  //end of calculateRadiation()  { ..


/* #define HIT_THE_FOLIAGE 1 */
/* #define NO_HIT 0 */
/* #define HIT_THE_WOOD -1 */

// 1)  EvaluateRadiationForCfTreeSegment_1s evaluates shading
// caused by all other segments on this conifer segment.  The shading caused by segments in the crown
// of tree itself is evaluated by ShadingEffectOfCfTreeSegment_1<TS,BUD> (it is the same as
// ShadingEffectOfCfTreeSegment in stl-lignum; the functor is duplicated in include/CalculateRadiation.h
// only for convenience).
// After that other trees are accounted for with voxelspace (voxel_space->getRoute() etc) and surrounding
//stand with border_forest->getBorderForestExtinction().
// This function evaluates the shading by surrounding trees and border forest separately (Qin_stand);
// EvaluateRadiationForCfTreeSegment_1 in include/CalculateLight.h does otherwise the same as this but does not
// evaluate Qin_stand.

// 2) EvaluateRadiationForCfTreeSegment_2 evaluates shading by all other segments by paiwise comparison (segments
// in own crown & other trees). It uses ShadingEffectOfCfTreeSegment_1<TS,BUD> to evaluate shading.
// This functor evaluates shading by own crown and shading by other trees (stand) separately and updates
// Qin_stand in TreeSegment.


/* //======================================================================================================= */

/* //This version of radiation evaluates radiation conditions for subject tree by pairwise */
/* // comparison */

/* //This functor EvaluateRadiationForCfTreeSegment evaluates shading */
/* //caused by all other segments on this conifer segment. This functor */
/* //uses functor ShadingEffectOfCfTreeSegment<TS,BUD> to go through all */
/* //segments to check the shading. */

/* //If the attributes voxel_space and border_forest are set (then */
/* //voxel_space != NULL), the attenuation of the bean in the voxel_space */
/* //and border_forest is taken into consideration. */

/* template <class TS, class BUD, class TREE> */
/* class EvaluateRadiationForCfTreeSegment_2 { */
/* public: */
/*     EvaluateRadiationForCfTreeSegment_2(const ParametricCurve& k, vector<TREE*>& vt): */
/*   K(k), vtree(vt) {only_self = false;} */
/*  EvaluateRadiationForCfTreeSegment_2(const ParametricCurve& k, vector<TREE*>& vt, bool o_s): */
/*       K(k), vtree(vt), only_self(o_s) {} */

/*       TreeCompartment<TS,BUD>* operator()(TreeCompartment<TS,BUD>* tc)const; */
/* private: */
/*   const ParametricCurve& K; */
/*   vector<TREE*>& vtree; */
/*   bool only_self;      */
/* }; */


/* //This functor EvaluateRadiationForCfTreeSegment evaluates shading */
/* //caused by all other segments on this conifer segment. This functor */
/* //uses functor ShadingEffectOfCfTreeSegment<TS,BUD> to go through all */
/* //segments to check the shading. */


/* template <class TS, class BUD, class TREE> */
/*   TreeCompartment<TS,BUD>* EvaluateRadiationForCfTreeSegment_2<TS,BUD,TREE>::operator() (TreeCompartment<TS, BUD>* tc)const */
/* { */
/*   if (TS* ts = dynamic_cast<TS*>(tc)){ */
/*     SetValue(*ts, LGAQin, 0.0); */
/*     SetValue(*ts, LGAQabs, 0.0); */

/*     //Radiation  conditions are not  evaluated if  the segment  has no */
/*     //foliage (in practice  there would be division by  0 in computing */
/*     //absorbed radiation) */
/*     if (GetValue(*ts, LGAWf) < R_EPSILON){ */
/* 	return tc; */
/*     } */

/*     Tree<TS,BUD>& tt = GetTree(*ts); */
/*     FirmamentWithMask& firmament = GetFirmament(tt); */
/*     int number_of_sectors = firmament.numberOfRegions(); */
/*     double a_dot_b = 0.0; */
/*     vector<double> radiation_direction(3); */

/*     vector<double> v(number_of_sectors,0.0);  */
/*     ShadingEffectOfCfTreeSegment_1<TS,BUD> s_e(ts,K,v); */


/*     //This  goes  through  the  tree  and computes  shading  based  on */
/*     //1)distance  light beam traverses  in foliage,  2)foliage density */
/*     //and 3) inclination light beam hits the segment. */

/*     //just go through all trees: first the others (stand) and then the tree itself. */
/*     //This is to evaluate the self-shading and stand-shading components. */

/*     //The target tree is the first tree in the vector; do it last in order */
/*     //to get the effect of surrounding stand */

/*     //In the case of only_self it is the only tree in the vector */
/*     if(only_self) { */
/*       TREE* t = vtree[0]; */
/*       ForEach(*t,s_e); */
/*     }  */
/*     else { */
/*       if(vtree.size() > 1) */
/* 	for(unsigned int k = 1; k < vtree.size(); k++) { */
/* 	  TREE* t = vtree[k]; */
/* 	  ForEach(*t,s_e); */
/* 	} */
/*     } */

/*     //Now the Qin after shading of others */
/*     vector<double> qis(number_of_sectors,0.0);  */
/*     vector<double>& ss = s_e.getS(); */
/*     for (int i = 0; i < number_of_sectors; i++){ */
/*       if (ss[i] != HIT_THE_WOOD){ */
/* 	MJ Io = firmament.diffuseRegionRadiationSum(i,radiation_direction); */
/* 	qis[i] = Io*exp(-ss[i]); */
/*       } */
/*     } */
/*     LGMdouble Qin_stand = 0.0; */
/*     Qin_stand = accumulate(qis.begin(),qis.end(),0.0); */
/*     ts->setQinStand(Qin_stand); */

/*     //Now the last tree i.e tree itself */
/*     //    ForEach(*(vtree[0]), s_e); */

/*     //implement  "Ip  =  Iope^(-Vp)",  s[i] =  radiation  coming  from */
/*     //direction i after this */
/*     vector<double>& s = s_e.getS(); */
/*     for (int i = 0; i < number_of_sectors; i++){ */
/*       if (s[i] == HIT_THE_WOOD){ */
/* 	s[i] = 0.0; */
/*       } */
/*       else { */
/* 	MJ Iop = firmament.diffuseRegionRadiationSum(i,radiation_direction); */
/* 	s[i] = Iop*exp(-s[i]); */
/*       } */
/*     } */
/*     //Total incoming radiation   */
/*     MJ Q_in = accumulate(s.begin(),s.end(),0.0); */

/*     //s contains now incoming radiation from each sector. Evaluate how */
/*     //much segment absorbs from incoming radation. */
/*     LGMdouble Lk, inclination, Rfk, Ack, extinction, sfk, Ask, Wfk; */
/*     Lk = Rfk = Ack =  extinction = sfk = Ask = Wfk = 0.0; */
/*     Lk = GetValue(*ts, LGAL);   //length is > 0.0, otherwise we would not bee here */
/*     Rfk = GetValue(*ts, LGARf);  //Radius to foliage limit  */
/*     Wfk = GetValue(*ts, LGAWf); //Foliage mass */
/*     //   sfk  = GetValue(tt, LGPsf); //Foliage m2/kg from tree */
/*     sfk  = GetValue(*ts, LGAsf); //Foliage m2/kg from segment!!! */

/*     for (int i = 0; i < number_of_sectors; i++){ */
/*       firmament.diffuseRegionRadiationSum(i,radiation_direction); */
/*       a_dot_b = Dot(GetDirection(*ts), PositionVector(radiation_direction)); */
/*       inclination = PI_DIV_2 - acos(fabs(a_dot_b)); */

/*       Ack = 2.0*Lk*Rfk*cos(inclination) + PI_VALUE*pow(Rfk,2.0)*sin(inclination); */

/*     extinction = (double)K(inclination); */

/*       if (Ack == 0.0){ */
/* 	cout << "ERROR EvaluateRadiationForCfTreeSegment: Ack == 0 (division by 0)" */
/* 	     << endl; */
/*       } */

/*       //implement I(k)p = Ip*Ask, Note  Ack must be greater than 0 (it */
/*       //should if there is any foliage) */
/*       Ask = (1.0 - exp(-extinction*((sfk*Wfk)/Ack)))*Ack; */
/*       s[i] *= Ask; */
/*     } */

/*     MJ Q_abs = accumulate(s.begin(),s.end(),0.0); */
/*     SetValue(*ts, LGAQabs, Q_abs); */
/*     SetValue(*ts, LGAQin, Q_in); */
/*   } */
/*   return tc; */
/* } */



/* //======================================================== */

/* class AccumulateOpticalDepth{ */
/*  public: */
/*  AccumulateOpticalDepth(LGMdouble side, LGMdouble a, LGMdouble b, Point loc, ParametricCurve kk, */
/* 			bool d_e, bool wd, bool cs, LGMdouble st) : */
/*   box_side_length(side), par_a(a), par_b(b), seg_loc(loc), K(kk), dir_effect(d_e), wood(wd), */
/*     constant_star(st),correct_star(cs) {box_volume = pow(box_side_length,3.0);} */
/*   double operator()(double o_d,VoxelMovement& vm){ */
/*     //    if((vm.af > R_EPSILON ||(wood && vm.wood_area > R_EPSILON)) && vm.n_segs_real > 0.0) { */
/*     if(vm.af > R_EPSILON ||(wood && vm.wood_area > R_EPSILON)) { */
/*       LGMdouble k; */
/*       if(constant_star > 0.0) */
/* 	k = constant_star; */
/*       else */
/* 	k = vm.STAR_mean; */

/*       if(correct_star) { */
/* 	k = max(0.0,-0.014+1.056*k); */
/*       } */

/*       //NOTE: here transformation STAR_eq --> STAR; documented in */
/*       //~/Riston-D/E/LIGNUM/Light/summer-09-test/STAR-vs-STAR_eq.pdf */

/*       //Effect of hit angle to the mean direction of shoots in the voxel box, documented in */
/*       //~/Riston-D/E/LIGNUM/Light/Article/vs-STARmean-all.pdf and */
/*       //~/Riston-D/E/LIGNUM/Light/Article/vs-STARmean-approximation.pdf */
/*       PositionVector mean_dir = vm.mean_direction; */
/*       LGMdouble mean_dir_length = mean_dir.length(); */
/*       LGMdouble effect = 1.0; */
/*       if(dir_effect) { */
/* 	if(mean_dir_length > 0.0){ */
/* 	  mean_dir.normalize(); */
/* 	  LGMdouble inclination  =  PI_DIV_2 - acos(fabs(Dot(mean_dir,beam_dir))); */
/* 	  effect =  K(inclination)/K(0.7); */
/* 	} */
/*       } */
/*       //this scales the effect depending on how parallel the segments are */
/*       /\*       if(mean_dir_length > 0.0) { *\/ */
/*       /\* 	mean_dir.normalize(); *\/ */
/*       /\* 	LGMdouble inclination  =  PI_DIV_2 - acos(fabs(Dot(mean_dir,beam_dir))); *\/ */
/*       /\* 	//	LGMdouble effect  = 1.13 - 0.24*pow(inclination,2.0); *\/ */
/*       /\* 	//	LGMdouble effect = 1.2-0.3*inclination*(inclination+0.3); *\/ */
/*       /\* 	LGMdouble u =  mean_dir_length/vm.n_segs_real; *\/ */
/*       /\* 	effect = 1.0 - u + u * K(inclination)/K(0.7); *\/ */
/*       /\* 	cout << " u " << u << endl; *\/ */
/*       /\*       } *\/ */

/*       o_d += effect * k * vm.af * vm.l / box_volume; */

/*       if(wood) { */
/* 	//Mean projection area of surface of a circular cylinder (excluding end disks) */
/* 	// is 1/4 of its area */
/* 	o_d += 0.25 * vm.wood_area * vm.l / box_volume; */
/*       } */

/*     } */
/*     return o_d; */
/*   } */

/*     //This is to take care of the direction of the beam of radiation */
/*     PositionVector beam_dir; */
/*  private: */
/*     LGMdouble box_side_length; */
/*     LGMdouble box_volume; */
/*     LGMdouble par_a, par_b; */
/*     Point seg_loc;   //location of segment */
/*     ParametricCurve K; */
/*     bool dir_effect;         //If direction effect of segments in box considered */
/*     bool wood;               //If woody parts are considered */
/*     LGMdouble constant_star; //If this > 0, k = constant_star else k = star_mean */
/*     bool correct_star;       //If star_eq -> star correction is done */
/* }; */

//Tassa valonlasenta voxelspacessa

// Calculation of radiation: voxel space & border forest with storing of
// Qin_stand into the segment (is otherwise as EvaluateRadiationForCfTreeSegment_1)

//This functor EvaluateRadiationForCfTreeSegment evaluates shading
//caused by all other segments on this conifer segment. This functor
//uses functor ShadingEffectOfCfTreeSegment<TS,BUD> to go through all
//segments to check the shading.

//If the attributes voxel_space and border_forest are set (then
//voxel_space != NULL), the attenuation of the bean in the voxel_space
//and border_forest is taken into consideration.



//This functor EvaluateRadiationForCfTreeSegment evaluates shading
//caused by all other segments on this conifer segment. This functor
//uses functor ShadingEffectOfCfTreeSegment<TS,BUD> to go through all
//segments to check the shading.

/* template <class TS, class BUD> */
/*   TreeCompartment<TS,BUD>* EvaluateRadiationForCfTreeSegment_3<TS,BUD>::operator() (TreeCompartment<TS, BUD>* tc)const */
/* { */
/*   if (TS* ts = dynamic_cast<TS*>(tc)){ */
/*     SetValue(*ts, LGAQin, 0.0); */
/*     SetValue(*ts, LGAQabs, 0.0); */
/*     //Radiation  conditions are not  evaluated if  the segment  has no */
/*     //foliage (in practice  there would be division by  0 in computing */
/*     //absorbed radiation) */
/*     if (GetValue(*ts, LGAWf) < R_EPSILON){ */
/* 	return tc; */
/*     } */

/*     Tree<TS,BUD>& tt = GetTree(*ts); */
/*     FirmamentWithMask& firmament = GetFirmament(tt); */
/*     int number_of_sectors = firmament.numberOfRegions(); */
/*     double a_dot_b = 0.0; */
/*     vector<double> radiation_direction(3); */
/*     Point middle = GetMidPoint(*ts); */

/*     vector<double> v(number_of_sectors,0.0);  */
/*     //    ShadingEffectOfCfTreeSegment_1<TS,BUD> s_e(ts,K,v); */
/*     //This  goes  through  the  tree  and computes  shading  based  on */
/*     //1)distance  light beam traverses  in foliage,  2)foliage density */
/*     //and 3) inclination light beam hits the segment. */
/*     //HUOMMMM   !!!!!!!! */
/*     //ForEach(tt,s_e); */
    
/*     //implement  "Ip  =  Iope^(-Vp)",  s[i] =  radiation  coming  from */
/*     //direction i after this */
/*     //    vector<double>& s = s_e.getS(); */
/*     vector<double> s(number_of_sectors, 0.0); */
/*     vector<double> qis(number_of_sectors, 0.0); */

/*     AccumulateOpticalDepth AOD(voxel_space->getXSideLength(), par_a, par_b, middle,K, */
/* 			       dir_effect, wood, correct_star, constant_star); */
/*     for (int i = 0; i < number_of_sectors; i++){ */
/*       MJ Iop = firmament.diffuseRegionRadiationSum(i,radiation_direction); */
      
/* 	//first attenuation in the voxel space */
/* 	LGMdouble transmission_voxel = 1.0; */
/* 	vector<VoxelMovement> vm; */
/* 	PositionVector dir(radiation_direction); */

/* 	voxel_space->getRoute(vm, middle, dir, K, false);  //this shoud return only the "box route" */

/* 	//with traveled lengths */
/* 	//calculate the extinction coeffient */
/* 	//Consider also the mean direction of shoots in box */
/* 	AOD.beam_dir = dir; */
    
/* 	LGMdouble optical_depth = accumulate(vm.begin(),vm.end(),0.0,AOD); */

/* 	//If tree itself is in calculation (dump_self = true) subtract the effect of own foliage */
/* 	//NOTE Assumes that all own foliage is in the first (= the one that contains the middle point */
/* 	//of the segment) voxel box */

/* 	if(dump_self) { */
/* 	  LGMdouble k =  max(0.0,-0.014+1.056*vm[0].STAR_mean); */
/* 	  optical_depth -= k * GetValue(*ts,LGAAf) * vm[0].l / voxel_space->getBoxVolume(); */

/* 	  if(optical_depth < 0.0) */
/* 	    optical_depth = 0.0; */
/* 	} */
/* 	if(optical_depth > R_EPSILON) */
/* 	  if(optical_depth < 20.0) */
/* 	    transmission_voxel = exp(-optical_depth); */
/* 	  else */
/* 	    transmission_voxel = 0.0; */
/* 	Iop *= transmission_voxel; */

/* 	//then attenuation in the BorderForest */
/* 	if(evaluate_border) */
/* 	  Iop *= border_forest->getBorderForestExtinction(middle, dir,k_border_conifer); */
 
/*       qis[i] = Iop; */
/*       s[i] = Iop; */


/* /\*       if (s[i] == HIT_THE_WOOD){ *\/ */
/* /\* 	s[i] = 0.0; *\/ */
/* /\*       } *\/ */
/* /\*       else *\/ */
/* /\* 	s[i] = Iop*exp(-s[i]); *\/ */
      
/*     } //End of no_sectors ... */


/*    //Total incoming radiation and radiation after stand */
/*     LGMdouble Qin_stand = accumulate(qis.begin(),qis.end(),0.0); */
/*     ts->setQinStand(Qin_stand); */

/*     MJ Q_in = accumulate(s.begin(),s.end(),0.0); */
    
/*     //s contains now incoming radiation from each sector. Evaluate how */
/*     //much segment absorbs from incoming radation. */
/*     LGMdouble Lk, inclination, Rfk, Ack, extinction, sfk, Ask, Wfk; */
/*     Lk = Rfk = Ack =  extinction = sfk = Ask = Wfk = 0.0; */
/*     Lk = GetValue(*ts, LGAL);   //length is > 0.0, otherwise we would not bee here */
/*     Rfk = GetValue(*ts, LGARf);  //Radius to foliage limit  */
/*     Wfk = GetValue(*ts, LGAWf); //Foliage mass */
/*     //sfk  = GetValue(tt, LGPsf); //Foliage m2/kg from tree */
/*     sfk  = GetValue(*ts, LGAsf); //Foliage m2/kg from segment!!! */

/*     for (int i = 0; i < number_of_sectors; i++){ */
/*       firmament.diffuseRegionRadiationSum(i,radiation_direction); */
/*       a_dot_b = Dot(GetDirection(*ts), PositionVector(radiation_direction)); */
/*       inclination = PI_DIV_2 - acos(fabs(a_dot_b)); */

/*       Ack = 2.0*Lk*Rfk*cos(inclination) + PI_VALUE*pow(Rfk,2.0)*sin(inclination); */
/*       extinction = (double)K(inclination); */

/*       if (Ack == 0.0){ */
/* 	cout << "ERROR EvaluateRadiationForCfTreeSegment: Ack == 0 (division by 0)" */
/* 	     << endl; */
/*       } */

/*       //implement I(k)p = Ip*Ask, Note  Ack must be greater than 0 (it */
/*       //should if there is any foliage) */
/*       Ask = (1.0 - exp(-extinction*((sfk*Wfk)/Ack)))*Ack; */
/*       s[i] *= Ask; */
/*     } */
/*     MJ Q_abs = accumulate(s.begin(),s.end(),0.0); */
/*     SetValue(*ts, LGAQabs, Q_abs); */
/*     SetValue(*ts, LGAQin, Q_in); */
/*   } */
/*   return tc; */
/* }  //end of EvaluateRadiationForCfTreeSegment()  { ... */



//===================================================================
// calculateRadiationOnlySelf()
//===================================================================

template<class TREE, class TS,class BUD>
  void GrowthLoop<TREE, TS,BUD>::calculateRadiationOnlySelf(){

  //1) Create the tree
  TREE* t = new TREE(Point(10.0,10.0,0.0),PositionVector(0,0,1),
		     "sf.fun","fapical.fun","fgo.fun",
		     "fsapwdown.fun","faf.fun","fna.fun", "fwd.fun",
		     "flr.fun");

  XMLDomTreeReader<ScotsPineSegment,ScotsPineBud> reader;
  reader.readXMLToTree(*t, input_tree_file);

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
    ifstream vf(voxelfile.c_str());
    LGMdouble vx,vy,vz; vx = vy = vz = 0.0;
    LGMdouble s1,s2,s3,b1,b2; s1 = s2 = s3 = 0.0; b1=b2=0.0;
    vf >> vx >> vy >> vz >> s1 >> s2 >> s3 >> b1 >> b2;
    vf.close();
    //vx,vy,vz define the voxel space  dimensions, s1, s2, s3 define the
    //voxel box dimensions
    if (verbose){
      cout << "Voxel Space: " << vx << " " << vy << " " << vz << " " 
	   << s1 << " " << s2 << " " << s3 << " " << b1 << " " << b2 <<endl;
 
      //NOTE: this overrides VoxelSpace.txt
      if(voxboxside > 0) {
	s1 = s2 = s3  = voxboxside;
      }
    }
    vs = new VoxelSpace(Point(0,0,0),Point(vx,vy,vz),
			s1,s2,s3,
			static_cast<int>(vx/s1),static_cast<int>(vy/s2),static_cast<int>(vz/s3),
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
    ForEach(*t, SegmentProductionBalance(phprodfile));
  }
  else {
    //int bs;                 // box STAR  0 = constant value,
    //           1 = value of box
    //           2 = value of box & mean direction effect
    //int corr;               // If STAR_eq -> STAR correction considered (0 = no, 1 = yes)
    //LGMdouble vox;          // Size of voxel (=length of side, assuming voxels are cubes)

    SegmentProductionBalanceInfo info;
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

    ForEach(*t, SegmentProductionBalance(phprodfile, info));
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
  void GrowthLoop<TREE,TS,BUD>::getTreesAndPositions() {  
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
