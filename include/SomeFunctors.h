#ifndef SOME_FUNCTORS_H
#define SOME_FUNCTORS_H
#include <sstream>
#include <Lignum.h>
#include <ScotsPine.h>
#include <VoxelSpace.h>


extern int ran3_seed;
typedef pair<pair<double,double>,double> ForestGap;

//Propagate up the Qin to newly created segments and buds
//Also set the relative light LGAip for the bud
class ForwardScotsPineQin{
public:
  double operator()(double& qin, TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const
  {
    if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
      double ts_qin = GetValue(*ts,LGAQin);
      if (GetValue(*ts,LGAage) == 0.0){//pass qin to this segment
	SetValue(*ts,LGAQin,qin);
      }
      else{ //computed
	qin = ts_qin;
      }
    }
    else if (ScotsPineBud* b = dynamic_cast<ScotsPineBud*>(tc)){
      SetValue(*b,LGAQin,qin);
      Tree<ScotsPineSegment,ScotsPineBud>& t = 
	dynamic_cast<Tree<ScotsPineSegment,ScotsPineBud>&>(GetTree(*b));
      double max_light = GetValue(t,TreeQinMax);
      double ip = qin/max_light;
      SetValue(*b,LGAip,ip);
    }
    return qin;
  }
};

//Specific leaf  area for Scots pine  is a function  of relative light
//and set once for the newly created segments after the calculation of
//the light climate and the use of ForwardScotsPineQin. 
class SetScotsPineSegmentSf{
public:
  void operator()(TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const
  {
    if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
      if (GetValue(*ts,LGAage) == 0.0)
	SetValue(*ts,LGAsf);
    }
  }
};

//Set the needle angle once when  the segment is created as a function
//of light
class SetScotsPineSegmentNeedleAngle{
public:
  void operator()(TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const
  {
    if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
      if (GetValue(*ts,LGAage) == 0.0){
	const ScotsPineTree& t = dynamic_cast<ScotsPineTree&>(GetTree(*ts));
	const ParametricCurve& fna = GetFunction(t,SPFNA);
	double nl = GetValue(t,LGPnl);
	double B = GetValue(t,TreeQinMax);
	double qin = GetValue(*ts,LGAQin);
	//Update foliage height and radius to foliage limit
	SetValue(*ts,LGAHf,nl*sin(fna(qin/B)));
	SetValue(*ts,LGARf,GetValue(*ts,LGAR)+GetValue(*ts,LGAHf));
      }
    }
  }
};

class summing {
 public:
  summing():d2(0.0),d2l(0.0),lsum(0.0),n_br(0){}
  double d2;
  double d2l;
  double lsum;
  int n_br;
};
  

class Branchmeans{
 public:
  summing& operator ()(summing&, TreeCompartment<ScotsPineSegment,
		       ScotsPineBud>* ts)const;
};

//Set  STAR  mean based  on  segment age.  See  Smolander  et al  1994
//TreePhys.    Simply   initialize   the  functor   with   appropriate
//ParametricCurve and use ForEach.
template <class TS, class BUD>
class SetStarMean{
public:
  SetStarMean(const ParametricCurve& f):fstar_mean(f){}
  TreeCompartment<ScotsPineSegment,ScotsPineBud>* operator() (TreeCompartment<TS,BUD>* tc) const{
    if (TS* ts = dynamic_cast<TS*>(tc)){
      LGMdouble starm = fstar_mean(GetValue(*ts,LGAage));
      SetValue(*ts,LGAstarm,starm);
    }
    return tc;
  }
private:
const ParametricCurve& fstar_mean;
};


//UnDumpScotsPineSegment  is   the  inverse  to  DumpScotsPineSegment,
//instead of adding foliage remove it

inline void UnDumpScotsPineSegment(VoxelBox &b, ScotsPineSegment& ts,double num_parts,
				   bool wood_voxel){	

  LGMdouble fmass = GetValue(ts, LGAWf) / (double)num_parts;  
  LGMdouble farea = GetValue(ts,LGAAf) / (double)num_parts;	
  b.subtractNeedleArea(farea);
  b.subtractNeedleMass(fmass);

  //Here  weight the  star mean  with  foliage area,  then in  the
  //VoxelBox::updateValues() remember to have: starm=starSum/AfTot
  b.subtractStarSum(GetValue(ts,LGAstarm)*farea);
  b.subtractWeight(farea);

  b.decreaseNumberOfSegments();

  b.addNumberOfSegmentsReal(-1.0/(double)num_parts); 

  if(wood_voxel) {
      LGMdouble r = GetValue(ts, LGAR);
      LGMdouble length = GetValue(ts, LGAL) / (double)num_parts;
      LGMdouble mass = GetValue(ts, LGAWood) / (double)num_parts;
      LGMdouble area = 2.0*PI_VALUE*r*length;
      b.subtractWoodMass(mass);
      b.subtractWoodArea(area);
   }
 }



template <class TS, class BUD>
  class  UnDumpScotsPineTreeFunctor{
 public:
  UnDumpScotsPineTreeFunctor(VoxelSpace &s, int parts, bool wood_vox):space(s),num_parts(parts),
    wood_voxel(wood_vox){}  
    TreeCompartment<TS,BUD>* 
      operator ()(TreeCompartment<TS,BUD>* tc)const{
      if (TS* cfts = dynamic_cast<TS*>(tc)){
	bool foliage = false;
	if (GetValue(*cfts,LGAWf) > R_EPSILON)
	  foliage = true;
	if(foliage || wood_voxel) {
	  Point p = GetPoint(*cfts);
	  PositionVector pv = GetDirection(*cfts);
	  LGMdouble length = GetValue(*cfts, LGAL);
	  //if the user wants 1 part (whole segment), the loop is executed
	  //once  and the midpoint  of the  segment is  used; if  the user
	  //wants 2 parts,  the loop is executed twice  and the points 1/3
	  //and 2/3 of the segment length are used, and so on
	  for (int i=1; i<(num_parts+1.0); i++){
	    Point p1 = p + (Point)(((double)i/((double)num_parts+1.0))*length*pv);
	    try{
	      VoxelBox& box = space.getVoxelBox(p1);
	      UnDumpScotsPineSegment(box,*cfts, num_parts,wood_voxel);
	    }
	    catch (OutOfVoxelSpaceException e){
	      ostringstream error;
	      Point b = e.getBox();
	      Point p = e.getPoint();
	      error << "DumpScotPineTree Functor Out of voxel space " << endl 
		    << "Box " << b.getX() << " " << b.getY() << " " << b.getZ() << " "
		    << "Point " << p.getX() << " " << p.getY() << " " << p.getZ() << endl;
	      LGMError(error.str());
	    }
  
	  }  //for(int i = 1 ....
	} //if(foliage || ...
      } 
      return tc;
    }
 private:
    VoxelSpace &space;
    double num_parts;//number of parts the segment is divided into (usually
    //1 but one may want to experiment)
    bool wood_voxel; //If woody part is considered also
};


//UnDumpScotsPineTree is the "inverse" to DumpCfTree in VoxelSpaceI.h and VoxelVoxI.h,
//it removes the values of this tree from the voxels
inline void UnDumpScotsPineTree(VoxelSpace &s, 
				Tree<ScotsPineSegment,ScotsPineBud> &tree,int num_parts,
				bool wood_voxel){
  UnDumpScotsPineTreeFunctor<ScotsPineSegment,ScotsPineBud> f(s,num_parts,wood_voxel);
  ForEach(tree, f);
}


inline
void SetScotsPineSegmentQabs(VoxelBox &b, ScotsPineSegment& ts, double num_parts)
{
  //Qin is needed for LGPsf
  SetValue(ts, LGAQin, GetValue(ts,LGAQin)+(b.getQin()/num_parts)); 
  //LGPAAf is LGPsf*LGAWf and LGPsf is now a function of light
  LGMdouble farea = GetValue(ts, LGAAf) / num_parts;
  LGMdouble qabs = 0.0;
  //Qabs computetd based on Qin, mean star and foliage area.
  qabs = b.getQin()*GetValue(ts,LGAstarm)*farea;
  if (qabs > b.getQin())
    cout << "QABS > QIN " <<  qabs << " " << b.getQin() << " "  << GetValue(ts,LGAstarm) << " "  << farea<< endl;
  //LGAAf, farea,  already takes account the num_parts
  SetValue(ts, LGAQabs,GetValue(ts,LGAQabs)+qabs);  
}

template <class TS, class BUD>
class SetScotsPineTreeQabsFunctor{
public:
  SetScotsPineTreeQabsFunctor(int parts):space(NULL),num_parts(parts){}
  TreeCompartment<TS,BUD>* operator()(TreeCompartment<TS,BUD>* tc)const;
  VoxelSpace* space;
private:
  double num_parts;
};

template <class TS, class BUD>
TreeCompartment<TS,BUD>* SetScotsPineTreeQabsFunctor<TS,BUD>::operator ()(TreeCompartment<TS,BUD>* tc)const
{
  if (TS* cfts = dynamic_cast<TS*>(tc))
    {
      LGMdouble fmass = GetValue(*cfts, LGAWf);
      //Reset Qabs Qin  here
      SetValue(*cfts, LGAQabs, 0.0);
      SetValue(*cfts, LGAQin, 0.0);
      if (fmass > R_EPSILON)
	{
	  Point p = GetPoint(*cfts);
	  PositionVector pv = GetDirection(*cfts);
	  LGMdouble length = GetValue(*cfts, LGAL);
	  //if the user wants 1 part (whole segment), the loop is executed
	  //once  and the midpoint  of the  segment is  used; if  the user
	  //wants 2 parts,  the loop is executed twice  and the points 1/3
	  //and 2/3 of the segment length are used, and so on
	  for (int i=1; i<(num_parts+1.0); i++)
	    {
	      Point p1 = p + (Point)(length * (i/(num_parts+1.0)) * pv);
	      VoxelBox box = space->getVoxelBox(p1);
	      SetScotsPineSegmentQabs(box, *cfts, num_parts);
	    }		
	}
    }
  return tc;    
}

template <class TS,class BUD>
void SetScotsPineTreeQabs(VoxelSpace &s, Tree<TS, BUD> &tree, int num_parts)
{
  SetScotsPineTreeQabsFunctor<TS,BUD> f(num_parts);
  f.space = &s;
  ForEach(tree, f);
}

    class SetSapwoodDemandAtJunction {
  public:
      double& operator()(double& oomega, TreeCompartment<ScotsPineSegment,
			 ScotsPineBud>* tc)const
      {
        if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment *>(tc))
          {
            double myomega = GetValue(*ts,LGAomega);
            if(myomega - R_EPSILON  > oomega) {
              const ParametricCurve& pc =
 		GetFunction(dynamic_cast<ScotsPineTree&>(GetTree(*ts)),SPSD);
              SetValue(*ts,SPAAsDown,pc(oomega));
              oomega = myomega;
            }
            else
              SetValue(*ts,SPAAsDown,1.0);
          }
        return oomega;
      }
  };

class BranchInformation {
public:
  BranchInformation(const string& filename): fname(filename)
  {
    ofstream f(fname.c_str() , ofstream::trunc);
    f << "Height Dist_top Diameter HwDiam Length Foliage_mass  Sapwood_mass Heartwood_mass Wood_mass"
         " Photosynthesis Respiration rel_pos" << endl;
    f.close();
  }

  TreeCompartment<ScotsPineSegment,ScotsPineBud>*
		operator()(TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc) const
  {
    if (Axis<ScotsPineSegment,ScotsPineBud>*ax =
	dynamic_cast<Axis<ScotsPineSegment,ScotsPineBud>*>(tc)){
     list<TreeCompartment<ScotsPineSegment,ScotsPineBud>*>& cmpls =
	GetTreeCompartmentList(*ax);
      if(cmpls.size() > 2) {
	TreeSegment<ScotsPineSegment,ScotsPineBud>* fs = GetFirstTreeSegment(*ax);
	if((GetValue(*fs, LGAomega) > 1.0) && (GetValue(*fs, LGAomega) < 3.0)) {
	  Tree<ScotsPineSegment,ScotsPineBud>& tt = GetTree(*fs);
	  LGMdouble h0 = GetValue(tt,LGAH);
	  LGMdouble hh = GetPoint(*fs).getZ();
	  LGMdouble d = 200.0*GetValue(*fs, LGAR);
	  LGMdouble dh = 200.0*GetValue(*fs, LGARh);
	  LGMdouble Wf = 2.0*GetBranchFoliage(*ax);
	  LGMdouble Ws = 2.0*GetBranchSapwoodMass(*ax);
	  LGMdouble Wh = 2.0*GetBranchHeartwoodMass(*ax);
	  LGMdouble Ww = 2.0*GetBranchWoodMass(*ax);
	  LGMdouble P = 2.0*GetBranchPhotosynthesis(*ax);
	  LGMdouble M = 2.0*GetBranchRespiration(*ax);

	  //Relative position of the branch (0 = top, 1 = crown base)
	  Tree<ScotsPineSegment,ScotsPineBud>& tree = GetTree(*fs);
	  DCLData dcl;
	  AccumulateDown(tree,dcl,AddBranchWf(),DiameterCrownBase<ScotsPineSegment,ScotsPineBud>());
	  LGMdouble Hc = dcl.HCrownBase();
	  LGMdouble rel_pos = (h0 - hh) /(h0 - Hc);
	  rel_pos = max<double>(0.0,min<double>(1.0,rel_pos));

	  ofstream f(fname.c_str() , ofstream::app);
	  f << hh << " " << h0 - hh << " " << d << " " << dh << " " << GetValue(*ax,LGAL) << " "
	    << Wf << " " << Ws << " " << Wh << " " << Ww << " " << P << " " << M << " " << rel_pos << endl;
	  f.close();
	}
      }

    }
    return tc;
  }

 private:
  string fname;
};



//Evaluate the surface area of the woody part (excluding end discs) of
//new segments (that have foliage)

class SurfaceAreaOfNewSegments{
public:
  double operator()(double& A, TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const
  {
    if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
      if (GetValue(*ts,LGAage) > R_EPSILON && GetValue(*ts,LGAWf) > R_EPSILON) {
	double L = GetValue(*ts, LGAL);
	double R = GetValue(*ts, LGAR);
	A += 2.0 * PI_VALUE * R * L;
      }
    }
    return A;
}

};

//Vertical distribution of  fip(ip) in segments The vector  v has Tree
//age  positions  representing  the  vertical  division  of  the  tree
//into intervals  of mean annual growth).
//Usage is simply:
//   pair<vector<pair<double,int> >,double> p;
//   Accumulate(tree,v,CollectVerticalDistributionOfFip());
//where the  p.first is the vector  for cumulative fips  and number of
//their observations and p.second the height interval When printing to
//the  file  remember  to  divide  cumulative "fips"  with  number  of
//observations.
class CollectVerticalDistributionOfFip{
public:
  pair<vector<pair<double,int> >,double>& 
  operator()(pair<vector<pair<double,int> >,double>& p,
	     TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)const
  {
    if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
      if (GetValue(*ts,LGAWf) > R_EPSILON){
	const ParametricCurve& fip = GetFunction(GetTree(*ts),LGMIP);
	double B = GetValue(GetTree(*ts),TreeQinMax);
	double qin = GetValue(*ts,LGAQin);
	double ip = qin/B;
	double result = fip(ip);
	double height = GetPoint(*ts).getZ();
	//Find the height category for the segment (p.second is the interval)
	double category = height/p.second;
	//Drop the decimal part, the integer part is the index for the
	//position in the vector v.
	unsigned int index = static_cast<unsigned int>(category);
	//cout << category << " " << index << " "  << height << " "  << p.second << endl;
	if (index < p.first.size()){
	  p.first[index].first += result;//cumulative f(ip) 
	  p.first[index].second += 1; //number of  observations (segments with
	                              //foliage)
	}
	else{
	  cout << "CollectVerticalDistributionOfFip: Segment height exceeds tree height" <<endl;
	  cout << "Ignoring " << height << " > "  << GetValue(GetTree(*ts),LGAH);
	}
      }
    }
    return p;
  }
};

class CollectSegmentDiameters{
public:
  list<double>& operator()(list<double>& ls,TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc)
  {
    if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
      double d = 2.0* GetValue(*ts,LGAR);
      ls.push_back(d);
    }
    return ls;
  }
};

class ResultsToFileInfo{
 public:
  int bs;                 // box STAR  0 = constant value,
                          //           1 = value of box
                          //           2 = value of box & mean direction effect
  int corr;               // If STAR_eq -> STAR correction considered (0 = no, 1 = yes)
  LGMdouble vox;          // Size of voxel (=length of side, assuming voxels are cubes)
  string location_file;
  string tree_file;
};

class ResultsToFile {
 public:
 ResultsToFile(const string& filename): fname(filename), qmax(1.0)
  {
    one_time = true;
    ofstream f(fname.c_str() , ofstream::trunc);
    f << "Height Dist_top Dist_stem Qin Qmax" << endl;
    f.close();
  }
 ResultsToFile(const string& filename,  ResultsToFileInfo in_fo): 
  fname(filename), info(in_fo)
  {
    one_time = false;
  }

  TreeCompartment<ScotsPineSegment,ScotsPineBud>*
    operator()(TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc) const
    {
      if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
	LGMdouble Wf = GetValue(*ts,LGAWf);
	if(Wf > 0.0) {
	  LGMdouble Qin = GetValue(*ts,LGAQin);
	  Tree<ScotsPineSegment,ScotsPineBud>& t = GetTree(*ts);
	  LGMdouble H = GetValue(t, LGAH);
	  Point pp = GetPoint(*ts);
	  LGMdouble ver_dist = H - pp.getZ();
	  Axis<ScotsPineSegment,ScotsPineBud>& ax = GetAxis(t);
	  Bud<ScotsPineSegment,ScotsPineBud>* lb = GetTerminatingBud(ax);
	  LGMdouble hor_dist = sqrt(pow(pp.getX()-GetPoint(*lb).getX(),2.0) +
				    pow(pp.getY()-GetPoint(*lb).getY(),2.0));

	  ofstream f(fname.c_str() , ofstream::app);
	  if(one_time) {
	    f << pp.getZ() << " " << ver_dist << " " << hor_dist << " " << Qin << " " << qmax << endl;
	  }
	  else {
	    f << info.location_file << " " << info.tree_file << " "
	      << info.vox << " " << info.bs << " " << info.corr << " "
	      << pp.getZ() << " " << ver_dist << " " << hor_dist << " " << Qin << " " << qmax << " " << endl;
	  }
	  f.close();
	}
      }
      return tc;
    }

  void setQmax(const LGMdouble qm) {qmax = qm;}

 private:
  string fname;
  bool one_time;
  ResultsToFileInfo info;
  LGMdouble qmax;
};


class ZeroWoodyRadius{
public:
  TreeCompartment<ScotsPineSegment,ScotsPineBud>* 
    operator()(TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc) const
  {
    if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
      SetValue(*ts,LGAR, 0.0001);
    }
    return tc;
  }
};

//Distributes trees TreeSegments that are in the box randomly                                                       

class RandomizeSegmentsInBox{
 public:
 RandomizeSegmentsInBox(const Point lli, const Point uri, ForestGap gi) :
  ll(lli), ur(uri), gap(gi), density(1,1,1), limits(1), small_box_size(3) {
    string str = "randomizeinbox.par";
    ifstream rf(str.c_str());
    if (!rf) {
      cout << "File randomizeinbox.par is missing with flag -randomizeInBox!"  << endl;
      exit(0);
    }
 
    string line;
    getline(rf,line);//The header
    getline(rf,line);//First data line
    istringstream iss(line);
    iss >> Nx >> Ny >> Nz;
    getline(rf,line);
    getline(rf,line);
    istringstream iss1(line);
    iss1 >> dens_min >> dens_max;
    rf.close();

    N_tot = Nx * Ny * Nz;
    small_box_size[0] = (ur.getX() - ll.getX()) / static_cast<double>(Nx);
    small_box_size[1] = (ur.getY() - ll.getY()) / static_cast<double>(Ny);
    small_box_size[2] = (ur.getZ() - ll.getZ()) / static_cast<double>(Nz);

    cout << Nx << " " << Ny << " " << Nz << endl;
    cout << dens_min << " " <<  dens_max << endl;

    //Now set the densities of the canopy space compartments
    if(dens_max <= dens_min) {
      cout << "-randomizeInBox dens_max <= dens_min !" << endl;
      exit(0);
    }

    density.resize(Nx, Ny, Nz);
    double dens_sum = 0.0;
    for(int i = 0; i < Nx; i++) {
      for(int j = 0; j < Ny; j++) {
	for(int k = 0; k < Nz; k++) {
	  density[i][j][k] = dens_min + ran3(&ran3_seed)*(dens_max - dens_min);
	  dens_sum += density[i][j][k];
	}
      }
    }
    for(int i = 0; i < Nx; i++) {
      for(int j = 0; j < Ny; j++) {
	for(int k = 0; k < Nz; k++) {
	  density[i][j][k] /= dens_sum;
	}
      }
    }

    limits.resize(N_tot);
    dens_sum = 0.0;
    vector<int> inds(3);
    int index = 0;
    for(int i = 0; i < Nx; i++) {
      for(int j = 0; j < Ny; j++) {
	for(int k = 0; k < Nz; k++) {
	  dens_sum += density[i][j][k];
	  limits[index] = dens_sum;
	  index++;
	}
      }
    }
  }


  void box_ind(const int n, const int ny, const int nz, vector<int>& out) const {
    div_t dr1, dr2;
    dr1 = div(n,ny*nz);
    dr2 = div(dr1.rem,nz);
    out[0] = dr1.quot;
    out[1] = dr2.quot;
    out[2] = dr2.rem;
    return;
  }

  int box_no(const vector<double>& lims, const double value) const {
    int ind = 0;
    while(lims[ind] < value) {
      ind++;
    }
    return ind;
  }


  bool inGap(const Point p, const ForestGap gp) const {
    return sqrt(pow(p.getX() - (gp.first).first, 2.0) + pow(p.getY() - (gp.first).second, 2.0))
      < gp.second;
  }

  bool inBigBox(const Point& p, const Point& ll, const Point& ur) const {
    bool result = false;
    if(p.getX() < ll.getX()) return result;
    if(p.getX() > ur.getX()) return result;
    if(p.getY() < ll.getY()) return result;
    if(p.getY() > ur.getY()) return result;
    if(p.getZ() < ll.getZ()) return result;
    if(p.getZ() > ur.getZ()) return result;

    result = true;
    return result;
  }

  TreeCompartment<ScotsPineSegment,ScotsPineBud>*
    operator()(TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc) const
    {
      if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
	//position in the small box
	double l = GetValue(*ts, LGAL);

	bool out_of_space = true;
	bool in_gap = true;
	int no_loop = 0;
	Point p;
	while((in_gap || out_of_space) && (no_loop < 30)) {

	  //first the small box
	  double ran = ran3(&ran3_seed);
	  int s_box = box_no(limits, ran);
	  vector<int> sb_ind(3);
	  box_ind(s_box, Ny, Nz, sb_ind);


	  Point sb_ll(ll.getX()+static_cast<double>(sb_ind[0])*small_box_size[0],
		      ll.getY()+static_cast<double>(sb_ind[1])*small_box_size[1],
		      ll.getZ()+static_cast<double>(sb_ind[2])*small_box_size[2]);
	  int no_loop = 0;
	  Point p;
	  while((in_gap || out_of_space) && (no_loop < 20)) {
	    LGMdouble x = sb_ll.getX() + ran3(&ran3_seed)*small_box_size[0];
	    LGMdouble y = sb_ll.getY() + ran3(&ran3_seed)*small_box_size[1];
	    LGMdouble z = sb_ll.getZ() + ran3(&ran3_seed)*small_box_size[2];

	    p = Point(x,y,z);

	    if(!inGap(p, gap)) {
	      in_gap = false;
	    }
	  
	    if ( (inBigBox(p - (Point)(0.5*l*GetDirection(*tc)), ll, ur)) &&
		 (inBigBox(p + (Point)(0.5*l*GetDirection(*tc)), ll, ur)) ) {
	      out_of_space = false;
	    }
	    no_loop++;                //this prevents infinite loop
	  }

	  SetPoint(*ts, p - (Point)(0.5*l*GetDirection(*tc)));
	}
      }
      return tc;
    }
 private:
  Point ll;         //lower left                                                                                  
  Point ur;         //upper right                                                                                 
  ForestGap gap;
  int Nx, Ny, Nz;  //division of forest canopy box for variable density
  int N_tot;        // = Nx*Ny*Nz
  double dens_min, dens_max;  //variation of density
  TMatrix3D<double> density;
  vector<double> limits;
  vector<double> small_box_size;   //x, y, and z side lengths of small boxes
};


  class PrintPositions{
  public:
    PrintPositions (ofstream& posfile, const int seedi): pfile(posfile), sd(seedi) {}
    TreeCompartment<ScotsPineSegment,ScotsPineBud>* operator()
      (TreeCompartment<ScotsPineSegment,ScotsPineBud>* tc) const
      {
	if (ScotsPineSegment* ts = dynamic_cast<ScotsPineSegment*>(tc)){
	  if(GetValue(*ts, LGAWf) > 0.0) {
	    Point p = GetMidPoint(*ts);
	    pfile << sd << " " << p.getX() << " " << p.getY() << " " << p.getZ() << endl;
	  }
	}
	return tc;
      }

  private:
    ofstream& pfile;
    int sd;
  };


#endif
 

