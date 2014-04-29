#ifndef CALCULATE_LIGHT_H
#define CALCULATE_LIGHT_H
#include <Lignum.h>
#include <ScotsPine.h>
#include <CopyDumpCfTree.h>
#include <VoxelSpace.h>
#include <SomeFunctors.h>
#include <vector>
#include <utility>
#include <BorderForest.h>
using namespace std;

using namespace Lignum;
using namespace sky;


#define HIT_THE_FOLIAGE 1
#define NO_HIT 0
#define HIT_THE_WOOD -1


// Functions in this file are to be used instead of radiation calculations in stl-lignum
// (stl-lignum/include/Shading.h * and ...I.h). Replace them eventually with this.
// _1 addition to the names are to distinquish the from functions of stl-lignum radiation
// calculations.

// EvaluateRadiationForCfTreeSegment_1 evaluates shading
// caused by all other segments on this conifer segment.  The shading caused by segments in the crown
// of tree itself is evaluated by ShadingEffectOfCfTreeSegment_1<TS,BUD> (it is the same as
// ShadingEffectOfCfTreeSegment in stl-lignum; the functor is duplicated here only for convenience).
// After that other trees are accounted for with voxelspace (voxel_space->getRoute() etc) and surrounding
// stand with border_forest->getBorderForestExtinction().
// This function does not evaluate the shading by surrounding trees (Qin_stand) and border forest separately 
// (slighthly faster);
// EvaluateRadiationForCfTreeSegment_1s (defined in ../GrowthLoopRadiationI.h) does that.

//If the attributes voxel_space and border_forest are set (then
//voxel_space != NULL), the attenuation of the beam in the voxel_space
//and border_forest is taken into consideration.

template <class TS, class BUD>
class EvaluateRadiationForCfTreeSegment_1 {
public:
    EvaluateRadiationForCfTreeSegment_1(const ParametricCurve& k) : K(k),
        evaluate_voxel_and_border(false) {}
    EvaluateRadiationForCfTreeSegment_1(const ParametricCurve& k,
                                        VoxelSpace* vs,BorderForest* bf, LGMdouble gk):
        K(k),voxel_space(vs), border_forest(bf), evaluate_voxel_and_border(true),
        green_voxel_ext_coeff(gk) {}

    TreeCompartment<TS,BUD>* operator()(TreeCompartment<TS,BUD>* tc)const;
private:
    const ParametricCurve& K;
    VoxelSpace* voxel_space;
    BorderForest* border_forest;
    bool evaluate_voxel_and_border;
    LGMdouble green_voxel_ext_coeff;
};


//This functor ShadingEffectOfCfTreeSegment<TS,BUD> evaluates shading caused
//by a conifer segment on this conifer segment (shaded_s)

template <class TS,class BUD>
class ShadingEffectOfCfTreeSegment_1 {
public:
    ShadingEffectOfCfTreeSegment_1(CfTreeSegment<TS,BUD>* ts, const ParametricCurve& K_in,
                                   vector<double>& sectors)
        :shaded_s(ts), K(K_in),S(sectors){}
    //ForEach functor to compute shadiness
    TreeCompartment<TS,BUD>*  operator()(TreeCompartment<TS,BUD>* tc)const;
    //Get vector for S (shadiness)
    vector<double>& getS(){return S;}
private:
    CfTreeSegment<TS,BUD>* shaded_s;
    //Avoid unnecessary constructor calls in generic algorithms
    const ParametricCurve& K;
    vector<double>& S;
};


//=======================================================================================================
//This version of radiation evaluates radiation conditions for subject tree by pairwise
// comparison

// EvaluateRadiationForCfTreeSegment_2 evaluates shading by all other segments by paiwise comparison (segments
// in own crown & other trees). It uses ShadingEffectOfCfTreeSegment_1<TS,BUD> to evaluate shading.
// This functor evaluates shading by own crown and shading by other trees (stand) separately and updates
// Qin_stand in TreeSegment.


//This functor EvaluateRadiationForCfTreeSegment evaluates shading
//caused by all other segments on this conifer segment. This functor
//uses functor ShadingEffectOfCfTreeSegment<TS,BUD> to go through all
//segments to check the shading.

//If the attributes voxel_space and border_forest are set (then
//voxel_space != NULL), the attenuation of the bean in the voxel_space
//and border_forest is taken into consideration.

template <class TS, class BUD, class TREE>
class EvaluateRadiationForCfTreeSegment_2 {
public:
    EvaluateRadiationForCfTreeSegment_2(const ParametricCurve& k, vector<TREE*>& vt):
        K(k), vtree(vt) {only_self = false;}
    EvaluateRadiationForCfTreeSegment_2(const ParametricCurve& k, vector<TREE*>& vt, bool o_s):
        K(k), vtree(vt), only_self(o_s) {}

    TreeCompartment<TS,BUD>* operator()(TreeCompartment<TS,BUD>* tc)const;
private:
    const ParametricCurve& K;
    vector<TREE*>& vtree;
    bool only_self;
};


//===============================================================================
// VoxelSpace calculation

template <class TS, class BUD>
class EvaluateRadiationForCfTreeSegment_3 {
public:
    EvaluateRadiationForCfTreeSegment_3(const ParametricCurve& k) : K(k),
        evaluate_border(false) {}
    EvaluateRadiationForCfTreeSegment_3(const ParametricCurve& k,
                                        VoxelSpace* vs,BorderForest* bf, bool border,
                                        LGMdouble a, LGMdouble b, bool sd,LGMdouble kbc,
                                        bool d_e, bool wd, bool cs, LGMdouble st,bool calculateDirectionalStar):
        K(k),voxel_space(vs), border_forest(bf), evaluate_border(border),
        par_a(a), par_b(b), dump_self(sd), k_border_conifer(kbc), dir_effect(d_e),
        wood(wd), correct_star(cs), constant_star(st),calculateDirectionalStar(calculateDirectionalStar) {}

    TreeCompartment<TS,BUD>* operator()(TreeCompartment<TS,BUD>* tc)const;

private:
    const ParametricCurve& K;
    VoxelSpace* voxel_space;
    BorderForest* border_forest;
    bool evaluate_border;
    LGMdouble par_a, par_b;
    bool dump_self;
    LGMdouble k_border_conifer;
    bool dir_effect;   //If effect of mean segment direction in boxes considered
    bool wood;         //If the needleless woody parts are included into calculation
    LGMdouble constant_star; //If this > 0, k = constant_star else k = star_mean
    bool correct_star;       //If star_eq -> star correction is done
    bool calculateDirectionalStar;



};

class AccumulateOpticalDepth{
public:
    AccumulateOpticalDepth(LGMdouble side, LGMdouble a, LGMdouble b, Point loc, ParametricCurve kk,
                           bool d_e, bool wd, bool cs, LGMdouble st,bool calculateDirectionalStar) :
        box_side_length(side), par_a(a), par_b(b), seg_loc(loc), K(kk), dir_effect(d_e), wood(wd), calculateDirectionalStar(calculateDirectionalStar),
        constant_star(st),correct_star(cs) {box_volume = pow(box_side_length,3.0);}
    double operator()(double o_d,VoxelMovement& vm){
        //    if((vm.af > R_EPSILON ||(wood && vm.wood_area > R_EPSILON)) && vm.n_segs_real > 0.0) {


        if(vm.af > R_EPSILON ||(wood && vm.wood_area > R_EPSILON)) {
            LGMdouble k;
            vector<LGMdouble>  kdir;

            if(calculateDirectionalStar){
                kdir = vm.starDir;


            }
            else{
                if(constant_star > 0.0)
                    k = constant_star;
                else
                    k = vm.STAR_mean;

            }

            if(correct_star) {
                k = max(0.0,-0.014+1.056*k);
            }

            //NOTE: here transformation STAR_eq --> STAR; documented in
            //~/Riston-D/E/LIGNUM/Light/summer-09-test/STAR-vs-STAR_eq.pdf

            //Effect of hit angle to the mean direction of shoots in the voxel box, documented in
            //~/Riston-D/E/LIGNUM/Light/Article/vs-STARmean-all.pdf and
            //~/Riston-D/E/LIGNUM/Light/Article/vs-STARmean-approximation.pdf
            PositionVector mean_dir = vm.mean_direction;
            LGMdouble mean_dir_length = mean_dir.length();
            LGMdouble effect = 1.0;
            if(dir_effect) {
                if(mean_dir_length > 0.0){
                    mean_dir.normalize();
                    LGMdouble inclination  =  PI_DIV_2 - acos(fabs(Dot(mean_dir,beam_dir)));
                    effect =  K(inclination)/K(0.7);
                }
            }
            //this scales the effect depending on how parallel the segments are
            /*       if(mean_dir_length > 0.0) { */
            /* 	mean_dir.normalize(); */
            /* 	LGMdouble inclination  =  PI_DIV_2 - acos(fabs(Dot(mean_dir,beam_dir))); */
            /* 	//	LGMdouble effect  = 1.13 - 0.24*pow(inclination,2.0); */
            /* 	//	LGMdouble effect = 1.2-0.3*inclination*(inclination+0.3); */
            /* 	LGMdouble u =  mean_dir_length/vm.n_segs_real; */
            /* 	effect = 1.0 - u + u * K(inclination)/K(0.7); */
            /* 	cout << " u " << u << endl; */
            /*       } */
            if(calculateDirectionalStar){
                LGMdouble secondfactor = effect * vm.af * vm.l / box_volume;
                std::transform(kdir.begin(),kdir.end(),kdir.begin(),std::bind1st(std::multiplies<LGMdouble>(),secondfactor));
                o_d  = std::accumulate(kdir.begin(),kdir.end(),0);

            }
            else{

                o_d += effect * k * vm.af * vm.l / box_volume;
            }


            if(wood) {
                //Mean projection area of surface of a circular cylinder (excluding end disks)
                // is 1/4 of its area
                o_d += 0.25 * vm.wood_area * vm.l / box_volume;
            }
            cout.precision(15);
            cout<<"this is the o_d value"<<o_d<<endl;
        }
        return o_d;
    }

    //This is to take care of the direction of the beam of radiation
    PositionVector beam_dir;
private:
    LGMdouble box_side_length;
    LGMdouble box_volume;
    LGMdouble par_a, par_b;
    Point seg_loc;   //location of segment
    ParametricCurve K;
    bool dir_effect;         //If direction effect of segments in box considered
    bool wood;               //If woody parts are considered
    LGMdouble constant_star; //If this > 0, k = constant_star else k = star_mean
    bool correct_star;       //If star_eq -> star correction is done
    bool calculateDirectionalStar;           // Directional Star values are required this is used
};


#undef HIT_THE_FOLIAGE
#undef NO_HIT
#undef HIT_THE_WOOD

#include <CalculateLightI.h>

#endif
