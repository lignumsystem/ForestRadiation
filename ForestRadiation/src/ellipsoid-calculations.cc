#include <CalculateLight.h>

//Ellipsoid


void ellipsoid_interception(const Point& p00, const vector<vector<LGMdouble> >& directions,
			    const vector<vector<LGMdouble> >& ellipsoids, vector<int>& ellipsoid_hits) {

  int number_of_directions = (int)directions.size();
  for(int i = 0; i < number_of_directions; i++) {
    vector<LGMdouble> r_dir = directions[i];    //Direction of line-of-sight
    vector<vector<LGMdouble> >::const_iterator E;

    for(E = ellipsoids.begin(); E != ellipsoids.end(); E++) {
      Point pe((*E)[0],(*E)[1],(*E)[2]);       //Center of ellipsoid
      Point p0 = p00 - pe;        //Starting point of line-of-sight in ellipsoid center coordinates
      LGMdouble rh = (*E)[3];     //Horizontal semi-axis of the ellipsoid
      LGMdouble rv = (*E)[4];     //Vertical semi-axis of the ellipsoid

      LGMdouble a = pow(r_dir[0]/rh,2.0) + pow(r_dir[1]/rh,2.0) + pow(r_dir[2]/rv,2.0);
      LGMdouble rh2 = pow(rh,2.0);
      LGMdouble b = 2.0*p0.getX()*r_dir[0]/rh2 + 2.0*p0.getY()*r_dir[1]/rh2
	          + 2.0*p0.getZ()*r_dir[2]/pow(rv,2.0);
      LGMdouble c = pow(p0.getX()/rh, 2.0) + pow(p0.getY()/rh,2.0)
	          + pow(p0.getZ()/rv,2.0) - 1.0;; 
      LGMdouble diskr = pow(b, 2.0) - 4.0*a*c;
      if(diskr >= 0.0) {
	(ellipsoid_hits[i])++;
      }
    }  //End of ellipsoids
 
  }

  return;
}
