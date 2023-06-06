 
#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/convex_hull_2.h>

// ok, no idea what any of this means.
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
// just for debugging
using std::cout;
using std::endl;

// this makes it a C wrapper
extern "C"{
    // calculate the convex hull of a set of 2D points
    void lo_cgal_chull2(const int &np,const double* r,int &nphull, double *&rhull){

        // Now I have to stuff these points into something that CGAL likes
        Point_2 *points_in_cgal_format;
        points_in_cgal_format = new Point_2[np]; 
        int l=0;
        for(int j=0;j<np;j++){
            points_in_cgal_format[j]=Point_2(r[l+0],r[l+1]);
            l=l+2;
        }

        // Not bad! First, make some space for the result:
        Point_2 *result;
        result = new Point_2[np];
        // And then call CGAL and let it do its thing
        Point_2 *ptr = CGAL::convex_hull_2( points_in_cgal_format, points_in_cgal_format+np, result );

        // This is how you figure out how many points there are on the hull.
        nphull=ptr-result;
        // Now I have to stuff these coordinates into my output array to send back to Fortran
        rhull= new double[nphull*2];
        l=-1;
        for(int j=0;j<nphull;j++){
            for(int k=0;k<2;k++){
                l=l+1;
                rhull[l]=result[j][k];
            }
        }
        // this is my feeble attempt at cleanup, seems to work
        delete points_in_cgal_format;        
        delete result;
    }
}

