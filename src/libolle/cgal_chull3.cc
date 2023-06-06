#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <vector>

// ok, no idea what any of this means.
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedron_3<K> Polyhedron_3;
typedef K::Point_3 Point_3;
// just for debugging
using std::cout;
using std::endl;

// this makes it a C wrapper
extern "C"{
    // calculate the convex hull of a set of 2D points
    void lo_cgal_chull3(const int &np,const double* r,int &nphull, double *&rhull){

        // Stuff the points into an iterator, whatever that is.
        std::vector<Point_3> points;
        // Stuff the points into this vector
        int l=0;
        for(int j=0;j<np;j++){
            points.push_back( Point_3(r[l+0],r[l+1],r[l+2]) );
            l+=3;
        }
        // define polyhedron to hold convex hull
        Polyhedron_3 poly;

        // compute convex hull of non-collinear points
        CGAL::convex_hull_3(points.begin(), points.end(), poly);
        // Get the number of point on the hull.
        nphull=poly.size_of_vertices();
        // Make some space in a normal array for the result
        rhull= new double[nphull*3];
        // Iterate through the polyhedron somehow.
        l=-1;
        for(Polyhedron_3::Vertex_iterator vert=poly.vertices_begin(); vert!=poly.vertices_end(); vert++){
            l++; rhull[l]=vert->point().x();
            l++; rhull[l]=vert->point().y();
            l++; rhull[l]=vert->point().z();
        }
        // Some cleanup perhaps
        //poly.destroy();
        // delete points
        // delete poly
    }
}

