#include <iostream>
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Homogeneous.h>
//#include <CGAL/Exact_integer.h>

//#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Simple_homogeneous.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Gmpq.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <vector>

// ok, no idea what any of this means.
//typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
//typedef CGAL::Homogeneous<CGAL::Exact_integer> K;
//typedef CGAL::Exact_predicates_exact_constructions_kernel K;
//typedef CGAL::Simple_cartesian<double> K;
//typedef CGAL::Gmpz  NT;
//typedef CGAL::Homogeneous<NT>  K;
typedef CGAL::Simple_cartesian<CGAL::Gmpq> K;
typedef CGAL::Polyhedron_3<K>                             Polyhedron_3;
typedef CGAL::Nef_polyhedron_3<K>                         Nef_Polyhedron_3;
typedef CGAL::Triangulation_3<K>                          Triangulation;
typedef CGAL::Tetrahedron_3<K>                            Tetrahedron;
typedef K::Point_3 Point_3;
// just for debugging
using std::cout;
using std::endl;

// this makes it a C wrapper
extern "C"{
    // The intersection between two convex hulls
    void lo_cgal_chull3_intersection(const int &np1,const int &np2, const double* r1, const double* r2, int &nphull, double *&rhull, double &vol){
        // Stuff the points into something that CGAL likes 
        std::vector<Point_3> points1;
        std::vector<Point_3> points2;
        // Stuff the points into this vector
        int l=0;
        for(int j=0;j<np1;j++){
            points1.push_back( Point_3(r1[l+0],r1[l+1],r1[l+2]) );
            l+=3;
        }
        l=0;
        for(int j=0;j<np2;j++){
            points2.push_back( Point_3(r2[l+0],r2[l+1],r2[l+2]) );
            l+=3;
        }
        // define polyhedron to hold convex hull
        Polyhedron_3 poly1;
        Polyhedron_3 poly2;
        // compute convex hull to get two polyhedrons
        CGAL::convex_hull_3(points1.begin(), points1.end(), poly1);
        CGAL::convex_hull_3(points2.begin(), points2.end(), poly2);
        // probably some sanity checks would be good
        //if( poly1.is_closed() != 1 ){}
        //if( poly2.is_closed() != 1 ){}
        // these should be converted to NEF polyhedrons. 
        Nef_Polyhedron_3 nefpoly1(poly1);
        Nef_Polyhedron_3 nefpoly2(poly2);
        // can it really be this simple?
        Nef_Polyhedron_3 intersection=nefpoly1*nefpoly2;
        // stuff it into a triangulation
        Triangulation tri;
        if ( intersection.is_simple() ){
            // Stuff the points
            for(Nef_Polyhedron_3::Vertex_const_iterator v = intersection.vertices_begin();
                v != intersection.vertices_end(); ++v){
                tri.insert(v->point());
            }
        } else {
            // probably something weird
            nphull=0;
        }
        // number of points
        nphull=tri.number_of_vertices();
        // if there are some, get the points and the volume
        if ( nphull > 0 ){
            // store the points
            rhull=new double[nphull*3];
            l=-1;
            for(Triangulation::Finite_vertices_iterator v=tri.finite_vertices_begin(); 
                v!=tri.finite_vertices_end(); ++v){
                l++; rhull[l]=CGAL::to_double(v->point().x());
                l++; rhull[l]=CGAL::to_double(v->point().y());
                l++; rhull[l]=CGAL::to_double(v->point().z());
            }
            // iterate over all cells
            vol=0;
            for(Triangulation::Finite_cells_iterator t=tri.finite_cells_begin();
                t != tri.finite_cells_end(); ++t){
                // yoink out the tetrahedron and measure the volume
                Tetrahedron tetr = tri.tetrahedron(t);
                vol += CGAL::to_double(tetr.volume());
                //
            }
        }
        // and by now it should be done
    }
}

