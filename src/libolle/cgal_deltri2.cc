
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, K> Vb;
typedef CGAL::Triangulation_data_structure_2<Vb> Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds> Delaunay;
typedef CGAL::Triangle_2<K> Triangle;
typedef K::Point_2 Point_2;

using std::cout;
using std::endl;

extern "C"{
    // Calculate the Delaunay triangulation of a set of 2D points.
    // returns an array of triangles, defined as indices
    void lo_cgal_deltri2(const int &np, const double *r, int &ntri, int *&tri, double *&triarea){
        // I think this is a vector with points and associated indices.
        std::vector< std::pair<Point_2,unsigned> > points;
        int l=0;
        for(int j=0;j<np;j++){
            points.push_back( std::make_pair( Point_2(r[l+0],r[l+1]), j ) );
            l=l+2;
        }
        // Do the actual triangulation
        Delaunay dt;
        dt.insert(points.begin(), points.end());
        // Return the number of triangles
        ntri=dt.number_of_faces();
        // Return the actual triangles.
        tri = new int[ntri*3];
        triarea = new double[ntri];
        l=-1;
        int ll=-1;
        for(Delaunay::Finite_faces_iterator fit = dt.finite_faces_begin(); fit != dt.finite_faces_end(); ++fit) {
            for(int j=0;j<3;j++){
                l++;
                tri[l]=fit->vertex(j)->info()+1;
            }
            Triangle triangle = dt.triangle( fit );
            ll++;
            triarea[ll]=triangle.area();
        }
        // and now it should all be done!
    }
}

