
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <vector>
#include <map>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_3<unsigned int, K> Vb;
typedef CGAL::Triangulation_data_structure_3<Vb> Tds;
typedef CGAL::Delaunay_triangulation_3<K, Tds> Delaunay;
typedef CGAL::Tetrahedron_3<K> Tetrahedron;
typedef K::Point_3 Point_3;

extern "C"{
    // Calculate the Delaunay triangulation of a set of 3D points.
    // returns an array of triangles, defined as a triplet of indices to the original points
    void lo_cgal_deltri3(const int &np, const double *r, int &ntet, int *&tet, double *&tetvol){
        // I think this is a vector with points and associated indices.
        std::vector< std::pair<Point_3,unsigned> > points;
        // Stuff the points into this vector
        int l=0;
        for(int j=0;j<np;j++){
            points.push_back( std::make_pair( Point_3(r[l+0],r[l+1],r[l+2]), j ) );
            l+=3;
        }
        // make the actual triangulation
        Delaunay dt;
        dt.insert(points.begin(), points.end());
        // get number of tetrahedrons
        ntet=dt.number_of_finite_cells();

        // Return the actual tetrahedrons, as well as their volumes
        tet = new int[ntet*4];
        tetvol = new double[ntet];
        int ll=-1;
        l=-1;
        for(Delaunay::Finite_cells_iterator fit = dt.finite_cells_begin(); fit != dt.finite_cells_end(); ++fit) {
            for(int j=0;j<4;j++){
                l++;
                tet[l]=fit->vertex(j)->info()+1;
            }
            Tetrahedron tetr = dt.tetrahedron( fit );
            ll++;
            tetvol[ll]=tetr.volume();
        }
        // Find the adjacent tetrahedrons

    }
}
