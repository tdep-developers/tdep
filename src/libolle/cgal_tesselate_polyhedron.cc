// include things for the polyhedron example
#include <vector>
#include <map>
#include <iostream>
#include <boost/iterator.hpp> // thing I had to add with cgal 4.11.2

// minimal thing
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
// to get the convex hull (can not use standard polyhedron)
#include <CGAL/convex_hull_3.h>
#include <CGAL/Mesh_polyhedron_3.h>

// for a domain with features
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>

// for the triangulating the domain
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

// New thingy?
typedef CGAL::Sequential_tag Concurrency_tag;

// really basic
typedef CGAL::Exact_predicates_inexact_constructions_kernel   Kernel;
typedef Kernel::Point_3                                       Point_3;

// this was also in the polyhderon with features example
typedef CGAL::Mesh_polyhedron_3<Kernel>::type                 Polyhedron;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<Kernel>  Mesh_domain;

// This is some strange stuff from the CGAL example. I do not understant at all.
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr,Mesh_domain::Corner_index,Mesh_domain::Curve_segment_index> C3t3;

// Criteria?
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;
// // To avoid verbose function and named parameters call
using namespace CGAL::parameters;

using std::cout;
using std::endl;

extern "C"{
    // Create a mesh of a convex polyhedron defined as the convex hull of the points r.
    void lo_cgal_tesselate_polyhedron(const int &np, const double* r, int &nvertices,
        double *&vertices, int &ntet, int *&tetrahedrons,
        const double &p_max_dihedral_angle,
        const double &p_edge_size,
        const double &p_facet_angle,
        const double &p_facet_size,
        const double &p_facet_distance,
        const double &p_cell_radius_edge_ratio,
        const double &p_cell_size
        ){

        // store the points in a vector of CGAL points
        std::vector<Point_3> vpoints;
        int l=0;
        for(int i=0;i<np*3;i+=3){
            vpoints.push_back( Point_3( r[i+0], r[i+1], r[i+2] ) );
        }
        // Get the polyhedron with a complex hull algorithm
        Polyhedron P;
        CGAL::convex_hull_3(vpoints.begin(), vpoints.end(), P);

        /*
            If I am not mistaken, I have now built the polyhedron in a way that CGAL likes.
            Then I have just adapted the provided example and use that to tesselate it.
        */

        // Create domain (no idea what that means)
        Mesh_domain domain(P);
        // Detect the sharp features of the domain. 190 is the max dihedral angle for detection.
        // Since I have trivial polygons, it's better to detect all of them.
        domain.detect_features(p_max_dihedral_angle);
        // Mesh criteria (I have to figure out what this means, and pass sensible settings)
        Mesh_criteria criteria(
            edge_size              = p_edge_size,
            facet_angle            = p_facet_angle,
            facet_size             = p_facet_size,
            facet_distance         = p_facet_distance,
            cell_radius_edge_ratio = p_cell_radius_edge_ratio,
            cell_size              = p_cell_size
        );
        //Mesh_criteria criteria(
        //    facet_angle            = p_facet_angle,
        //    facet_size             = p_facet_size,
        //    facet_distance         = p_facet_distance,
        //    cell_radius_edge_ratio = p_cell_radius_edge_ratio,
        //    cell_size              = p_cell_size
        //);

        // Generate the actual mesh
        //C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria);
        C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                            lloyd(time_limit=30,do_freeze=false,convergence=0.00001),
                                            exude(time_limit=30));
        // Generate the mesh in a cooler way!

        /*
            Ok. I got a mesh. Now I have to dump it to fortran somehow. First I have to figure
            how many vertices there are. It seems I extract the triangulation? Or do I just alias
            it to some other name? Not sure, but it seems to work.
        */

        const C3t3::Triangulation &tr=c3t3.triangulation();
        // now I can figure out how many vertices there are
        nvertices=tr.number_of_vertices();
        vertices = new double[nvertices*3];

        /*
            I iterate over all the vertices, store the coordinates in my array and also generate
            a map so that I can get the indices from the triangulation structure later
        */

        // the map to indices
        std::map<C3t3::Triangulation::Vertex_handle, int> V;
        // the iterator
        C3t3::Triangulation::Finite_vertices_iterator vt;
        l=0;
        int ll=-1;
        for(vt = tr.finite_vertices_begin(); vt!=tr.finite_vertices_end(); vt++){
            // I think this creates a map to simple integers. I index the fortran way
            l++;
            V[vt]=l;
            ll++; vertices[ll]=vt->point().x();
            ll++; vertices[ll]=vt->point().y();
            ll++; vertices[ll]=vt->point().z();
        }

        // I should know how many tetrahedrons there are
        ntet=c3t3.number_of_cells_in_complex();
        tetrahedrons = new int[ntet*4];

        // With the vertices labelled, I can get the tetrahedrons.
        C3t3::Cells_in_complex_iterator it;
        l=-1;
        for(it = c3t3.cells_in_complex_begin(); it!=c3t3.cells_in_complex_end(); it++){
            for(int i=0; i<4; i++){
                l++;
                tetrahedrons[l]=V[it->vertex(i)];
            }
        }
        // Output the cgal way. To check that I passed things correctly to fortran.
        // std::ofstream medit_file("out.mesh");
        // c3t3.output_to_medit(medit_file);
    }
}


//  // No idea how this works. I hope it builds a polyhedron from my definition of one.
//  template<class HDS>
//  class polyhedron_builder : public CGAL::Modifier_base<HDS> {
//      public:
//      std::vector<double> &coords;
//      std::vector<int>    &tris;
//      const int           &maxnf; // den lade jag till!
//
//      //polyhedron_builder( std::vector<double> &_coords, std::vector<int> &_tris ) : coords(_coords), tris(_tris) {}
//      polyhedron_builder( std::vector<double> &_coords, std::vector<int> &_tris,const int &_maxnf ) : coords(_coords), tris(_tris), maxnf(_maxnf) {}
//
//      void operator()( HDS& hds) {
//          typedef typename HDS::Vertex   Vertex;
//          typedef typename Vertex::Point Point;
//
//          // create a cgal incremental builder
//          CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
//          //B.begin_surface( coords.size()/3, tris.size()/3 );
//          B.begin_surface( coords.size()/3, tris.size()/maxnf );
//
//          // add the polyhedron vertices. Ok I think.
//          for( int i=0; i<(int)coords.size(); i+=3 ){
//              B.add_vertex( Point( coords[i+0], coords[i+1], coords[i+2] ) );
//          }
//
//          // my version of adding faces
//          for( int i=0; i<(int)tris.size(); i+=maxnf ){
//              B.begin_facet();
//              // loop over correct number of vertices
//              for( int j=0; j<maxnf; j++){
//                  if( tris[i+j]>-1 ){
//                      B.add_vertex_to_facet( tris[i+j] );
//                  }
//              }
//              B.end_facet();
//          }
//          // finish up the surface
//          B.end_surface();
//      }
//  };
