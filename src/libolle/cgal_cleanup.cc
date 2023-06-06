
// these should be in a generic interface. I wonder what that is called in c++
extern "C"{ void lo_cgal_cleanup_int_pointer(int *&A){ delete A; }}
extern "C"{ void lo_cgal_cleanup_double_pointer(double *&A){ delete A; }}

