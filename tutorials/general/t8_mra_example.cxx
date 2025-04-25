#include <t8.h>                                     /* General t8code header, always include this. */
#include <t8_cmesh.hxx>                             /* cmesh definition and basic interface. */
#include <t8_forest/t8_forest_general.h>            /* forest definition and basic interface. */
#include <t8_vtk/t8_vtk_writer.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx> /* linear geometry of the cmesh */
#include <t8_forest/t8_forest_io.h>                 /* forest io interface. */
#include <t8_schemes/t8_default/t8_default.hxx> /* default refinement scheme. */
#include <t8_forest/t8_forest_geometrical.h>        /* geometrical information */
#include "t8_mra/vecmat.hxx"
#include "t8_mra/basis_functions.hxx"
#include "t8_mra/mask_coefficients.hxx"
#include "t8_mra/dunavant.hxx"
#include "t8_mra/t8_gsl.hxx"
#include "t8_mra/t8_mra.hxx"
#include "t8_mra/t8_unordered_dense.hxx"
#include "t8_unordered_dense.hxx"
#include <cmath>
#include <vector>
#include <sc_statistics.h>
#include <t8_refcount.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_profiling.h>
#include <t8_forest/t8_forest_private.h>
#include <t8_forest/t8_forest_types.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_forest/t8_forest_ghost.h>
#include <t8_forest/t8_forest_adapt.h>
#include <t8_forest/t8_forest_balance.h>
#include <t8_cmesh/t8_cmesh_offset.h>
#include <t8_cmesh/t8_cmesh_trees.h>
//#include <t8_element_c_interface.h>
#include <iostream>
#include <t8_eclass.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri_bits.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_dtri_connectivity.h>
#include <t8_schemes/t8_default/t8_default_tri/t8_default_tri.hxx>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
// #include <gsl/gsl_math.h>
// #include <gsl/gsl_interp2d.h>
// #include <gsl/gsl_spline2d.h>
#include <fstream>
#include <time.h>
// 3. Declare function prototypes or class declarations (if needed)

// 4. Optionally, include any namespaces you want to use (avoid if using standard practice)
using namespace std;

/* We build 4 kinds of cmeshes in the following by hand. */

/* This cmesh is [0,1]² with 8 triangles. */
t8_cmesh_t
t8_cmesh_new_basic (sc_MPI_Comm comm)
{

  /* 1. Defining an array with all vertices */
  /* Just all vertices of all trees. partly duplicated */
  double vertices[72] = { 0,0,0,0.5,0,0,0.5,0.5,0, //triangle 2
                          0,0,0,0,0.5,0,0.5,0.5,0, //triangle 1
                          0.5,0,0,1,0,0,1,0.5,0,   //triangle 3
                          0.5,0,0,0.5,0.5,0,1,0.5,0, //triangle 4
                          0,0.5,0,0.5,0.5,0,0.5,1,0, //triangle 5
                          0,0.5,0,0,1,0,0.5,1,0,     //triangle 6
                          0.5,0.5,0,1,0.5,0,1,1,0,   //triangle 7
                          0.5,0.5,0,0.5,1,0,1,1,0,   //triangle 8
  };

  /* 2. Initialization of the mesh */
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);

  /* 3. Definition of the geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);

  /* 4. Definition of the classes of the different trees */
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 3, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 4, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 5, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 6, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 7, T8_ECLASS_TRIANGLE);


  /* 5. Classification of the vertices for each tree */
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 3);
  t8_cmesh_set_tree_vertices (cmesh, 1, vertices + 9, 3);
  t8_cmesh_set_tree_vertices (cmesh, 2, vertices + 18, 3);
  t8_cmesh_set_tree_vertices (cmesh, 3, vertices + 27, 3);
  t8_cmesh_set_tree_vertices (cmesh, 4, vertices + 36, 3);
  t8_cmesh_set_tree_vertices (cmesh, 5, vertices + 45, 3);
  t8_cmesh_set_tree_vertices (cmesh, 6, vertices + 54, 3);
  t8_cmesh_set_tree_vertices (cmesh, 7, vertices + 63, 3);

  /* 6. Definition of the face neighbors between the different trees */
  t8_cmesh_set_join (cmesh, 0, 1, 1, 1, 0);
  t8_cmesh_set_join (cmesh, 0, 3, 0, 2, 0);
  t8_cmesh_set_join (cmesh, 1, 4, 0, 2, 0);
  t8_cmesh_set_join (cmesh, 2, 3, 1, 1, 0);
  t8_cmesh_set_join (cmesh, 3, 6, 0, 2, 0);
  t8_cmesh_set_join (cmesh, 4, 5, 1, 1, 0);
  t8_cmesh_set_join (cmesh, 4, 7, 0, 2, 0);
  t8_cmesh_set_join (cmesh, 6, 7, 1, 1, 0);



  /* 7. Commit the mesh */
  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}

/* This cmesh is [0,1]² with 2 triangles. for debugging purposes to compare to results from Florian Sieglar */
t8_cmesh_t
t8_cmesh_new_debugging (sc_MPI_Comm comm)
{

  /* 1. Defining an array with all vertices */
  /* Just all vertices of all trees. partly duplicated */
  double vertices[18] = { 0,0,0,0,1,0,1,0,0, //triangle 1
                          1,1,0,1,0,0,0,1,0, //triangle 2
  };

  /* 2. Initialization of the mesh */
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);

  /* 3. Definition of the geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);

  /* 4. Definition of the classes of the different trees */
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);


  /* 5. Classification of the vertices for each tree */
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 3);
  t8_cmesh_set_tree_vertices (cmesh, 1, vertices + 9, 3);

  /* 6. Definition of the face neighbors between the different trees */
  t8_cmesh_set_join (cmesh, 0, 1, 0, 0, 0);

  /* 7. Commit the mesh */
  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}


/* This cmesh is an octagon with 14 elements */
t8_cmesh_t
t8_cmesh_new_octagon (sc_MPI_Comm comm)
{
  double a=2.4142135623731; //1+sqrt(2)
  /* 1. Defining an array with all vertices */
  /* Just all vertices of all trees. partly duplicated */
  double vertices[126] = { -1,-a,0,-1,-1,0,-a,-1,0, //triangle 1
                          -1,-a,0,1,-a,0,1,-1,0, //triangle 2
                          -1,-a,0,-1,-1,0,1,-1,0,   //triangle 3
                          1,-a,0,1,-1,0,a,-1,0, //triangle 4
                          -a,-1,0,-1,-1,0,-1,1,0, //triangle 5
                          -a,-1,0,-a,1,0,-1,1,0,     //triangle 6
                          -1,-1,0,1,-1,0,1,1,0,   //triangle 7
                          -1,-1,0,-1,1,0,1,1,0,   //triangle 8
                          1,-1,0,a,-1,0,1,1,0, //triangle 9
                          a,-1,0,1,1,0,a,1,0, //triangle 10
                          -a,1,0,-1,1,0,-1,a,0,     //triangle 11
                          -1,1,0,1,a,0,1,1,0,   //triangle 12
                          -1,1,0,-1,a,0,1,a,0,   //triangle 13
                          1,a,0,1,1,0,a,1,0,   //triangle 14
  };

  /* 2. Initialization of the mesh */
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);

  /* 3. Definition of the geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);

  /* 4. Definition of the classes of the different trees */
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 3, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 4, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 5, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 6, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 7, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 8, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 9, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 10, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 11, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 12, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 13, T8_ECLASS_TRIANGLE);

  /* 5. Classification of the vertices for each tree */
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 3);
  t8_cmesh_set_tree_vertices (cmesh, 1, vertices + 9, 3);
  t8_cmesh_set_tree_vertices (cmesh, 2, vertices + 18, 3);
  t8_cmesh_set_tree_vertices (cmesh, 3, vertices + 27, 3);
  t8_cmesh_set_tree_vertices (cmesh, 4, vertices + 36, 3);
  t8_cmesh_set_tree_vertices (cmesh, 5, vertices + 45, 3);
  t8_cmesh_set_tree_vertices (cmesh, 6, vertices + 54, 3);
  t8_cmesh_set_tree_vertices (cmesh, 7, vertices + 63, 3);
  t8_cmesh_set_tree_vertices (cmesh, 8, vertices + 72, 3);
  t8_cmesh_set_tree_vertices (cmesh, 9, vertices + 81, 3);
  t8_cmesh_set_tree_vertices (cmesh, 10, vertices + 90, 3);
  t8_cmesh_set_tree_vertices (cmesh, 11, vertices + 99, 3);
  t8_cmesh_set_tree_vertices (cmesh, 12, vertices + 108, 3);
  t8_cmesh_set_tree_vertices (cmesh, 13, vertices + 117, 3);

  /* 6. Definition of the face neighbors between the different trees */
  t8_cmesh_set_join (cmesh, 0, 2, 2, 2, 0);
  t8_cmesh_set_join (cmesh, 0, 4, 0, 2, 1);
  t8_cmesh_set_join (cmesh, 1, 2, 1, 1, 0);
  t8_cmesh_set_join (cmesh, 1, 3, 0, 2, 0);
  t8_cmesh_set_join (cmesh, 2, 6, 0, 2, 0);
  t8_cmesh_set_join (cmesh, 3, 8, 0, 2, 0);
  t8_cmesh_set_join (cmesh, 4, 5, 1, 1, 0);
  t8_cmesh_set_join (cmesh, 4, 7, 0, 2, 0);
  t8_cmesh_set_join (cmesh, 5, 10, 0, 2, 0);
  t8_cmesh_set_join (cmesh, 6, 7, 1, 1, 0);
  t8_cmesh_set_join (cmesh, 6, 8, 0, 1, 1);
  t8_cmesh_set_join (cmesh, 7, 11, 0, 1, 1);
  t8_cmesh_set_join (cmesh, 8, 9, 0, 2, 0);
  t8_cmesh_set_join (cmesh, 9, 13, 0, 0, 0);
  t8_cmesh_set_join (cmesh, 10, 12, 0, 2, 0);
  t8_cmesh_set_join (cmesh, 11, 12, 2, 1, 1);
  t8_cmesh_set_join (cmesh, 11, 13, 0, 2, 0);

  /* 7. Commit the mesh */
  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}


/* This cmesh is a complex polygonal shape with 8 triangles. */
t8_cmesh_t
t8_cmesh_new_complex_polygonal_shape (sc_MPI_Comm comm)
{

  /* 1. Defining an array with all vertices */
  /* Just all vertices of all trees. partly duplicated */
  double vertices[72] = { 0.4,0.5,0,1,0,0,1,0.6,0, //triangle 1
                          0.4,0.2,0,0.4,0.5,0,1,0,0, //triangle 2
                          0,0,0,0.4,0.2,0,0.4,0.5,0,   //triangle 3
                          0,0,0,0.2,0.4,0,0.4,0.5,0, //triangle 4
                          0.4,0.8,0,0.2,0.4,0,0.4,0.5,0, //triangle 5
                          0.4,0.8,0,0.2,0.4,0,0,0.6,0,     //triangle 6
                          0.4,0.8,0,0.4,1,0,0,0.6,0,   //triangle 7
                          0.4,0.8,0,0.4,1,0,1,0.8,0,   //triangle 8
  };

  /* 2. Initialization of the mesh */
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);

  /* 3. Definition of the geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);

  /* 4. Definition of the classes of the different trees */
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 3, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 4, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 5, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 6, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 7, T8_ECLASS_TRIANGLE);


  /* 5. Classification of the vertices for each tree */
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 3);
  t8_cmesh_set_tree_vertices (cmesh, 1, vertices + 9, 3);
  t8_cmesh_set_tree_vertices (cmesh, 2, vertices + 18, 3);
  t8_cmesh_set_tree_vertices (cmesh, 3, vertices + 27, 3);
  t8_cmesh_set_tree_vertices (cmesh, 4, vertices + 36, 3);
  t8_cmesh_set_tree_vertices (cmesh, 5, vertices + 45, 3);
  t8_cmesh_set_tree_vertices (cmesh, 6, vertices + 54, 3);
  t8_cmesh_set_tree_vertices (cmesh, 7, vertices + 63, 3);

  /* 6. Definition of the face neighbors between the different trees */
  t8_cmesh_set_join (cmesh, 0, 1, 2, 0, 0);
  t8_cmesh_set_join (cmesh, 1, 2, 2, 0, 0);
  t8_cmesh_set_join (cmesh, 2, 3, 1, 1, 0);
  t8_cmesh_set_join (cmesh, 3, 4, 0, 0, 0);
  t8_cmesh_set_join (cmesh, 4, 5, 2, 2, 0);
  t8_cmesh_set_join (cmesh, 5, 6, 1, 1, 0);
  t8_cmesh_set_join (cmesh, 6, 7, 2, 2, 0);



  /* 7. Commit the mesh */
  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}


/* This cmesh is an L-shape with 4 triangles. */
t8_cmesh_t
t8_cmesh_new_l_shape (sc_MPI_Comm comm)
{

  /* 1. Defining an array with all vertices */
  /* Just all vertices of all trees. partly duplicated */
  double vertices[36] = { 0.5,0.5,0,1,0,0,1,0.5,0, //triangle 1
                          0.5,0.5,0,1,0,0,0,0,0, //triangle 2
                          0.5,0.5,0,0,1,0,0,0,0,   //triangle 3
                          0.5,0.5,0,0,1,0,0.5,1,0, //triangle 4
  };

  /* 2. Initialization of the mesh */
  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);

  /* 3. Definition of the geometry */
  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);

  /* 4. Definition of the classes of the different trees */
  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 3, T8_ECLASS_TRIANGLE);


  /* 5. Classification of the vertices for each tree */
  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 3);
  t8_cmesh_set_tree_vertices (cmesh, 1, vertices + 9, 3);
  t8_cmesh_set_tree_vertices (cmesh, 2, vertices + 18, 3);
  t8_cmesh_set_tree_vertices (cmesh, 3, vertices + 27, 3);

  /* 6. Definition of the face neighbors between the different trees */
  t8_cmesh_set_join (cmesh, 0, 1, 2, 2, 0);
  t8_cmesh_set_join (cmesh, 1, 2, 1, 1, 0);
  t8_cmesh_set_join (cmesh, 2, 3, 2, 2, 0);



  /* 7. Commit the mesh */
  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}

/* This cmesh is (x-axis)[0,360]x(y-axis)[-90,90] with 4 triangles (specifically for the era5 data). */
t8_cmesh_t
t8_cmesh_new_earth (sc_MPI_Comm comm)
{


  double vertices[36] = { 0,-90,0,180,-90,0,180,90,0,//triangle 1
                          360,90,0,180,-90,0,360,-90,0, //triangle 2
                          0,-90,0,0,90,0,180,90,0, //triangle 3
                          360,90,0,180,-90,0,180,90,0   //triangle 4
  };

  t8_cmesh_t cmesh;
  t8_cmesh_init (&cmesh);


  t8_cmesh_register_geometry<t8_geometry_linear> (cmesh);


  t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 2, T8_ECLASS_TRIANGLE);
  t8_cmesh_set_tree_class (cmesh, 3, T8_ECLASS_TRIANGLE);



  t8_cmesh_set_tree_vertices (cmesh, 0, vertices, 3);
  t8_cmesh_set_tree_vertices (cmesh, 1, vertices + 9, 3);
  t8_cmesh_set_tree_vertices (cmesh, 2, vertices + 18, 3);
  t8_cmesh_set_tree_vertices (cmesh, 3, vertices + 27, 3);


  t8_cmesh_set_join (cmesh, 0, 2, 1, 1, 0);
  t8_cmesh_set_join (cmesh, 0, 3, 0, 0, 0);
  t8_cmesh_set_join (cmesh, 1, 3, 2, 2, 0);

  t8_cmesh_commit (cmesh, comm);

  return cmesh;
}



/* These are six test functions. */
double F(double x, double y) {
  if ((x == -1.) && (y == -1.)) return 1.;
  return sin(2.*M_PI*x)*sin(2.*M_PI*y);
}

// 5. Main program entry point
int main(int argc, char** argv) {
  int mpiret;
  sc_MPI_Comm comm;
  /* Initialize MPI. This has to happen before we initialize sc or t8code. */
  mpiret = sc_MPI_Init (&argc, &argv);
  /* Error check the MPI return value. */
  SC_CHECK_MPI (mpiret);
  /* Initialize the sc library, has to happen before we initialize t8code. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
  /* Initialize t8code with log level SC_LP_PRODUCTION. See sc.h for more info on the log levels. */
  t8_init (SC_LP_PRODUCTION);
  /* We will use MPI_COMM_WORLD as a communicator. */
  comm = sc_MPI_COMM_WORLD;
  /* Creation of a basic two dimensional cmesh. */
  t8_cmesh_t cmesh = t8_cmesh_new_debugging (comm);
  const t8_scheme *scheme = t8_scheme_new_default ();
  t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, 0, 0, comm);
  const char *prefix_forest = "t8_mra_forest";
  t8_forest_write_vtk (forest, prefix_forest);
  levelgrid_map<t8_data_per_element_1d_gh> grid(4);
  //levelgrid_map<t8_data_per_element_1d_gh>* grid = new levelgrid_map<t8_data_per_element_1d_gh>(4);
  struct t8_data_per_element_1d * element_data;
  element_data=t8_create_element_data (grid, forest,F, 10, 4);
  // 5.5. Return a value indicating successful execution
  /* Finalize the sc library */
  /* Destroy the forest. */
  t8_forest_unref (&forest);
  // Do something with the grid
  //delete grid;
  t8_global_productionf (" [step5] Destroyed forest.\n");
  T8_FREE (element_data);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return 0;  // Return 0 to indicate successful execution
}
