#include <t8.h>                                     /* General t8code header, always include this. */
#include <t8_cmesh.hxx>                             /* cmesh definition and basic interface. */
#include <t8_forest/t8_forest_general.h>            /* forest definition and basic interface. */
#include <t8_vtk/t8_vtk_writer.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.hxx> /* linear geometry of the cmesh */
#include <t8_forest/t8_forest_io.h>                 /* forest io interface. */
#include <t8_schemes/t8_default/t8_default.hxx> /* default refinement scheme. */
#include <t8_forest/t8_forest_geometrical.h>        /* geometrical information */
#include "vecmat.hxx"
#include "basis_functions.hxx"
#include "mask_coefficients.hxx"
#include "dunavant.hxx"
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
#include "t8_unordered_dense.hxx"
#include "t8_gsl.hxx"
// #include <gsl/gsl_math.h>
// #include <gsl/gsl_interp2d.h>
// #include <gsl/gsl_spline2d.h>
#include <fstream>
#include <time.h>

typedef double (*func)(double, double);
/* declaring the matrices */
mat M0, M1, M2, M3, N0, N1, N2, N3;
/* defines the maximum refinement level (should currently be less than approx. 10 depending on how many basecells there are) */
const int max_level=4;
typedef int8_t t8_dtri_cube_id_t;
// Precompute powers of 4 up to level 29(T8_DTRI_MAXLEVEL)
#define MAX_LEVEL T8_DTRI_MAXLEVEL
//BIT Length for bitwise level multi index
const int PATH_BITS = 2;   // Each path segment is 2 bits
const int LEVEL_BITS = 5;  // Level is encoded in 5 bits
const int BASECELL_BITS = 21; // Basecell is encoded in 21 bits
uint64_t pow4[MAX_LEVEL + 1];

// Function to initialize the powers of 4
void initialize_pow4() {
    pow4[0] = 1;
    for (int i = 1; i <= MAX_LEVEL; ++i) {
        pow4[i] = pow4[i - 1] * 4;
    }
}

/* To store the data for the 3d case (e.g. in AuswertungSinglescale) */
struct double_3d_array{
  double dim_val[3];
};


/* We can drop the volume of the element and the level, hasFather, haschilds, childs ids, father ids*/
/* This struct stores the element data */
struct t8_data_per_element_1d
  {
    double u_coeff[M_mra]; //single-scale coefficients for all dof/ basis polynomials
    uint64_t lmi;
  };

/* This struct stores the element data */
struct t8_data_per_element_3d
  {
    double u_coeff_d1[M_mra]; //single-scale coefficients for all dof/ basis polynomials
    double u_coeff_d2[M_mra];
    double u_coeff_d3[M_mra];
    uint64_t lmi;
  };

/* This struct stores the element data */
struct adapt_data_1d_wavelet_based
  {
    double gamma;
    double C_thr;
    int current_level;
    levelgrid_map<t8_data_per_element_1d_gh>* grid_map_ptr;
    func my_func;
  };

/* This struct stores the element data */
struct adapt_data_3d_wavelet_based
  {
    double gamma;
    double C_thr;
    int current_level;
    levelgrid_map<t8_data_per_element_3d_gh>* grid_map_ptr;
    func my_func1;
    func my_func2;
    func my_func3;
  };

/* This struct stores the element data */
struct adapt_data_1d_wavelet_free
  {
    double gamma;
    double C_thr;
    int current_level;
    levelgrid_map<t8_data_per_element_waveletfree_1d_gh>* grid_map_ptr;
    func my_func;
  };

/* This struct stores the element data */
struct adapt_data_3d_wavelet_free
  {
    double gamma;
    double C_thr;
    int current_level;
    levelgrid_map<t8_data_per_element_waveletfree_3d_gh>* grid_map_ptr;
    func my_func1;
    func my_func2;
    func my_func3;
  };

/* This struct stores the element data */
struct adapt_data_1d_wavelet_based_prediction
  {
    double gamma;
    double C_thr;
    int current_level;
    levelgrid_map<t8_data_per_element_1d_gh>* grid_map_ptr_before;
    levelgrid_map<t8_data_per_element_1d_gh>* grid_map_ptr_after;
    func my_func;
  };

/* This struct stores the element data */
struct adapt_data_3d_wavelet_based_prediction
  {
    double gamma;
    double C_thr;
    int current_level;
    levelgrid_map<t8_data_per_element_3d_gh>* grid_map_ptr_before;
    levelgrid_map<t8_data_per_element_3d_gh>* grid_map_ptr_after;
    func my_func1;
    func my_func2;
    func my_func3;
  };

/* This struct stores the element data */
struct adapt_data_1d_wavelet_free_prediction
  {
    double gamma;
    double C_thr;
    int current_level;
    levelgrid_map<t8_data_per_element_waveletfree_1d_gh>* grid_map_ptr_before;
    levelgrid_map<t8_data_per_element_waveletfree_1d_gh>* grid_map_ptr_after;
    func my_func;
  };

/* This struct stores the element data */
struct adapt_data_3d_wavelet_free_prediction
  {
    double gamma;
    double C_thr;
    int current_level;
    levelgrid_map<t8_data_per_element_waveletfree_3d_gh>* grid_map_ptr_before;
    levelgrid_map<t8_data_per_element_waveletfree_3d_gh>* grid_map_ptr_after;
    func my_func1;
    func my_func2;
    func my_func3;
  };

  /* Compute the cube-id of t's ancestor of level "level" in constant time.
   * If "level" is greater then t->level then the cube-id 0 is returned. */
  static t8_dtri_cube_id_t
  compute_cubeid (const t8_dtri_t *t, int level)
  {
    t8_dtri_cube_id_t id = 0;
    t8_dtri_coord_t h;

    /* TODO: assert that 0 < level? This may simplify code elsewhere */

    T8_ASSERT (0 <= level && level <= T8_DTRI_MAXLEVEL);
    h = T8_DTRI_LEN (level);
    if (level == 0) {
      return 0;
    }

    id |= ((t->x & h) ? 0x01 : 0);
    id |= ((t->y & h) ? 0x02 : 0);
  #ifdef T8_DTRI_TO_DTET
    id |= ((t->z & h) ? 0x04 : 0);
  #endif

    return id;
  }

  /* A routine to compute the type of t's ancestor of level "level", if its type at an intermediate level is already
   * known. If "level" equals t's level then t's type is returned. It is not allowed to call this function with "level"
   * greater than t->level. This method runs in O(t->level - level).
   */
  static t8_dtri_type_t
  compute_type_ext (const t8_dtri_t *t, int level, t8_dtri_type_t known_type, int known_level)
  {
    int8_t type = known_type;
    t8_dtri_cube_id_t cid;
    int i;

    T8_ASSERT (0 <= level && level <= known_level);
    T8_ASSERT (known_level <= t->level);
    if (level == known_level) {
      return known_type;
    }
    if (level == 0) {
      /* TODO: the type of the root tet is hardcoded to 0
       *       maybe once we want to allow the root tet to have different types */
      return 0;
    }
    for (i = known_level; i > level; i--) {
      cid = compute_cubeid (t, i);
      /* compute type as the type of T^{i+1}, that is T's ancestor of level i+1 */
      type = t8_dtri_cid_type_to_parenttype[cid][type];
    }
    return type;
  }

  /* A routine to compute the type of t's ancestor of level "level". If "level" equals t's level then t's type is
   * returned. It is not allowed to call this function with "level" greater than t->level. This method runs in
   * O(t->level - level).
   */
  static t8_dtri_type_t
  compute_type (const t8_dtri_t *t, int level)
  {
    return compute_type_ext (t, level, t->type, t->level);
  }

  static void get_point_order(int *first,int *second, int *third, t8_dtri_cube_id_t cube_id){
      if(*first==0 && *second==1 && *third==2){
        if(cube_id==0){
          *first=0;
          *second=1;
          *third=2;
        }
        else if(cube_id==1){
          *first=2;//1
          *second=0;//2
          *third=1;//0
        }
        else if(cube_id==2){
          *first=1;//2
          *second=2;//0
          *third=0;//1
        }
        else if(cube_id==3){
          *first=0;
          *second=2;
          *third=1;
        }
      }
      else if(*first==2 && *second==0 && *third==1){
        if(cube_id==0){
          *first=0;
          *second=1;
          *third=2;
        }
        else if(cube_id==1){
          *first=2;
          *second=0;
          *third=1;
        }
        else if(cube_id==2){
          *first=1;
          *second=2;
          *third=0;
        }
        else if(cube_id==3){
          *first=2;//1
          *second=1;//0
          *third=0;//2
        }
      }
      else if(*first==1 && *second==2 && *third==0){
        if(cube_id==0){
          *first=0;
          *second=1;
          *third=2;
        }
        else if(cube_id==1){
          *first=2;
          *second=0;
          *third=1;
        }
        else if(cube_id==2){
          *first=1;
          *second=2;
          *third=0;
        }
        else if(cube_id==3){
          *first=1;//2
          *second=0;//1
          *third=2;//0
        }
      }
      else if(*first==0 && *second==2 && *third==1){
        if(cube_id==0){
          *first=0;
          *second=2;
          *third=1;
        }
        else if(cube_id==1){
          *first=1;
          *second=0;
          *third=2;
        }
        else if(cube_id==2){
          *first=2;
          *second=1;
          *third=0;
        }
        else if(cube_id==3){
          *first=2;
          *second=0;
          *third=1;
        }
      }
      else if(*first==1 && *second==0 && *third==2){
        if(cube_id==0){
          *first=0;
          *second=2;
          *third=1;
        }
        else if(cube_id==1){
          *first=1;
          *second=0;
          *third=2;
        }
        else if(cube_id==2){
          *first=2;
          *second=1;
          *third=0;
        }
        else if(cube_id==3){
          *first=0;
          *second=1;
          *third=2;
        }
      }
      else if(*first==2 && *second==1 && *third==0){
        if(cube_id==0){
          *first=0;
          *second=2;
          *third=1;
        }
        else if(cube_id==1){
          *first=1;
          *second=0;
          *third=2;
        }
        else if(cube_id==2){
          *first=2;
          *second=1;
          *third=0;
        }
        else if(cube_id==3){
          *first=1;
          *second=2;
          *third=0;
        }
      }
  }


  // static void get_point_order(int *first, int *second, int *third, t8_dtri_cube_id_t cube_id) {
  //     // Create a 2D array (6x4) for the lookup table
  //     // The table stores the results for each permutation of {*first, *second, *third} and cube_id
  //     // Permutation indices for {*first, *second, *third}: (0,1,2), (0,2,1), (1,0,2), (1,2,0), (2,0,1), (2,1,0)
  //
  //     // Table format: [permutation_index][cube_id] -> [first, second, third]
  //     int lookup[6][4][3] = {
  //         // For permutation (0,1,2)
  //         {{0, 1, 2}, {2, 0, 1}, {1, 2, 0}, {0, 2, 1}},
  //         // For permutation (0,2,1)
  //         {{0, 1, 2}, {2, 0, 1}, {1, 2, 0}, {2, 1, 0}},
  //         // For permutation (1,2,0)
  //         {{0, 1, 2}, {2, 0, 1}, {1, 2, 0}, {1, 0, 2}},
  //         // For permutation (0,2,1)
  //         {{0, 2, 1}, {1, 0, 2}, {2, 1, 0}, {2, 0, 1}},
  //         // For permutation (2,0,1)
  //         {{0, 2, 1}, {1, 0, 2}, {2, 1, 0}, {0, 1, 2}},
  //         // For permutation (2,1,0)
  //         {{0, 1, 2}, {1, 0, 2}, {2, 1, 0}, {0, 2, 1}},
  //     };
  //
  //     // Map the values of *first, *second, *third to an index
  //     int permutation_index = -1;
  //     if (*first == 0 && *second == 1 && *third == 2) {
  //         permutation_index = 0;
  //     } else if (*first == 0 && *second == 2 && *third == 1) {
  //         permutation_index = 1;
  //     } else if (*first == 1 && *second == 2 && *third == 0) {
  //         permutation_index = 2;
  //     } else if (*first == 0 && *second == 1 && *third == 2) {
  //         permutation_index = 3;
  //     } else if (*first == 1 && *second == 0 && *third == 2) {
  //         permutation_index = 4;
  //     } else if (*first == 2 && *second == 1 && *third == 0) {
  //         permutation_index = 5;
  //     }
  //
  //     if (permutation_index != -1) {
  //         // Use the lookup table to directly set *first, *second, *third for the correct cube_id
  //         *first = lookup[permutation_index][cube_id][0];
  //         *second = lookup[permutation_index][cube_id][1];
  //         *third = lookup[permutation_index][cube_id][2];
  //     }
  // }

  static void invert_order(int *first, int *second, int *third){
    int first_new=*first;
    int second_new=*second;
    int third_new=*third;
    if(first_new==0){
      *first=0;
    }
    else if(second_new==0){
      *first=1;
    }
    else if(third_new==0){
      *first=2;
    }
    if(first_new==1){
      *second=0;
    }
    else if(second_new==1){
      *second=1;
    }
    else if(third_new==1){
      *second=2;
    }
    if(first_new==2){
      *third=0;
    }
    else if(second_new==2){
      *third=1;
    }
    else if(third_new==2){
      *third=2;
    }
  }


  // static void invert_order(int *first, int *second, int *third) {
  //     int values[3] = {*first, *second, *third};
  //
  //     // Direct mapping of the values
  //     *first = (values[0] == 0) ? 0 : (values[1] == 0) ? 1 : 2;
  //     *second = (values[0] == 1) ? 0 : (values[1] == 1) ? 1 : 2;
  //     *third = (values[0] == 2) ? 0 : (values[1] == 2) ? 1 : 2;
  // }


  int get_correct_order_children(int type,int child_id,int first, int second, int third){
    if(type==1){
      if(first==0&&second==1&&third==2){
        if(child_id==0){
          return 1;//3;
        }
        else if(child_id==1){
          return 0;//0;
        }
        else if(child_id==2){
          return 2;//1;
        }
        else if(child_id==3){
          return 3;//2;
        }
      }
      else if(first==2&&second==0&&third==1){
        if(child_id==0){
          return 1;//2;
        }
        else if(child_id==1){
          return 3;//3;
        }
        else if(child_id==2){
          return 0;//1;
        }
        else if(child_id==3){
          return 2;//0;
        }
      }
      else if(first==1&&second==2&&third==0){
        if(child_id==0){
          return 1;//3;
        }
        else if(child_id==1){
          return 2;//1;
        }
        else if(child_id==2){
          return 3;//2;
        }
        else if(child_id==3){
          return 0;//0;
        }
      }
      else if(first==0&&second==2&&third==1){
        if(child_id==0){
          return 1;//1;
        }
        else if(child_id==1){
          return 0;//3;
        }
        else if(child_id==2){
          return 3;//2;
        }
        else if(child_id==3){
          return 2;//0;
        }
      }
      else if(first==1&&second==0&&third==2){
        if(child_id==0){
          return 1;//2;
        }
        else if(child_id==1){
          return 2;//1;
        }
        else if(child_id==2){
          return 0;//3;
        }
        else if(child_id==3){
          return 3;//0;
        }
      }
      else if(first==2&&second==1&&third==0){
        if(child_id==0){
          return 1;//3;
        }
        else if(child_id==1){
          return 3;//2;
        }
        else if(child_id==2){
          return 2;//2;
        }
        else if(child_id==3){
          return 0;//0;
        }
      }
    }
    else{
      if(first==0&&second==1&&third==2){
        if(child_id==0){
          return 2;//3;
        }
        else if(child_id==1){
          return 0;//0;
        }
        else if(child_id==2){
          return 1;//1;
        }
        else if(child_id==3){
          return 3;//2;
        }
      }
      else if(first==2&&second==0&&third==1){
        if(child_id==0){
          return 2;//2;
        }
        else if(child_id==1){
          return 3;//3;
        }
        else if(child_id==2){
          return 0;//1;
        }
        else if(child_id==3){
          return 1;//0;
        }
      }
      else if(first==1&&second==2&&third==0){
        if(child_id==0){
          return 2;//3;
        }
        else if(child_id==1){
          return 1;//1;
        }
        else if(child_id==2){
          return 3;//2;
        }
        else if(child_id==3){
          return 0;//0;
        }
      }
      else if(first==0&&second==2&&third==1){
        if(child_id==0){
          return 2;//1;
        }
        else if(child_id==1){
          return 0;//3;
        }
        else if(child_id==2){
          return 3;//2;
        }
        else if(child_id==3){
          return 1;//0;
        }
      }
      else if(first==1&&second==0&&third==2){
        if(child_id==0){
          return 2;//2;
        }
        else if(child_id==1){
          return 1;//1;
        }
        else if(child_id==2){
          return 0;//3;
        }
        else if(child_id==3){
          return 3;//0;
        }
      }
      else if(first==2&&second==1&&third==0){
        if(child_id==0){
          return 2;//3;
        }
        else if(child_id==1){
          return 3;//2;
        }
        else if(child_id==2){
          return 1;//2;
        }
        else if(child_id==3){
          return 0;//0;
        }
      }
    }
  }


//   int get_correct_order_children(int type, int child_id, int first, int second, int third) {
//     // Lookup tables for the two types (type == 1 or type == 2)
//     static const int lookup_type_1[6][4] = {
//         {1, 0, 2, 3}, // first=0, second=1, third=2
//         {1, 3, 0, 2}, // first=2, second=0, third=1
//         {1, 2, 3, 0}, // first=1, second=2, third=0
//         {1, 0, 3, 2}, // first=0, second=2, third=1
//         {1, 2, 0, 3}, // first=1, second=0, third=2
//         {1, 3, 2, 0}  // first=2, second=1, third=0
//     };
//
//     static const int lookup_type_2[6][4] = {
//         {2, 0, 1, 3}, // first=0, second=1, third=2
//         {2, 3, 0, 1}, // first=2, second=0, third=1
//         {2, 1, 3, 0}, // first=1, second=2, third=0
//         {2, 0, 3, 1}, // first=0, second=2, third=1
//         {2, 1, 0, 3}, // first=1, second=0, third=2
//         {2, 3, 1, 0}  // first=2, second=1, third=0
//     };
//
//     // Determine the correct lookup table based on the 'type'
//     const int (*lookup)[4] = (type == 1) ? lookup_type_1 : lookup_type_2;
//
//     // Map first, second, and third to an index (0-5)
//     int index = (first == 0 && second == 1 && third == 2) ? 0 :
//                 (first == 2 && second == 0 && third == 1) ? 1 :
//                 (first == 1 && second == 2 && third == 0) ? 2 :
//                 (first == 0 && second == 2 && third == 1) ? 3 :
//                 (first == 1 && second == 0 && third == 2) ? 4 : 5;
//
//     // Return the result from the lookup table
//     return lookup[index][child_id];
// }

// Function to initialize an LMI with the given basecell, path set to 0, and level set to 0
uint64_t initialize_lmi(uint64_t basecell) {
    // Path will be initialized to 0
    uint64_t path = 0;

    // Level will be set to 0
    uint64_t level = 0;

    // Combine basecell, level, and path into the LMI
    uint64_t lmi = (path << (LEVEL_BITS + BASECELL_BITS)) | (level << BASECELL_BITS) | basecell;

    return lmi;
}

uint64_t calculate_lmi(uint64_t basecell,const t8_element_t *elem,t8_eclass_t tree_class){
  uint64_t lmi=initialize_lmi(basecell);
  int first=0;
  int second=1;
  int third=2;
  for(int level=0;level<element_get_level(elem);level++){
    int ancestor_id=element_get_ancestor_id (tree_class, elem,level+1);
    t8_dtri_type_t type=compute_type (((t8_dtri_t *)elem),level);
    int correct_child=get_correct_order_children((int)type,ancestor_id,first,second,third);
    lmi=get_jth_child_lmi_binary(lmi, correct_child);
    get_point_order(&first,&second, &third, t8_dtri_type_cid_to_beyid[type][ancestor_id]);
  }
  return lmi;
}



  inline bool isZero(double x)
  {
    const double epsilon = 1e-15;
    return std::abs(x) <= epsilon;
  }

  // Function to decrease the level and reset corresponding path bits to zero, returning the new LMI
  uint64_t get_parents_lmi_binary(uint64_t lmi) {
      // Extract the basecell (21 bits) - the lowest bits
      uint64_t basecell = lmi & ((1ULL << BASECELL_BITS) - 1);
      lmi >>= BASECELL_BITS;  // Shift to remove the basecell part

      // Extract the level (5 bits)
      uint64_t level = lmi & ((1ULL << LEVEL_BITS) - 1);
      lmi >>= LEVEL_BITS;  // Shift to remove the level part

      // Extract the path (38 bits)
      uint64_t path = lmi;  // Remaining part is the path

      // Decrease the level by 1 if it is greater than 0
      if (level > 0) {
          level--;  // Decrease the level

          // Reset the path bits for the decreased level (for this we just shift out the last 2 bits for the current level)
          path >>= PATH_BITS;  // Shift to remove the last path bits

          // Re-encode the LMI with the decreased level, reset path, and basecell
          // Reinsert the level (shifted to its position)
          lmi = (path << (LEVEL_BITS + BASECELL_BITS)) | (level << BASECELL_BITS) | basecell;  // Reassemble the LMI
      }

      return lmi;
  }

  // Function to increase the level by 1 and insert a new path segment (j) at the corresponding level, returning the new LMI
  uint64_t get_jth_child_lmi_binary(uint64_t lmi, uint64_t j) {
      // Extract the basecell (21 bits) from the lowest bits
      uint64_t basecell = lmi & ((1ULL << BASECELL_BITS) - 1);
      lmi >>= BASECELL_BITS;  // Remove the basecell part

      // Extract the level (5 bits)
      uint64_t level = lmi & ((1ULL << LEVEL_BITS) - 1);
      lmi >>= LEVEL_BITS;  // Remove the level part

      // Extract the path (38 bits) as the remaining part
      uint64_t path = lmi;  // Remaining part after removing level and basecell

      // Increase the level by 1
      level++;  // Increment the level

      // Insert the new path segment (j) at the corresponding position (2 bits)
      path = (path << 2) | (j & 3);  // Add the new path segment by shifting and OR-ing with j (ensuring j is between 0 and 3)

      // Re-encode the LMI with the increased level, new path segment, and basecell in the correct order:
      // 1. Shift the path into its position
      // 2. Insert the level (shifted into the correct position)
      // 3. Insert the basecell (shifted into the correct position)
      lmi = (path << (LEVEL_BITS + BASECELL_BITS)) | (level << BASECELL_BITS) | basecell;

      return lmi;
  }

  // Function to initialize the struct
  void initialize_t8_data_per_element_1d(struct t8_data_per_element_1d_gh *data) {
      // Initialize u_coeff and d_coeff to zero
      for (int i = 0; i < M_mra; i++) {
          data->u_coeff[i] = 0.0;  // Set all u_coeff to 0
      }
      for (int i = 0; i < 3 * M_mra; i++) {
          data->d_coeff[i] = 0.0;  // Set all d_coeff to 0
      }

      // Set significant to false
      data->significant = false;

      // Set the bit fields
      data->first = 0;
      data->second = 1;
      data->third = 2;
  }

  // Function to initialize t8_data_per_element_3d_gh struct
  void initialize_t8_data_per_element_3d(struct t8_data_per_element_3d_gh *data) {
      // Initialize u_coeff_d1, u_coeff_d2, u_coeff_d3 to zero
      for (int i = 0; i < M_mra; i++) {
          data->u_coeff_d1[i] = 0.0;
          data->u_coeff_d2[i] = 0.0;
          data->u_coeff_d3[i] = 0.0;
      }

      // Initialize d_coeff_d1, d_coeff_d2, d_coeff_d3 to zero
      for (int i = 0; i < 3 * M_mra; i++) {
          data->d_coeff_d1[i] = 0.0;
          data->d_coeff_d2[i] = 0.0;
          data->d_coeff_d3[i] = 0.0;
      }

      // Set significant to false
      data->significant = false;

      // Set the bit fields
      data->first = 0;
      data->second = 1;
      data->third = 2;
  }


  // Function to initialize t8_data_per_element_waveletfree_1d_gh struct
  void initialize_t8_data_per_element_waveletfree_1d(struct t8_data_per_element_waveletfree_1d_gh *data) {
      // Initialize u_coeff to zero
      data->u_coeff={0};
      // for (int i = 0; i < M_mra; i++) {
      //     data->u_coeff[i] = 0.0;
      // }
      data->d_coeff_wavelet_free = {0};
      // // Initialize d_coeff_wavelet_free to zero
      // for (int i = 0; i < M_mra; i++) {
      //     for (int j = 0; j < 4; j++) {
      //         data->d_coeff_wavelet_free[i][j] = 0.0;
      //     }
      // }

      // Set significant to false
      data->significant = false;

      // Set the bit fields
      data->first = 0;
      data->second = 1;
      data->third = 2;
  }

  // Function to initialize t8_data_per_element_waveletfree_3d_gh struct
  void initialize_t8_data_per_element_waveletfree_3d(struct t8_data_per_element_waveletfree_3d_gh *data) {
      // Initialize u_coeff_d1, u_coeff_d2, u_coeff_d3 to zero
      for (int i = 0; i < M_mra; i++) {
          data->u_coeff_d1[i] = 0.0;
          data->u_coeff_d2[i] = 0.0;
          data->u_coeff_d3[i] = 0.0;
      }

      // Initialize d_coeff_wavelet_free_d1, d_coeff_wavelet_free_d2, d_coeff_wavelet_free_d3 to zero
      for (int i = 0; i < M_mra; i++) {
          for (int j = 0; j < 4; j++) {
              data->d_coeff_wavelet_free_d1[i][j] = 0.0;
              data->d_coeff_wavelet_free_d2[i][j] = 0.0;
              data->d_coeff_wavelet_free_d3[i][j] = 0.0;
          }
      }

      // Set significant to false
      data->significant = false;

      // Set the bit fields
      data->first = 0;
      data->second = 1;
      data->third = 2;
  }


  /* The data that we want to store for each element.
   * In this example we want to store the element's level and volume. */
struct t8_data_per_element_1d *
  t8_create_element_data (levelgrid_map<t8_data_per_element_1d_gh>& grid_hierarchy, t8_forest_t forest,func F, const int rule, const int max_lev)
  {
    int order_num;
    double *wtab;
    double *xytab;
    double *xytab_ref;
    order_num = dunavant_order_num(rule);
    wtab = T8_ALLOC (double, order_num);
    xytab = T8_ALLOC (double, 2*order_num);
    xytab_ref = T8_ALLOC (double, 2*order_num);
    mat A;
    vector<int> r;
    dunavant_rule(rule, order_num, xytab_ref, wtab);
    t8_locidx_t num_local_elements;
    t8_locidx_t num_ghost_elements;
    struct t8_data_per_element_1d *element_data;
    //initial_grid_hierarchy.lev_arr[max_level].forest_arr

    /* Check that forest is a committed, that is valid and usable, forest. */
    T8_ASSERT (t8_forest_is_committed (forest));

    /* Get the number of local elements of forest. */
    num_local_elements = t8_forest_get_local_num_elements (forest);
    /* Get the number of ghost elements of forest. */
    num_ghost_elements = t8_forest_get_num_ghosts (forest);

    /* Now we need to build an array of our data that is as long as the number
     * of elements plus the number of ghosts. You can use any allocator such as
     * new, malloc or the t8code provide allocation macro T8_ALLOC.
     * Note that in the latter case you need
     * to use T8_FREE in order to free the memory.
     */
    element_data = T8_ALLOC (struct t8_data_per_element_1d, num_local_elements + num_ghost_elements);//hier
    /* Note: We will later need to associate this data with an sc_array in order to exchange the values for
     *       the ghost elements, which we can do with sc_array_new_data (see t8_step5_exchange_ghost_data).
     *       We could also have directly allocated the data here in an sc_array with
     *       sc_array_new_count (sizeof (struct data_per_element), num_local_elements + num_ghost_elements);
     */

    /* Let us now fill the data with something.
     * For this, we iterate through all trees and for each tree through all its elements, calling
     * t8_forest_get_element_in_tree to get a pointer to the current element.
     * This is the recommended and most performant way.
     * An alternative is to iterate over the number of local elements and use
     * t8_forest_get_element. However, this function needs to perform a binary search
     * for the element and the tree it is in, while t8_forest_get_element_in_tree has a
     * constant look up time. You should only use t8_forest_get_element if you do not know
     * in which tree an element is.
     */
    {
      t8_locidx_t itree, num_local_trees;
      t8_locidx_t current_index;
      t8_locidx_t ielement, num_elements_in_tree;
      t8_eclass_t tree_class;
      const t8_scheme *eclass_scheme;
      const t8_element_t *element;

      /* Get the number of trees that have elements of this process. */
      num_local_trees = t8_forest_get_num_local_trees (forest);
      // long long int basecell_num_digits_offset=countDigit(t8_forest_get_num_global_trees (initial_grid_hierarchy.lev_arr[level].forest_arr)-1)-1;
      for (itree = 0, current_index = 0; itree < num_local_trees; ++itree) {
        /* This loop iterates through all local trees in the forest. */
        /* Each tree may have a different element class (quad/tri/hex/tet etc.) and therefore
         * also a different way to interpret its elements. In order to be able to handle elements
         * of a tree, we need to get its eclass_scheme, and in order to so we first get its eclass. */
        tree_class = t8_forest_get_tree_class (forest, itree);
        eclass_scheme = t8_forest_get_scheme(forest);
        /* Get the number of elements of this tree. */
        num_elements_in_tree = t8_forest_get_tree_num_elements (forest, itree);
        for (ielement = 0; ielement < num_elements_in_tree; ++ielement, ++current_index) {
          /* This loop iterates through all the local elements of the forest in the current tree. */
          /* We can now write to the position current_index into our array in order to store
           * data for this element. */
          /* Since in this example we want to compute the data based on the element in question,
           * we need to get a pointer to this element. */

          element = t8_forest_get_element_in_tree (forest, itree, ielement);

          /* We want to store the elements level and its volume as data. We compute these
           * via the eclass_scheme and the forest_element interface. */
          //element_data[current_index].level = eclass_scheme->t8_element_level (element);
          double volume = t8_forest_element_volume (forest, itree, element);
          t8_gloidx_t base_element=t8_forest_global_tree_id (forest, itree);
          element_data[current_index].lmi=initialize_lmi((uint64_t)base_element);
          t8_global_productionf ("The lmi aus elem data is: %" PRIu64 "\n",element_data[current_index].lmi);
          t8_data_per_element_1d_gh data;
          initialize_t8_data_per_element_1d(&data);
          grid_hierarchy.insert(0, element_data[current_index].lmi, data);
          t8_global_productionf ("Initial grid hierarchy data u coeff: %f \n",grid_hierarchy.get(0, element_data[current_index].lmi).u_coeff[0]);
          double verts[3][3] = { 0 };
            t8_forest_element_coordinate (forest, itree,element,  0,
                                  verts[0]);
            t8_forest_element_coordinate (forest, itree,element,  1,
                                  verts[1]);
            t8_forest_element_coordinate (forest, itree,element,  2,
                                  verts[2]);

          A.resize(3,3);
          r.resize(3);
          A(0,0)=verts[0][0];
          A(0,1)=verts[1][0];
          A(0,2)=verts[2][0];
          A(1,0)=verts[0][1];
          A(1,1)=verts[1][1];
          A(1,2)=verts[2][1];
          A(2,0)=1;
          A(2,1)=1;
          A(2,2)=1;
          A.lr_factors(A,r);
          double eckpunkte[6] = {
            verts[0][0], verts[0][1],
            verts[1][0], verts[1][1],
            verts[2][0], verts[2][1]};
          reference_to_physical_t3 (eckpunkte, order_num, xytab_ref, xytab);
          for (int i = 0; i < M_mra; ++i) {
            double quad = 0.;
            for (int order = 0; order < order_num; ++order) {
              double x = xytab[order*2];
              double y = xytab[1+order*2];
              vec tau(3); tau(0) = x; tau(1) = y; tau(2) = 1.;
              A.lr_solve(A, r, tau);
              quad += wtab[order] * F(x,y) * sqrt(1./(2.*volume)) * skalierungsfunktion(i,tau(0),tau(1));
              }
            quad *= volume;
            element_data[current_index].u_coeff[i] = quad;
            t8_global_productionf ("berechnet quad also u coeff: %f \n",quad);
            //Access the data at the specified level and key
            // Now, access the data at level 1, key 42
            t8_data_per_element_1d_gh& data_ref = grid_hierarchy.get(0, element_data[current_index].lmi);  // Access the data by reference
            data_ref.u_coeff[i]= quad;  // Modify the `significant` flag
            t8_global_productionf ("access u coeff: %f \n",grid_hierarchy.get(0, element_data[current_index].lmi).u_coeff[i]);
            }
          }
      }
    }
    T8_FREE(wtab);
    T8_FREE(xytab);
    T8_FREE(xytab_ref);
    return element_data;
  }

  /* The adaptation callback function. This function will be called once for each element
   * and the return value decides whether this element should be refined or not.
   *   return > 0 -> This element should get refined.
   *   return = 0 -> This element should not get refined.
   * If the current element is the first element of a family (= all level l elements that arise from refining
   * the same level l-1 element) then this function is called with the whole family of elements
   * as input and the return value additionally decides whether the whole family should get coarsened.
   *   return > 0 -> The first element should get refined.
   *   return = 0 -> The first element should not get refined.
   *   return < 0 -> The whole family should get coarsened.
   *
   * \param [in] forest       The current forest that is in construction.
   * \param [in] forest_from  The forest from which we adapt the current forest (in our case, the uniform forest)
   * \param [in] which_tree   The process local id of the current tree.
   * \param [in] lelement_id  The tree local index of the current element (or the first of the family).
   * \param [in] ts           The refinement scheme for this tree's element class.
   * \param [in] is_family    if 1, the first \a num_elements entries in \a elements form a family. If 0, they do not.
   * \param [in] num_elements The number of entries in \a elements elements that are defined.
   * \param [in] elements     The element or family of elements to consider for refinement/coarsening.
   */
  int
  t8_mra_bottom_up_init_callback (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t which_tree, t8_locidx_t lelement_id,
                           t8_eclass_scheme_c *ts, const int is_family, const int num_elements, t8_element_t *elements[])
  {
    const struct adapt_data_1d_wavelet_based *adapt_data = (const struct adapt_data_1d_wavelet_based *) t8_forest_get_user_data (forest);
    const t8_element_t *elem = t8_forest_get_element_in_tree (forest_from, which_tree, lelement_id);
    /* You can use T8_ASSERT for assertions that are active in debug mode (when configured with --enable-debug).
     * If the condition is not true, then the code will abort.
     * In this case, we want to make sure that we actually did set a user pointer to forest and thus
     * did not get the NULL pointer from t8_forest_get_user_data.
     */
    int rule=10;
    T8_ASSERT (adapt_data != NULL);
    int order_num;
    double *wtab;
    double *xytab;
    double *xytab_ref;
    order_num = dunavant_order_num(rule);
    wtab = T8_ALLOC (double, order_num);
    xytab = T8_ALLOC (double, 2*order_num);
    xytab_ref = T8_ALLOC (double, 2*order_num);
    mat A;
    vector<int> r;
    dunavant_rule(rule, order_num, xytab_ref, wtab);
    // adapt_data->gamma;
    // adapt_data->C_thr;
    // adapt_data->current_level;
    levelgrid_map<t8_data_per_element_1d_gh>& grid_map = *(adapt_data->grid_map_ptr);

    const t8_eclass_t tree_class = t8_forest_get_tree_class (forest_from, which_tree);
    int num_children = element_get_num_children(tree_class,elem);  // The expected number of children
    t8_element_t *children;
    double *u_coeff_children;
    int *child_order;
    uint64_t *children_lmi;
    children_lmi=T8_ALLOC(uint64_t,num_children);
    children=T8_ALLOC(t8_element_t,num_children);
    u_coeff_children=T8_ALLOC(double,num_children);
    child_order=T8_ALLOC(int,num_children);
    // Call the function to get the children of the element
    element_get_children(tree_class, elem, num_children, children);
    bool significant=0;
    // 1D Wavelet-Based Grid Data
    element_data_children = T8_ALLOC (struct t8_data_per_element_1d_gh,num_children);
    uint64_t parent_lmi=calculate_lmi((uint64_t)t8_forest_global_tree_id (forest_from, which_tree),elem,tree_class);
    int first_parent=grid_map.get(element_get_level(elem), parent_lmi).first;
    int second_parent=grid_map.get(element_get_level(elem), parent_lmi).second;
    int third_parent=grid_map.get(element_get_level(elem), parent_lmi).third;
    for (ichild = 0; ichild < num_children; ichild++){
      double verts[3][3] = { 0 };
      int first_copy = first_parent;
      int second_copy = second_parent;
      int third_copy = third_parent;
      invert_order(&first_copy, &second_copy,&third_copy);
      child_order[ichild]=get_correct_order_children((((t8_dtri_t *)elem)->type),ichild,first_parent,second_parent, third_parent);
      uint64_t child_lmi= get_jth_child_lmi_binary(parent_lmi, child_order[ichild]);
      first_copy = first_parent;
      second_copy = second_parent;
      third_copy = third_parent;
      double volume=t8_forest_element_volume (forest_from, which_tree, children[ichild]);
      get_point_order(&first_copy, &second_copy, &third_copy,
                       t8_dtri_type_cid_to_beyid[((t8_dtri_t *)elem)->type][ichild]);
      t8_forest_element_coordinate (forest_from, which_tree,children[ichild],  0,
                            verts[first_copy]);
      t8_forest_element_coordinate (forest_from, which_tree,children[ichild],  1,
                            verts[second_copy]);
      t8_forest_element_coordinate (forest_from, which_tree,children[ichild],  2,
                            verts[third_copy]);

      A.resize(3,3);
      r.resize(3);
      A(0,0)=verts[0][0];
      A(0,1)=verts[1][0];
      A(0,2)=verts[2][0];
      A(1,0)=verts[0][1];
      A(1,1)=verts[1][1];
      A(1,2)=verts[2][1];
      A(2,0)=1;
      A(2,1)=1;
      A(2,2)=1;
      A.lr_factors(A,r);
      double eckpunkte[6] = {
        verts[0][0], verts[0][1],
        verts[1][0], verts[1][1],
        verts[2][0], verts[2][1]};
      reference_to_physical_t3 (eckpunkte, order_num, xytab_ref, xytab);
      for (int i = 0; i < M_mra; ++i) {
        double quad = 0.;
        for (int order = 0; order < order_num; ++order) {
          double x = xytab[order*2];
          double y = xytab[1+order*2];
          vec tau(3); tau(0) = x; tau(1) = y; tau(2) = 1.;
          A.lr_solve(A, r, tau);
          quad += wtab[order] * adapt_data->F(x,y) * sqrt(1./(2.*volume)) * skalierungsfunktion(i,tau(0),tau(1));
          }
        quad *= volume;
        element_data_children[child_order[ichild]].first=first_copy;
        element_data_children[child_order[ichild]].second=second_copy;
        element_data_children[child_order[ichild]].third=third_copy;
        element_data_children[child_order[ichild]].u_coeff[i]=quad;
        element_data_children[child_order[ichild]].significant=false;
        element_data_children[child_order[ichild]].d_coeff={0};
    }
  }
    // Two scale transformation
    for (int i = 0; i < 3*M_mra; ++i) {
      double d_sum = 0.;
      for (int j = 0; j < M_mra; ++j) {
        double v0=element_data_children[child_order[0]].u_coeff[i];
        double v1=element_data_children[child_order[1]].u_coeff[i];
        double v2=element_data_children[child_order[2]].u_coeff[i];
        double v3=element_data_children[child_order[3]].u_coeff[i];
        d_sum += N0(i,j)*v0;
        d_sum += N1(i,j)*v1;
        d_sum += N2(i,j)*v2;
        d_sum += N3(i,j)*v3;
      }
      grid_map.get(element_get_level(elem), parent_lmi).d_coeff[i]=d_sum;
    }
    for (int i = 0; i < 3*M_mra; ++i) {
      if (abs(grid_map.get(element_get_level(elem), parent_lmi).d_coeff[i]) > sqrt(2.0*volume)*adapt_data->C_thr){
        /* Do not change this element. */
        T8_FREE(element_data_children);
        T8_FREE(u_coeff_children);
        T8_FREE(children_lmi);
        T8_FREE(children);
        T8_FREE(wtab);
        T8_FREE(xytab_ref);
        T8_FREE(xytab);
        T8_FREE(child_order);
        return 1;
      }
  }
    /* Do not change this element. */
    T8_FREE(element_data_children);
    T8_FREE(u_coeff_children);
    T8_FREE(children_lmi);
    T8_FREE(children);
    T8_FREE(wtab);
    T8_FREE(xytab_ref);
    T8_FREE(xytab);
    T8_FREE(child_order);
    return 0;
  }

  // Correct
int main() {
    InitialisiereKoeff(p_mra, M0, M1, M2, M3, N0, N1, N2, N3);
    return 0;
}
