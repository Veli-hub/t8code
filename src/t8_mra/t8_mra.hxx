// t8 mra hxx
// ich brauche folgendes:
// element data hier ohne Details
// die ganzen FUnktionen wie vorher
// Auswertung MS
// AUswertung SS
//
//
// Bottom Up Init:
// input: cmesh, ctresh, gamma, grid hierarchy, soll forest direkt adapten
// und grid hierarchy Daten eintragen
// thresholding:
// Ctresh, grid hierarchy: Soll Gitterhierarchie und forest anpassen
// Gitterhierarchie nur signifikant setzen
// Prediction:
// grid hierarchy anpassen und forest adaptieren
// nehmen als Input immer forest und grid hierarchy
//
// Element data ohne Detail coefficients
//
// lmi Funktionen

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
struct t8_data_per_element_1d *
  t8_create_element_data (levelgrid_map<t8_data_per_element_1d_gh>& grid_hierarchy, t8_forest_t forest,func F, const int rule, const int max_lev);
