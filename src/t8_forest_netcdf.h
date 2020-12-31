/*
  This file is part of t8code.
  t8code is a C library to manage a collection (a forest) of multiple
  connected adaptive space-trees of general element classes in parallel.

  Copyright (C) 2015 the developers

  t8code is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  t8code is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with t8code; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/* TODO: document this file */

#ifndef T8_FOREST_NETCDF_H
#define T8_FOREST_NETCDF_H

#include <t8_forest.h>
#include <t8_netcdf.h>

T8_EXTERN_C_BEGIN ();

void                t8_forest_write_netcdf (t8_forest_t forest,
                                            const char *file_prefix,
                                            const char *file_title, int dim,
                                            int num_extern_netcdf_vars,
                                            t8_netcdf_variable_t *
                                            ext_variables[]);

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_NETCDF_H */
