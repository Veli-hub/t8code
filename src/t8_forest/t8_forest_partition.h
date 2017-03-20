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

/** \file t8_forest_partition.h
 * We define the partition routine to partition a forest of trees in this file.
 */

/* TODO: begin documenting this file: make doxygen 2>&1 | grep t8_forest_partition */

#ifndef T8_FOREST_PARTITION_H
#define T8_FOREST_PARTITION_H

#include <t8.h>
#include <t8_forest.h>

T8_EXTERN_C_BEGIN ();
/* TODO: document */
void                t8_forest_partition (t8_forest_t forest);

/** Create the element_offset array of a partitioned forest.
 * \param [in,out]  forest The forest.
 * \a forest must be committed before calling this function.
 */
void                t8_forest_partition_create_offsets (t8_forest_t forest);

/** Create the array of global_first_element ids of a partitioned forest.
 * \param [in,out]  forest The forest.
 * \a forest must be committed before calling this function.
 */
void                t8_forest_partition_create_first_elements (t8_forest_t
                                                               forest);

/** Create the array tree offsets of a partitioned forest.
 * This arrays stores at position p the global id of the first tree of this process.
 * Or if this tree is shared, it stores -(global_id) - 1.
 * \param [in,out]  forest  The forest.
 * \a forest must be committed before calling this function.
 */
void                t8_forest_partition_create_tree_offsets (t8_forest_t
                                                             forest);

T8_EXTERN_C_END ();

#endif /* !T8_FOREST_PARTITION_H! */
