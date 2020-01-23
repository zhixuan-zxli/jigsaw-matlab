
/* 
------------------------------------------------------------
 * data-structures for discrete complexes in R^d. 
------------------------------------------------------------
 *
 * This program may be freely redistributed under the 
 * condition that the copyright notices (including this 
 * entire header) are not removed, and no compensation 
 * is received through use of the software.  Private, 
 * research, and institutional use is free.  You may 
 * distribute modified versions of this code UNDER THE 
 * CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE 
 * TO IT IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE 
 * ORIGINAL AUTHOR, BOTH SOURCE AND OBJECT CODE ARE 
 * MADE FREELY AVAILABLE WITHOUT CHARGE, AND CLEAR 
 * NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution 
 * of this code as part of a commercial system is 
 * permissible ONLY BY DIRECT ARRANGEMENT WITH THE 
 * AUTHOR.  (If you are not directly supplying this 
 * code to a customer, and you are instead telling them 
 * how they can obtain it for free, then you are not 
 * required to make any arrangement with me.) 
 *
 * Disclaimer:  Neither I nor: Columbia University, The
 * Massachusetts Institute of Technology, The 
 * University of Sydney, nor The National Aeronautics
 * and Space Administration warrant this code in any 
 * way whatsoever.  This code is provided "as-is" to be 
 * used at your own risk.
 *
------------------------------------------------------------
 *
 * Last updated: 08 December, 2019
 *
 * Copyright 2013-2019
 * Darren Engwirda
 * de2363@columbia.edu
 * https://github.com/dengwirda/
 *
------------------------------------------------------------
 */

#   pragma once

#   ifndef __MESHTYPE__
#   define __MESHTYPE__
 
    namespace mesh 
    {
/*-------------------------- classification of mesh cells */
    char_type EMPTY_tag = +0 ;
    char_type POINT_tag = 10 ;
    char_type EDGE2_tag = 20 ;
    char_type TRIA3_tag = 30 ;
    char_type QUAD4_tag = 40 ;
    char_type TRIA4_tag = 50 ;
    char_type HEXA8_tag = 60 ;
    char_type WEDG6_tag = 70 ;
    char_type PYRA5_tag = 80 ;
    }
 
#   include "containers.hpp"

#   include "hashfunc.hpp"

#   include "mesh_type/mesh_complex_type_k.hpp"

#   include "mesh_type/tria_complex_1.hpp"
#   include "mesh_type/tria_complex_2.hpp"
#   include "mesh_type/tria_complex_3.hpp"

#   endif//__MESHTYPE__



