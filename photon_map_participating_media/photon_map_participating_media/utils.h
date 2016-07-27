//
//  utils.h
//  photon_map_participating_media
//
//

#ifndef __UTILS_H__
#define __UTILS_H__


#include "Vector.h"

#define epslion  1e-5f
#define Infinity 1e20f

#define AIR_REFR     1.0
#define GLASS_REFR   1.6

#define MAX_2(a,b)   (a)>(b)?(a):(b)

#define MIN_2(a,b)   (a)>(b)?(b):(a)

#define MAX_3(a,b,c) MAX_2(MAX_2(a,b),c)

#define MIN_3(a,b,c) MIN_2(MIN_2(a,b),c)


enum Material_t
{
    DIFF,
    SPEC,
    REFR
};  // material types


extern int tone_map(double x);
extern void get_ortho_aixs(Vector &x, Vector &y, Vector &z);
#endif
