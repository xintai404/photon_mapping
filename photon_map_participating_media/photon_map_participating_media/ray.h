//
//  ray.h
//  photon_map_participating_media
//
//


#ifndef __RAY_H__
#define __RAY_H__

#include "vector.h"

class Ray
{
public:
    Vector dir;
    Vector ori;
    
    Ray() {}
    Ray(Vector p, Vector d): dir(d),ori(p){}
};




#endif