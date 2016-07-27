//
//  intersection.cpp
//  photon_map_participating_media
//

//

#include <stdio.h>

#include "ray.h"
#include "sphere.h"
#include "utils.h"
#include "scene.h"




bool intersection(const Ray &r,double &t,int &id, double max=Infinity)
{
    int n = sizeof(sph) / sizeof(Sphere);
    double d, inf = max; t = max;
    for(int i=0;i<n;i++)
    {
        d=sph[i].intersect(r);
        if(d<t)
        {
            t=d;
            id=i;
        }
    }
    return t<inf;
}