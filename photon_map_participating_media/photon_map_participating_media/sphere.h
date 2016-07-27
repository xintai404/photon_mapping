//
//  sphere.h
//  photon_map_participating_media
//
//
#ifndef __SPHERE_H__
#define __SPHERE_H__

#include "ray.h"
#include "utils.h"
#include <math.h>
#include <stdio.h>

class Sphere
{
public:
    Sphere(double r,Vector p,Vector e,Vector c, Material_t m)
    {
        radius = r;
        pos= p;
        col = c;
        ma = m;
        emi = e;
    }
    
    double intersect(const Ray &r, double *near=NULL, double *far=NULL) const
    { // returns distance
        Vector op=pos-r.ori;
        double t;
        double b,delta;
        
        b=r.dir.dot(op);
        delta =b*b-op.dot(op)+radius*radius;
        if (delta < epslion)
            return Infinity;
        else
            delta = sqrt(delta);
        
        if((near != NULL) &&(far != NULL) )
        {
            if((b-delta)<epslion)
                *near = 0;
            else
                *near = b-delta;
            if((b+delta)<epslion)
                *far = 0;
            else
                *far = b+delta;
                
        }
        if((t=b-delta)> epslion)
        {
            return t;
        }
        else
        {
            if((t=b+delta)> epslion)
                return t;
            else
                return Infinity;
        }

        
    }
    
    double radius;
    Vector pos, col;//postion, color
    Vector emi;     // emission color
    Material_t ma;
};



#endif