//
//  vector.h
//  photon_map_participating_media
//
//  Created by kaichen on 10/7/14.
//  Copyright (c) 2014 kaichen. All rights reserved.
//

#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <math.h>
#include <iostream>
#include <cassert>


class Vector
{
public:
    Vector(const Vector &l)
    {
        d[0]=l.d[0];
        d[1]=l.d[1];
        d[2]=l.d[2];
    }
    Vector(double a = 0, double b = 0, double c = 0)
    {
        d[0] = a;
        d[1] = b;
        d[2] = c;
    }
    
    Vector operator+(const Vector &l) const
    {
        return Vector(d[0]+l.d[0], d[1]+l.d[1], d[2]+l.d[2]);
    }
    
    Vector operator+(double l) const
    {
        return Vector(d[0] + l, d[1] + l, d[2] + l);
    }
    
    Vector operator-(const Vector &l) const
    {
        return Vector(d[0]-l.d[0], d[1]-l.d[1], d[2]-l.d[2]);
    }

    Vector operator-(double l) const
    {
        return Vector(d[0] - l, d[1] - l, d[2] - l);
    }
    
    Vector& operator=(const Vector& l)
    {
        d[0]=l.d[0];
        d[1]=l.d[1];
        d[2]=l.d[2];
        
        return *this;
    }
    Vector& operator+=(const Vector &l)
    {
        d[0] +=l.d[0];
        d[1] +=l.d[1];
        d[2] +=l.d[2];
        
        return *this;
    }
    
    Vector& operator-=(const Vector &l)
    {
        d[0] -=l.d[0];
        d[1] -=l.d[1];
        d[2] -=l.d[2];
        
        return *this;
    }
    
    Vector operator*(double l) const
    {
        return Vector(d[0] * l, d[1] * l, d[2] * l);
    }
    
    Vector mul(const Vector &l) const
    {
        return Vector(d[0] * l.d[0], d[1] * l.d[1] , d[2] * l.d[2]);
    }
    
    
    double length()
    {
        return sqrt(d[0]*d[0]+ d[1]*d[1] + d[2]*d[2]);
    }
    
    double length2()
    {
        return (d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
    }
    
    Vector normal()
    {
        return (*this) * (1.0 / length());
    }
    
    double dot(const Vector &b) const
    {
        return d[0] * b.d[0] + d[1] * b.d[1] + d[2] * b.d[2];
    }
    
    Vector cross(const Vector &b) const
    {
        return Vector(d[1]*b.d[2]-d[2]*b.d[1], d[2]*b.d[0]-d[0]*b.d[2], d[0]*b.d[1]-d[1]*b.d[0]);
    }
    Vector scale(double i)
    {
        d[0] *= i;
        d[1] *= i;
        d[2] *= i;
        return *this;
    }
    
    double operator[] (const size_t index) const
    {
        assert(0<=index && index<=3);
        return d[index];
    }
    
    double d[3];
};




#endif
