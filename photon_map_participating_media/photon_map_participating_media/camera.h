//
//  camera.h
//  photon_map_participating_media
//
//  Created by kaichen on 10/21/14.
//  Copyright (c) 2014 kaichen. All rights reserved.
//

#include <stdio.h>
#include "vector.h"


#define CAM_POS   Vector(50,48,295.6)
#define CAM_DIR   (Vector(0,-0.042612,-1).normal())

class Camera
{
public:
    Camera(double a, double b)
    {
        w = a;
        h = b;
        pos = CAM_POS;
        dir = CAM_DIR;
    }
    
    void setCamDir(Vector &d);
    void setCamPos(Vector &p);
    
    Vector getCamDir(void);
    Vector getCamPos(void);
    
    void setCamX(void);
    void setCamY(void);
    
    Vector getCamX(void);
    Vector getCamY(void);
    
    Vector getEyePos(Vector &d);
    Vector getEyeDir(double x, double y);
    
    
    void init(void);
private:
    Vector pos;
    Vector dir;
    Vector _x;
    Vector _y;
    
    double w;
    double h;
};