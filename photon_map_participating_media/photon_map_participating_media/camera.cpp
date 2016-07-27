//
//  camera.cpp
//  photon_map_participating_media
//
//

#include "camera.h"


void Camera::setCamDir(Vector &d)
{
    this->dir = d;
    return;
}
void Camera::setCamPos(Vector &p)
{
    this->pos = p;
    return;
}

Vector Camera::getCamDir(void)
{
    return this->dir;
}

Vector Camera::getCamPos(void)
{
    return this->pos;
}

void Camera::setCamX(void)
{
    double temp = w/h;
    _x= Vector(temp*.5135);
}

void Camera::setCamY(void)
{
    Vector temp = _x.cross(dir);
    temp = temp.normal();
    _y  = temp*0.5135;
}

Vector Camera::getCamX()
{
    return _x;
}

Vector Camera::getCamY()
{
    return _y;
}

Vector Camera::getEyePos(Vector &d)
{
    return pos + d*140;
}


Vector Camera::getEyeDir(double x, double y)
{
    Vector d;
    d = _x * ((x + 0.5) / w - 0.5) + _y * (-(y + 0.5) / h + 0.5);
    return d+ dir;

}

void Camera::init(void)
{
    setCamX();
    setCamY();
    
    return;
}








