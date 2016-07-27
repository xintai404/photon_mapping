//
//  utils.cpp
//  photon_map_participating_media
//
//

#include "math.h"
#include "Vector.h"

int tone_map(double x)
{
    //tone mapping
    return int(pow(1-exp(-x),1/2.2)*255+.5);
}

void get_ortho_aixs(Vector &x, Vector &y, Vector &z)
{
    x=((fabs(z.d[0])>.1?Vector(0,1):Vector(1)).cross(z)).normal();
    y=z.cross(x);
    return;
}