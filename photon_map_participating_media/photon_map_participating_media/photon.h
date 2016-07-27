//
//  photon.h
//  photon_map_participating_media
//
//

#ifndef __photon_map_H__
#define __photon_map_H__


// copy from book
#include "vector.h"
#include "ray.h"


#define MAX_PHOTON_DEPTH     15
class Photon
{
public:
    Vector pos;
    Vector power;
    Vector dir;    // postion, power, direction
    int plane;
};


class NearestPhotons
{
public:
    int max;
    int found;
    int got_heap;
    Vector pos;
    double *dist2;
    Photon **index;
};

class Photon_map
{
public:
    Photon_map() {};
    ~Photon_map();
    
    void init(int max_photon);
    void store(Vector power,  Vector pos,  Vector dir);
    void scale_photon_power(double scale);
    void balance(void);
    void irradiance_estimate(Vector& radiance, Vector& pos, Vector& normal, double max_dist, int nphotons,int flag=0) const;
    void locate_photons(NearestPhotons* np, int index) const;
    void photon_dir(Vector& dir,  Photon *p) const;
    
    void photon_print() const;

private:
    void balance_segment(Photon **pbal, Photon **porg, int index, int start, int end);
    void median_split(Photon **p, int start, int end, int median, int axis);
    
    Photon *photons;
    int stored_photons;
    int half_stored_photons;
    int max_photons;
    int prev_scale;
    
    float bbox_min[3];
    float bbox_max[3];
};
extern void photon_emission1(int photon_num, unsigned short seed);
extern void photon_emission2(int photon_num, unsigned short seed);
extern Vector surface_radiance(const Ray &r, int depth,int photon_cal);
extern Vector volume_radiance(const Ray &r, int depth, int photon_cal);
extern bool photon_surface_interact(Photon_map* map, Ray &r,Vector &power, unsigned short *seed);
#endif 
