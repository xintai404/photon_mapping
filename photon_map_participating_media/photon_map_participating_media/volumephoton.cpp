//
//  participatingmeida.cpp
//  photon_map_participating_media
//
//

#include <math.h>
#include "photon.h"
#include "ray.h"
#include "vector.h"
#include "sphere.h"
#include "noise.h"


//#define HOMO

//#define HETE

//#define ISO
//#define ANISO

#if HETE

#define MARCH_STEP   1.4
const double sigma_s = 0.015;
const double sigma_a = 0.009;

#endif

#if HOMO

#define MARCH_STEP   3
const double sigma_s = 0.009;
const double sigma_a = 0.006;
#endif

double schlick_k = 0.8;
double hg_g      = -0.5;

const double sigma_t = sigma_s+sigma_a;

Photon_map volume_photon_map;

extern Photon_map global_photon_map;
extern Sphere sph[];
extern bool intersection(const Ray &r,double &t,int &id, double max=Infinity);
extern void sample_light_direction(double &x, double &y, double &z, unsigned short *sed);
extern Vector refl_dir_get(Ray r, Vector norm);
extern Vector cosin_dis_sample(Vector &n, double r1, double r2);



Sphere patipating_medium(200, Vector(50,40,85), Vector(), Vector(), DIFF);


double get_turbulance(Vector h)
{
    double p[3];
    p[0] = h.d[0];
    p[1] = h.d[1];
    p[2] = h.d[2];
    
    return abs(noise3(p));
    //return abs(PerlinNoise3D(p[0],p[1],p[2],2,2,16));
}

double sample_homo_dist(double num, double t)
{
    //d = -log(num)/sigma_t
    return -log(num)/t;
}


double sample_hete_dist(double num, Ray &r, double t)
{
    double step = 0.05;
    double dist = 0.0f;
    double cumlativ_step = step;
    double init = t* step;
    double dist_para = -log(num);
    double t_var = t;
    Vector pos = r.ori + r.dir*step;
    while(dist_para > init)
    {
        cumlativ_step +=step;
        pos = pos + r.dir*step;
        t_var = t + get_turbulance(pos);
        init += t_var *step;
    }
    dist = cumlativ_step;
    return dist;
}

void test_sample_hete_dist()
{
    unsigned short seed[]={0,0,222};
    Ray r(Vector(0,0,0), Vector(1,1,1));
    for(int i = 0; i < 100; i++)
    {
        printf("dist %f\n", sample_hete_dist(erand48(seed),r,sigma_t));
    }
}
void test_sample_dist()
{
    unsigned short seed[]={0,0,222};
    for(int i = 0; i < 100; i++)
    {
        printf("dist %f\n", sample_homo_dist(erand48(seed),sigma_t));
    }
}


Vector sample_homo_scatter_dir(Vector &z, double n1, double n2)
{
    double t = 2*M_PI*n1;

    Vector x,y;
    get_ortho_aixs(x,y,z);
    Vector dir = x*cos(t)*sqrt(n2)+y*sin(t)*sqrt(n2)+z*sqrt(1-n2);
    return dir;
}

Vector sample_schlick_phase_dir( double n1, double n2)
{
    double t = 2*M_PI*n1;

    // costheta = (2e+k-1)/(2k*e-k+1)
    double cost = (2*n2+schlick_k-1)/(2*schlick_k*n2-schlick_k+1);
    double sint = sqrt(1- cost*cost);
    
    return Vector(cos(t)*cost, sin(t)*cost, sint);
}

Vector sample_hg_phase_dir(double n1, double n2)
{
    double t = 2*M_PI*n1;
    
    
    //costheta =(1/abs(g))*(1+g*g-(1-g*g)^2/(1-g+2g*n2)^2)
    double temp = (1-hg_g*hg_g)/(1-hg_g+2*hg_g*n2);
    double cost = (1./abs(hg_g))*(1+hg_g*hg_g-temp*temp);
    double sint = sqrt(1-cost*cost);
    return Vector(cos(t)*cost, sin(t)*cost, sint);
}

void sample_light_position(double &x, double &z,unsigned short *seed)
{
    double e = erand48(seed);
    double t = 2*M_PI*e;
    
    x = cos(t);
    z = sin(t);
}



//photon interact with surface for participating media
bool photon_surface_interact(Ray &r,Vector &power, unsigned short *seed)
{
    double t,p;
    int id;
    Vector hit, norm, amb;
    Vector nl;
    
    
    if(!intersection(r,t,id))return 1;
    Sphere obj = sph[id];
    
    
    amb=obj.col;
    //hit point pos
    hit=r.ori+r.dir*t;
    // hit point normal direction
    norm=(hit-obj.pos).normal();
    
    // inside flag
    nl=norm.dot(r.dir)<0?norm:norm*-1.0;
    p = (amb.d[0]+amb.d[1] + amb.d[2])/3;
    if (obj.ma == DIFF)
    {
        double r1=2.0*M_PI*erand48(seed);
        double r2=erand48(seed);
        Vector d;

        d =cosin_dis_sample(nl, r1, r2);
        Ray diff_ray(hit,d);
        if (erand48(seed)>p) return 1;// absorbation
        else// diffuse
        {
            r = diff_ray;
            power = power;
            return 0;
        }
    }
    else if (obj.ma == SPEC)
    {
        Vector reflect_dir;
        reflect_dir = refl_dir_get(r,norm);
        Ray specular_ray(hit,reflect_dir);
        
        r = specular_ray;
        power = power;
        return 0;
    }
    else
    {
        Vector reflect_dir = refl_dir_get(r,norm);
        Ray rel_r(hit,reflect_dir);
        /* ray from inside or outside */
        bool into = (norm.dot(nl)>0.0)?1:0;
        /* air refractive index*/
        double na = AIR_REFR;
        /* glass refractive index */
        double ng = GLASS_REFR;
        
        double ratio = into?na/ng:ng/na;
        double cos_theta_i = r.dir.dot(nl);
        
        double cos_theta_r,tmp;
        if( (tmp =1-ratio*ratio*(1-cos_theta_i*cos_theta_i))<0)
        {
            r = rel_r;
            power = power;
            return 0;
        }
        
        cos_theta_r = sqrt(tmp);
        // refraction direction
        Vector ref_d = (r.dir*ratio - norm*((into?1:-1)*(cos_theta_r + cos_theta_i*ratio))).normal();
        double a=ng-na;
        double b=ng+na;
        double f0=a*a/(b*b);
        double c = 1-(into?-cos_theta_i:ref_d.dot(norm));
        double fresnel=f0+(1-f0)*c*c*c*c*c;
        
        Ray ref_r(hit,ref_d);
        //Russian Roulette
        if(erand48(seed)<fresnel)//reflection
        {
            r = rel_r;
            power = power;
            return 0;
        }
        else//refraction
        {
            r = ref_r;
            power = power;
            return 0;
        }
    }
}

void volume_photon_tracing(Ray r, double depth, Vector power, unsigned short *seed,int refra_flg)
{
    int id;
    double t;
    double near, far;
    Vector hit;
    if(++depth>15) return;
    
    bool meet_medium = patipating_medium.intersect(r, &near, &far);
    // homogenerous media
    if(meet_medium && !refra_flg)
    {
        double dist;
#if HETE // hete
        Ray new_r;
        new_r.ori = r.ori+r.dir*near;
        new_r.dir = r.dir;
        dist = sample_hete_dist(erand48(seed), new_r, sigma_t);
        
        
        
        hit = r.ori +r.dir*near + r.dir*dist;
        // recalculate sigma_s and sigma_t
        double cur_sigma_s, cur_sigma_t, cur_sigma_a;
        cur_sigma_s = sigma_s+ get_turbulance(hit);
        cur_sigma_t = sigma_t+ get_turbulance(hit)*6;
        double p = cur_sigma_s/cur_sigma_t;
#endif
        
#if HOMO // homo
        dist = sample_homo_dist(erand48(seed),sigma_t);
        hit = r.ori +r.dir*near + r.dir*dist;
        double p = sigma_s/sigma_t;
#endif
        
        if(!intersection(r,t,id,near+dist))
        {
            volume_photon_map.store(power,hit,r.dir);
            //interact with media
            if(erand48(seed)<p)
            {
                double n1 = erand48(seed);
                double n2 = erand48(seed);
#if ISO                //isotropic
                Vector scatter_dir = sample_homo_scatter_dir(r.dir, n1,n2);
#endif
                
#if ANISO               //anisotropic
                Vector scatter_dir = sample_schlick_phase_dir(n1,n2);
                //Vector scatter_dir = sample_hg_phase_dir(n1,n2);
#endif
                
                Ray scatter_r = Ray(hit,scatter_dir);
                volume_photon_tracing(scatter_r, depth, power*(1./p), seed,0);
                return;
            }
            else
            {
                //absorb
                //double absorb = exp(-sigma_a*dist);
                //r.ori = hit;
                //volume_photon_tracing(r, depth, power*absorb, seed);
                return;
            
            }
        }
        else
        {
            Sphere obj = sph[id];
            if(obj.ma ==REFR)
            {
                photon_surface_interact(r, power, seed);
                if(refra_flg) power = power * 1.5;
                refra_flg = refra_flg?0:1;
                volume_photon_tracing(r, depth, power,seed,refra_flg);
                return;
            }
        }
    }
    else
    {
        if(intersection(r,t,id))
        {
            Sphere obj = sph[id];
            if(obj.ma ==REFR)
            {
                photon_surface_interact(r, power, seed);
                if(refra_flg) power = power*1.5;
                refra_flg = refra_flg?0:1;
                volume_photon_tracing(r, depth, power,seed,refra_flg);
                return;
            }
        }
    }

}

void photon_emission1(int photon_num, unsigned short seed)
{
    int ne = 0;
    double x,y,z;
    unsigned short sed[3] = {0,0,seed};
    Ray r;
    Vector power;
    double p_x, p_z;
    
    volume_photon_map.init(photon_num*30);
    while(ne < photon_num)
    {
        ne++;
        sample_light_direction(x,y,z,sed);
        sample_light_position(p_x,p_z,sed);
        r.dir = Vector(x,y,z);
#if HETE
        power = Vector(1000,1000,1000)*6;
#endif
        
#if HOMO
        power = Vector(1000,1000,1000)*10;
#endif
        
        //radius
        r.ori.d[0] = 50+2*p_x;
        r.ori.d[2] = 85+2*p_z;
        r.ori.d[1] = 80;
        //r.ori=Vector(50,80.6,85);
        
        volume_photon_tracing(r,0,power,sed,0);
    }
    
    volume_photon_map.scale_photon_power(1./photon_num);
    volume_photon_map.balance();
    
    //volume_photon_map.photon_print();
}


void volume_ray_march(const Ray &r, double near, double far, double step, Vector &power, int photon_cal)
{
    
    if((near + step) > far)
    {
        step = far - near;
    }

    for(double d = far-step; d>= near; d-=step)
    {
        Vector in_scatter;
        Vector n;
        Vector pos = r.ori + r.dir*(d);
        
#if HETE  // non-homo
        // reculate sigma_a and sigma_t
        double cur_sigma_s, cur_sigma_t, cur_sigma_a;
        cur_sigma_s = sigma_s+ get_turbulance(pos);
        cur_sigma_t = sigma_t + get_turbulance(pos)*6;
        double absorb = exp(-step*cur_sigma_t);
        volume_photon_map.irradiance_estimate(in_scatter, pos, n, 10, 50,1);
        power = power*absorb + in_scatter*step*(cur_sigma_t/cur_sigma_t);
#endif
        
        
#if  HOMO  //normal homo
        double absorb = exp(-step*sigma_t);
#if !VISUAL
        volume_photon_map.irradiance_estimate(in_scatter, pos, n, 100, 200,1);
#else
        volume_photon_map.irradiance_estimate(in_scatter, pos, n, MAX_2(1,photon_cal/20), photon_cal,1);
#endif
        power = power*absorb + in_scatter*step*(sigma_s/sigma_t);
#endif
    }
    
}


Vector volume_radiance(const Ray &r, int depth, int photon_cal)
{
    double t;                               // distance to intersection
    int id=0;                               // id of intersected object
    if((++depth>=5)||!intersection(r,t,id))
        return Vector();
    const Sphere &obj = sph[id];        // the hit object
    Vector hit=r.ori+r.dir*t;               // hit position
    Vector n=(hit-obj.pos).normal();        // hit position normal
    Vector nl=n.dot(r.dir)<0?n:n*-1;
    Vector amb=obj.col;
    Vector power;
    double step_random = random()*1.0f/RAND_MAX+0.7f;
    
    if (obj.ma == DIFF){//diffuse
        Vector irridance;
#if !VISUAL
        global_photon_map.irradiance_estimate(irridance, hit, n, Infinity, photon_cal,0);
#else
        global_photon_map.irradiance_estimate(irridance, hit, n, MAX_2(1,photon_cal/20), photon_cal,0);
#endif
        irridance = irridance*(1.0);       //scale irridance
        

        //return obj.emi+amb.mul(irridance);
        power = obj.emi+amb.mul(irridance);
        volume_ray_march(r, 0, t,MARCH_STEP*step_random,power,photon_cal);
        return power;
    }
    else if (obj.ma == SPEC)//reflection
    {
        Vector rel_dir = refl_dir_get(r,n);
        power= obj.emi+amb.mul((volume_radiance(Ray(hit,rel_dir),depth,photon_cal)));
        volume_ray_march(r, 0, t,MARCH_STEP*step_random,power,photon_cal);
        return power;
    }
    else// reflection & refraction
    {
        Vector rel_dir = refl_dir_get(r,n);
        Ray rel_r(hit,rel_dir);
        bool into = (n.dot(nl)>0.0)?1:0;
        double na = AIR_REFR;  //air refractive index
        double ng = GLASS_REFR;   //gloass refractive index
        double ratio = into?na/ng:ng/na;
        double cos_theta_i = r.dir.dot(nl);
        
        double cos_theta_r,tmp;
        if( (tmp =1-ratio*ratio*(1-cos_theta_i*cos_theta_i))<0)
        {
            power = obj.emi+amb.mul(volume_radiance(rel_r,depth,photon_cal));
            volume_ray_march(r, 0, t,MARCH_STEP*step_random,power,photon_cal);
            return power;
        }
        
        cos_theta_r = sqrt(tmp);
        // refraction direction
        Vector ref_d = (r.dir*ratio - n*((into?1:-1)*(cos_theta_r +cos_theta_i*ratio))).normal();
        
        //schlick fresnel approximation F = f0+(1-f0)Pow((1-N*R),5)
        double a=ng-na;
        double b=ng+na;
        double f0=a*a/(b*b);
        double c = 1-(into?-cos_theta_i:ref_d.dot(n));
        double fresnel=f0+(1-f0)*c*c*c*c*c;
        Ray ref_r(hit,ref_d);

        power = obj.emi+((volume_radiance(rel_r,depth,photon_cal)*fresnel).mul(amb) + volume_radiance(ref_r,depth,photon_cal)*(1-fresnel)).mul(amb);
        volume_ray_march(r, 0, t,MARCH_STEP*step_random,power,photon_cal);
        return power;
    }
}






















