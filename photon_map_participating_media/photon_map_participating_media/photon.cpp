//
//  photon.cpp
//  photon_map_participating_media
//
//
//photon_map structure is adapted from outside source 

#include <stdio.h>

#include "photon.h"
#include "utils.h"
#include "ray.h"
#include "math.h"
#include "sphere.h"

#include <algorithm>

extern Photon_map volume_photon_map;
extern Sphere sph[];
extern bool intersection(const Ray &r,double &t,int &id, double max=Infinity);

Photon_map global_photon_map;


void Photon_map:: init(int max_photon)
{
    stored_photons = 0;
    prev_scale = 1;
    max_photons = max_photon;
    
    photons = new Photon[max_photons+1];
    
    bbox_min[0]=bbox_min[1]=bbox_min[2] = Infinity;
    bbox_max[0]=bbox_max[1]=bbox_max[2] = epslion;
    
}

Photon_map::~Photon_map()
{
    delete [] photons;
}

void Photon_map::photon_dir(Vector& dir, Photon *p) const
{
    dir = p->dir;
}

void Photon_map::irradiance_estimate(Vector& radiance, Vector& pos, Vector& normal, double max_dist, int nphotons,int flag) const
{
    Vector pdir;
    double tmp;
    NearestPhotons np;
    
    radiance = Vector();
    
    np.dist2 = new double[nphotons+1];
    np.index = new Photon*[nphotons+1];
    
    np.pos = pos;
    np.max = nphotons;
    np.found = 0;
    np.got_heap = 0;
    np.dist2[0] = max_dist*max_dist;
    
    locate_photons(&np, 1);
    /*
    int i = 0;
    while(i < np.found)
    {
        i++;
        Photon *photon = np.index[i];
        printf("np photon : %d\n",i);
        printf("power: %e, %e, %e\n",photon->power.d[0],photon->power.d[1],photon->power.d[2]);
        printf("pos  : %e, %e, %e\n", photon->pos.d[0],photon->pos.d[1],photon->pos.d[2]);
        printf("dir  : %e, %e, %e\n", photon->dir.d[0],photon->dir.d[1],photon->dir.d[2]);
        printf("plane: %d\n", photon->plane);
    }
     */

    for(int i =1; i<= np.found; i++)
    {
        Photon *p = np.index[i];
        photon_dir(pdir,p);
        radiance = radiance+ p->power;
        //printf("power %e, %e %e\n",p->power.d[0],p->power.d[1],p->power.d[2]);
        
    }
    if(flag)// volume irradiance estimate
    {
        double d = sqrt(np.dist2[0]);
        tmp = (3.0f/4/M_PI)/(d*d*d);
    }
    else// surface irradiance estimate
    {
        tmp = (1.0f/M_PI)/(np.dist2[0]);
    }
    
    radiance.d[0] *= (tmp);
    radiance.d[1] *= (tmp);
    radiance.d[2] *= (tmp);
    
    delete [] np.dist2;
    delete [] np.index;
    
}

void Photon_map:: locate_photons(NearestPhotons *np, int index) const
{
    Photon *p = &photons[index];
    double dist1, dist2;
    Vector temp;
    
    if(index < half_stored_photons)
    {
        dist1 = np->pos[p->plane] - p->pos[p->plane];
    
        if(dist1 > epslion)
        {
            locate_photons(np, 2*index+1);
            if(dist1*dist1 < np->dist2[0])
            {
                locate_photons(np, 2*index);
            }
        }
        else
        {
            locate_photons(np, 2*index);
            if(dist1 * dist1 < np->dist2[0])
            {
                locate_photons(np, 2*index+1);
            }
        }
    }
    
    temp = p->pos - np->pos;
    dist2 = temp.length2();
    
    if(dist2 < np->dist2[0])
    {
        // found a photon
        if(np->found < np->max)
        {
            np->found++;
            np->dist2[np->found] = dist2;
            np->index[np->found] = p;
        }
        else
        {
            // bulid maximum heap
            int j, parent;
            if(np->got_heap == 0)
            {
                double dst2 = 0.0f;
                Photon *phot;
                int half_found = np->found>>1;
                for(int k = half_found; k >= 1; k--)
                {
                    parent = k;
                    phot = np->index[k];
                    dst2 = np->dist2[k];
                    while(parent <= half_found)
                    {
                        j = parent*2;
                        if(j<np->found && np->dist2[j]<np->dist2[j+1])
                        {
                            j++;
                        }
                        if(dst2 >= np->dist2[j])
                        {
                            break;
                        }
                        np->dist2[parent] = np->dist2[j];
                        np->index[parent] = np->index[j];
                        parent = j;
                    }
                    np->dist2[parent] = dst2;
                    np->index[parent] = phot;
                }
                np->got_heap = 1;
            }
            // delete maximum heap head node, and insert a new one
            parent = 1;
            j = 2;
            while( j <= np->found)
            {
                if(j < np->found && np->dist2[j] < np->dist2[j+1])
                    j++;
                if(dist2 > np-> dist2[j])
                    break;
                np->dist2[parent] = np->dist2[j];
                np->index[parent] = np->index[j];
                parent = j;
                j += j;
            }
            np->index[parent] = p;
            np->dist2[parent] = dist2;
            np->dist2[0] = np->dist2[1];

        }
    }
}

void Photon_map::store(Vector power, Vector pos, Vector dir)
{
    if(stored_photons > max_photons)
    {
        //printf("photon map is full\n");
        return;
    }
    
    stored_photons++;
    Photon * node = &photons[stored_photons];
    
    node->pos = pos;
    node->power = power;
    node->dir = dir;
    
    bbox_min[0] = MIN_2(node->pos.d[0], bbox_min[0]);
    bbox_min[1] = MIN_2(node->pos.d[1], bbox_min[1]);
    bbox_min[2] = MIN_2(node->pos.d[2], bbox_min[2]);
    
    bbox_max[0] = MAX_2(node->pos.d[0], bbox_max[0]);
    bbox_max[1] = MAX_2(node->pos.d[1], bbox_max[1]);
    bbox_max[2] = MAX_2(node->pos.d[2], bbox_max[2]);
    
    
}


void Photon_map::scale_photon_power(double scal)
{
    for(int i = prev_scale; i <= stored_photons;i++)
    {
        photons[i].power.scale(scal);
    }
    
    prev_scale = stored_photons;
}

void Photon_map::balance(void)
{
    if(stored_photons>1)
    {
        Photon **pa1 = new Photon*[stored_photons+1];
        Photon **pa2 = new Photon*[stored_photons+1];
        
        for(int i = 0; i <= stored_photons; i++)
        {
            pa2[i] = &photons[i];
        }
        balance_segment(pa1,pa2, 1, 1, stored_photons);
        delete [] pa2;

    
        int d, j=1, foo = 1;
        Photon foo_photon= photons[j];
        
        for(int i = 1; i <= stored_photons; i++)
        {
            d = pa1[j] - photons;
            pa1[j] = NULL;
            if(d != foo)
                photons[j] = photons[d];
            else
            {
                photons[j] = foo_photon;
                if(i < stored_photons)
                {
                    for(; foo<=stored_photons;foo++)
                    {
                        if(pa1[foo] != NULL)
                            break;
                    }
                    foo_photon = photons[foo];
                    j = foo;
                }
                continue;
            }
            j = d;
        }
        delete [] pa1;
    }
    half_stored_photons = stored_photons/2;
}

bool cmpx(const Photon* a,const Photon* b)
{
    return a->pos.d[0] < b->pos.d[0];
}

bool cmpy(const Photon* a,const Photon* b)
{
    return a->pos.d[1] < b->pos.d[1];
}

bool cmpz(const Photon* a,const Photon* b)
{
    return a->pos.d[2] < b->pos.d[2];
}

void Photon_map::median_split(Photon **p, int start, int end , int median, int axis)
{
    std::sort(p+start,p+end,axis==0?cmpx:axis==1?cmpy:cmpz);
    /*
    int left = start;
    int right  = end;
    Photon *tmp;
    while(right > left)
    {
        double v = p[right]->pos[axis];
        int i = left -1;
        int j = right;
        
        for(;;)
        {
            while(p[++i]->pos[axis] < v)
                ;
            while(p[--j]->pos[axis] > v && j > left)
                ;
            if( i >= j)
                break;
            tmp = p[i];
            p[i] = p[j];
            p[j] = tmp;
        }
        
        tmp = p[right];
        p[right] = p[i];
        p[i] = tmp;
        
        if(i >=  median)
            right = i - 1;
        if(i <= median)
            left = i + 1;
    }
     */
}


void Photon_map:: balance_segment(Photon **pbal, Photon **porg, int index, int start, int end)
{
    int median = 1;
    double tmp;
    
    while((4*median) <= (end-start+1))
    {
        median += median;
    }
    
    if((3*median) <= (end - start + 1))
    {
        median += median;
        median += start - 1;
    }
    else
    {
        median = end - median + 1;
    }
    
    int axis = 2;
    if((bbox_max[0] - bbox_min[0])>(bbox_max[1]-bbox_min[1])&&(bbox_max[0]-bbox_min[0])>(bbox_max[2]-bbox_min[2]))
    {
        axis = 0;
    }
    else if((bbox_max[1]-bbox_min[1])>(bbox_max[2]-bbox_min[2]))
    {
        axis = 1;
    }
    
    median_split(porg, start, end, median, axis);
    
    pbal[index] = porg[median];
    pbal[index]->plane = axis;
    
    if(median > start)
    {
        if( start < median - 1)
        {
            tmp = bbox_max[axis];
            bbox_max[axis] = pbal[index]->pos[axis];
            balance_segment(pbal, porg, 2*index, start, median - 1);
            bbox_max[axis] = tmp;
        }
        else
        {
            pbal[2*index] = porg[start];
        }
    }
    
    if(median < end)
    {
        if(median + 1 < end)
        {
            tmp = bbox_min[axis];
            bbox_min[axis] = pbal[index]->pos[axis];
            balance_segment(pbal, porg, 2*index+1, median+1, end);
            bbox_min[axis] = tmp;
        }
        else
        {
            pbal[2*index+1] = porg[end];
        }
    }
}

void Photon_map::photon_print(void) const
{
    /*
    printf("num_photon %d\n",max_photons);
    printf("stored_photon %d\n",stored_photons);
    printf("bbmin %e %e %e\n",bbox_min[0],bbox_min[1],bbox_min[2]);
    printf("bbmax %e %e %e\n",bbox_max[0],bbox_max[1],bbox_max[2]);
     */
    int i = 0;
    Photon photon;
    while(i < stored_photons)
    {
        i++;
        photon = photons[i];
        printf("photon%d\n",i);
        printf("power: %e, %e, %e\n",photon.power.d[0],photon.power.d[1],photon.power.d[2]);
        printf("pos  : %e, %e, %e\n", photon.pos.d[0],photon.pos.d[1],photon.pos.d[2]);
        printf("dir  : %e, %e, %e\n", photon.dir.d[0],photon.dir.d[1],photon.dir.d[2]);
       // printf("plane: %d\n", photon.plane);
    }
}

void kd_tree_test()
{
    Photon_map photon_m;
    Vector power(1.0,1.0,1.0);
    Vector dir(1.0,1.0,1.0);
    Vector pos,irrad,normal;
    int i;
    photon_m.init(100);

    // store photons
    for ( i = 0; i < 100; i++ ) {
        pos.d[0] = (double)i;
        pos.d[1] = (double)i;
        pos.d[2] = (double)i;
        
        photon_m.store(power,pos,dir);
    }
    
    
    // balance kd-tree
    photon_m.balance();
    
    // estimate irradiance
    pos.d[0] = 4.0;
    pos.d[1] = 4.0;
    pos.d[2] = 4.0;
    normal.d[0] = -1.0;
    normal.d[1] = -1.0;
    normal.d[2] = -1.0;
    double max_dist = 25.00;
    int nphotons = 10;
    
    NearestPhotons np;
    np.dist2 = new double[101];
    np.index = new Photon*[101];
    
    np.pos = pos;
    np.max = 10;
    np.found = 0;
    np.got_heap = 0;
    np.dist2[0] = max_dist*max_dist;
    
    photon_m.locate_photons(&np, 1);
    i = 0;
    while(i < np.found)
    {
        i++;
        Photon *photon = np.index[i];
        printf("np photon : %d\n",i);
        printf("power: %e, %e, %e\n",photon->power.d[0],photon->power.d[1],photon->power.d[2]);
        printf("pos  : %e, %e, %e\n", photon->pos.d[0],photon->pos.d[1],photon->pos.d[2]);
        printf("dir  : %e, %e, %e\n", photon->dir.d[0],photon->dir.d[1],photon->dir.d[2]);
        printf("plane: %d\n", photon->plane);
    }
    

    photon_m.irradiance_estimate(irrad,pos,normal,max_dist,nphotons);
    printf("dist2 %f\n", np.dist2[0]);
    printf(" irrad %f %f %f\n",irrad[0], irrad[1], irrad[2]);

}



Vector refl_dir_get(Ray r, Vector norm)
{
    Vector reflection;
    reflection = r.dir- norm*2.0*norm.dot(r.dir);
    return reflection;
}



Vector cosin_dis_sample(Vector &n, double r1, double r2)
{
    double c=sqrt(r2);
    Vector z=n, x, y;
    get_ortho_aixs(x,y,z);
    Vector d= (x*cos(r1)*c + y*sin(r1)*c + z*sqrt(1-r2)).normal();

    return d;
}


void photon_trace(Ray r, int depth, Vector power,unsigned short *seed)
{
    double t,p;
    int id;
    Vector hit, norm, amb;
    Vector nl;

    if(!intersection(r,t,id)||(++depth>5))return;
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
        //diffuse direction
        d =cosin_dis_sample(nl, r1, r2);
        //store photon into photon map
        global_photon_map.store(power,hit,r.dir);

        Ray diff_ray(hit,d);
        if (erand48(seed)>p) return;// absorbation
        else// diffuse
        {
            photon_trace(diff_ray,depth, amb.mul(power)*(1./p), seed);
        }
    }
    else if (obj.ma == SPEC)
    {
        Vector reflect_dir;
        reflect_dir = refl_dir_get(r,norm);
        Ray specular_ray(hit,reflect_dir);
        //Vector spec_factor = amb.mul(power);
        photon_trace(specular_ray, depth, amb.mul(power),seed);
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
            return photon_trace(rel_r,depth,power,seed);
        
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
            photon_trace(rel_r,depth,amb.mul(power),seed);
        else//refraction
            photon_trace(ref_r,depth,amb.mul(power),seed);
    }
}

void accept_rej_sample(double &x, double &y, double &z)
{
    unsigned short sed[3] = {0,0,1000};
    do
    {
        x = erand48(sed)-.5;
        y = erand48(sed)-.5;
        z = erand48(sed)-.5;
        
    }while(x*x+y*y+z*z>1);
    
    return;
}

void sample_light_direction(double &x, double &y, double &z, unsigned short *sed)
{
    double phi,theta;
    phi   =2.*M_PI*erand48(sed);
    theta =2.*acos(sqrt(1.-erand48(sed)));
    
    x = cos(phi)*sin(theta);
    y = cos(theta);
    z = sin(phi)*sin(theta);
}

void photon_emission2(int photon_num, unsigned short seed)
{
    int ne = 0;
    double x,y,z;
    unsigned short sed[3] = {0,0,seed};
    Ray r;
    Vector power;
    
    global_photon_map.init(photon_num*15);
    
    while(ne < photon_num)
    {
        ne++;
        sample_light_direction(x,y,z,sed);
        r.dir = Vector(x,y,z);
        power = Vector(1000,1000,1000)*18;
        
        //radius 10
        r.ori.d[0] = 50+2*cos(2*M_PI*erand48(sed));
        r.ori.d[2] = 85+2*sin(2*M_PI*erand48(sed));
        r.ori.d[1] = 80.6;
        //r.ori=Vector(50,80.6,85);
        
        photon_trace(r,0,power,sed);
    }
    global_photon_map.scale_photon_power(1./photon_num);
    global_photon_map.balance();
    //global_photon_map.photon_print();
    
}


Vector surface_radiance(const Ray &r, int depth, int photon_cal)
{
    double t;                               // distance to intersection
    int id=0;                               // id of intersected object
    if((++depth>=15)||!intersection(r,t,id))
        return Vector();
    const Sphere &obj = sph[id];        // the hit object
    Vector hit=r.ori+r.dir*t;               // hit position
    Vector n=(hit-obj.pos).normal();        // hit position normal
    Vector nl=n.dot(r.dir)<0?n:n*-1;
    Vector amb=obj.col;

    if (obj.ma == DIFF){//diffuse
        Vector irridance;
        global_photon_map.irradiance_estimate(irridance, hit, n, MAX_2(1,photon_cal/20), photon_cal);
        irridance = irridance*(1.0);       //scale irridance
        return obj.emi+amb.mul(irridance);
    }

    else if (obj.ma == SPEC)//reflection
    {
        Vector rel_dir = refl_dir_get(r,n);
        return obj.emi+amb.mul((surface_radiance(Ray(hit,rel_dir),depth,photon_cal)));
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
                return obj.emi+amb.mul(surface_radiance(rel_r,depth,photon_cal));
        
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
        return  obj.emi+((surface_radiance(rel_r,depth,photon_cal)*fresnel).mul(amb) + surface_radiance(ref_r,depth,photon_cal)*(1-fresnel)).mul(amb);
    }
}



