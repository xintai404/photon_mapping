//
//  scene.cpp
//  photon_map_participating_media
//
//

#include "scene.h"

Sphere sph[8] = { // Scene: radius, position, emission, color, material
    Sphere(1e5, Vector( 1e5+1,40.8,81.6), Vector(),Vector(1.0,0.0,0.0),DIFF),//Left
    Sphere(1e5, Vector(-1e5+99,40.8,81.6),Vector(),Vector(0.0,0.0,1.0),DIFF),//Rght
    Sphere(1e5, Vector(50,40.8, 1e5),     Vector(),Vector(0.75,0.75,0.75),DIFF),//Back
    Sphere(1e5, Vector(50,40.8,-1e5+170), Vector(),Vector(),           DIFF),//Frnt
    Sphere(1e5, Vector(50, 1e5, 81.6),    Vector(),Vector(0.75,0.75,0.75),DIFF),//Botm
    Sphere(1e5, Vector(50,-1e5+81.6,81.6),Vector(),Vector(0.75,0.75,0.75),DIFF),//Top
    Sphere(16.5, Vector(73,16.5,47),      Vector(),Vector(1.0,1.0,1.0), SPEC),//Mirr
    Sphere(13.5, Vector(34,36.5,88),      Vector(),Vector(1.0,1.0,1.0), REFR), //Glas
};


// following scene copy from kevin Beason
#if 0
//double R=60;
double R=120;     // radius
double T=30*M_PI/180.;
double D=R/cos(T);     //distance
// double D=60;     //distance
// double R=D*sqrt(2);
double Z=62;
Vector C=Vector(0.275, 0.612, 0.949);
Sphere sph[] = {//Scene: radius, position, emission, color, material
    
    Sphere(R, Vector(50,28,Z)+Vector( cos(T),sin(T),0)*D,    C*6e-2,Vector(1,1,1)*.996, SPEC), //red
    Sphere(R, Vector(50,28,Z)+Vector(-cos(T),sin(T),0)*D,    C*6e-2,Vector(1,1,1)*.996, SPEC), //grn
    Sphere(R, Vector(50,28,Z)+Vector(0,-1,0)*D,              C*6e-2,Vector(1,1,1)*.996, SPEC), //blue
    Sphere(R, Vector(50,28,Z)+Vector(0,0,-1)*R*2*sqrt(2./3.),C*0e-2,Vector(1,1,1)*.996, SPEC), //back
    //  Sphere(1e5, Vec(50,28,Z)+Vec(0,0,1e5+170),   Vec(1,1,1)*0,Vec(1,1,1)*.996, SPEC), //front
    //  Sphere(2*R*2*sqrt(2./3.)-R*2*sqrt(2./3.)/3., Vec(50,28,Z)+Vec(0,0,-R*2*sqrt(2./3.)/3.),   Vec(1,1,1)*0,Vec(1,1,1)*.3333, SPEC), //front
    Sphere(2*2*R*2*sqrt(2./3.)-R*2*sqrt(2./3.)/3., Vector(50,28,Z)+Vector(0,0,-R*2*sqrt(2./3.)/3.),   Vector(1,1,1)*0,Vector(1,1,1)*.5, SPEC), //front
};
#endif
#if 0
double D=50;
double R=40;
Sphere sph[]={// = {//Scene: radius, position, emission, color, material
    Sphere(150, Vector(50+75,28,62), Vector(1,1,1)*0e-3, Vector(1,.9,.8)*.93, REFR),
    Sphere(28,  Vector(50+5,-28,62), Vector(1,1,1)*1e1, Vector(1,1,1)*0, DIFF),
    Sphere(300, Vector(50,28,62), Vector(1,1,1)*0e-3, Vector(1,1,1)*.93, SPEC)
};

#endif

#if 0
Vector Cen(50,40.8,-860);
Sphere sph[] = {//Scene: radius, position, emission, color, material
    // center 50 40.8 62
    // floor 0
    // back  0
    
    Sphere(1600, Vector(1,0,2)*3000, Vector(1,.9,.8)*1.2e1*1.56*2,Vector(), DIFF), // sun
    Sphere(1560, Vector(1,0,2)*3500,Vector(1,.5,.05)*4.8e1*1.56*2, Vector(),  DIFF), // horizon sun2
    //   Sphere(10000,Cen+Vec(0,0,-200), Vec(0.0627, 0.188, 0.569)*6e-2*8, Vec(.7,.7,1)*.25,  DIFF), // sky
    Sphere(10000,Cen+Vector(0,0,-200), Vector(0.00063842, 0.02001478, 0.28923243)*6e-2*8, Vector(.7,.7,1)*.25,  DIFF), // sky
    
    Sphere(100000, Vector(50, -100000, 0),  Vector(),Vector(.3,.3,.3),DIFF), // grnd
    Sphere(110000, Vector(50, -110048.5, 0),  Vector(.9,.5,.05)*4,Vector(),DIFF),// horizon brightener
    Sphere(4e4, Vector(50, -4e4-30, -3000),  Vector(),Vector(.2,.2,.2),DIFF),// mountains
    //  Sphere(3.99e4, Vec(50, -3.99e4+20.045, -3000),  Vec(),Vec(.7,.7,.7),DIFF),// mountains snow
    
    Sphere(26.5,Vector(22,26.5,42),   Vector(),Vector(1,1,1)*.596, SPEC), // white Mirr
    Sphere(13,Vector(75,13,82),   Vector(),Vector(.96,.96,.96)*.96, REFR),// Glas
    Sphere(22,Vector(87,22,24),   Vector(),Vector(.6,.6,.6)*.696, REFR)    // Glas2
};
#endif
#if 0
Vector tc(0.0588, 0.361, 0.0941);
Vector sc = Vector(1,1,1)*.7;
Sphere sph[] = {//Scene: radius, position, emission, color, material
    // center 50 40.8 62
    // floor 0
    // back  0
    //  Sphere(1e5, Vec(50, 1e5+100, 0),  Vec(1,1,1)*1,Vec(),DIFF), //lite
    //  Sphere(1e5, Vec(50, -1e5, 0),  Vec(),Vec(.3,.3,.1),DIFF), //grnd
    //  Sphere(1e5, Vec(50, 1e5+100, 0),  Vec(0.761, 0.875, 1.00)*1.3,Vec(),DIFF),
    //  //lite
    Sphere(1e5, Vector(50, 1e5+130, 0),  Vector(1,1,1)*1.3,Vector(),DIFF), //lite
    Sphere(1e2, Vector(50, -1e2+2, 47),  Vector(),Vector(1,1,1)*.7,DIFF), //grnd
    
    Sphere(1e4, Vector(50, -30, 300)+Vector(-sin(50*M_PI/180),0,cos(50*M_PI/180))*1e4, Vector(), Vector(1,1,1)*.99,SPEC),// mirr L
    Sphere(1e4, Vector(50, -30, 300)+Vector(sin(50*M_PI/180),0,cos(50*M_PI/180))*1e4,  Vector(), Vector(1,1,1)*.99,SPEC),// mirr R
    Sphere(1e4, Vector(50, -30, -50)+Vector(-sin(30*M_PI/180),0,-cos(30*M_PI/180))*1e4,Vector(), Vector(1,1,1)*.99,SPEC),// mirr FL
    Sphere(1e4, Vector(50, -30, -50)+Vector(sin(30*M_PI/180),0,-cos(30*M_PI/180))*1e4, Vector(), Vector(1,1,1)*.99,SPEC),// mirr
    
    
    Sphere(4, Vector(50,6*.6,47),   Vector(),Vector(.13,.066,.033), DIFF),//"tree"
    Sphere(16,Vector(50,6*2+16*.6,47),   Vector(), tc,  DIFF),//"tree"
    Sphere(11,Vector(50,6*2+16*.6*2+11*.6,47),   Vector(), tc,  DIFF),//"tree"
    Sphere(7, Vector(50,6*2+16*.6*2+11*.6*2+7*.6,47),   Vector(), tc,  DIFF),//"tree"
    
    Sphere(15.5,Vector(50,1.8+6*2+16*.6,47),   Vector(), sc,  DIFF),//"tree"
    Sphere(10.5,Vector(50,1.8+6*2+16*.6*2+11*.6,47),   Vector(), sc,  DIFF),//"tree"
    Sphere(6.5, Vector(50,1.8+6*2+16*.6*2+11*.6*2+7*.6,47),   Vector(), sc,  DIFF),//"tree"
};
#endif
#if 0
Sphere sph[] = { // Scene: radius, position, emission, color, material
    Sphere(1e5, Vector( 1e5+1,40.8,81.6), Vector(),Vector(.75,.25,.25),DIFF),//Left
    Sphere(1e5, Vector(-1e5+99,40.8,81.6),Vector(),Vector(.25,.25,.75),DIFF),//Rght
    Sphere(1e5, Vector(50,40.8, 1e5),     Vector(),Vector(.75,.75,.75),DIFF),//Back
    Sphere(1e5, Vector(50,40.8,-1e5+170), Vector(),Vector(),           DIFF),//Frnt
    Sphere(1e5, Vector(50, 1e5, 81.6),    Vector(),Vector(.75,.75,.75),DIFF),//Botm
    Sphere(1e5, Vector(50,-1e5+81.6,81.6),Vector(),Vector(.75,.75,.75),DIFF),//Top
    Sphere(16.5, Vector(73,16.5,47),      Vector(),Vector(1,1,1)*.999, SPEC),//Mirr
    Sphere(16.5, Vector(40,36.5,88),      Vector(),Vector(1,1,1)*.999, REFR), //Glas
    Sphere(600, Vector(50,681.6-.27,81.6),Vector(1,1,1),  Vector(), DIFF) //Lite
    
};
#endif