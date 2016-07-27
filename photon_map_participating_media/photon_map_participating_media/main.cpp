//
//  main.cpp
//  photon_map_participating_media
//
//

#include <iostream>
#include <glut/glut.h>
#include <OpenGL/OpenGL.h>
#include <math.h>   
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "utils.h"
#include "photon.h"
#include "ray.h"
#include "vector.h"

#include "scene.h"
#include "camera.h"

#include "noise.h"


#define MAX(x, y) ((x > y) ? x : y)
#define HEIGHT  384
#define WIDTH   512
Vector *d;
int samples;
unsigned char *pixel;


void my_display(void)
{
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glClear(GL_COLOR_BUFFER_BIT);
    glDrawPixels(WIDTH, HEIGHT, GL_RGBA, GL_UNSIGNED_BYTE, pixel);
    //glFlush();
    glutSwapBuffers();
}

void transfer_to_pixels(Vector *c, unsigned char* pixel)
{
    int w = WIDTH, h =HEIGHT;
    for(int i =0; i <w; i++)
        for(int j = 0; j < h; j++)
        {
            if(c[i+j*w].d[0] > 1) c[i+j*w].d[0] = 1.0;
            if(c[i+j*w].d[1] > 1) c[i+j*w].d[1] = 1.0;
            if(c[i+j*w].d[2] > 1) c[i+j*w].d[2] = 1.0;
            int index = (i+j*w)*4;
            pixel[index+0] = (c[i+j*w].d[0])*255;
            pixel[index+1] = (c[i+j*w].d[1])*255;
            pixel[index+2] = (c[i+j*w].d[2])*255;
            pixel[index+3] = 255;
        } 
}

void idle(void)
{
    static int photon_cal = 1;
    static int fir = 1;
    Camera cam(WIDTH,HEIGHT);
    cam.init();
    if(photon_cal < samples)
    {
#if MEDIUM
        photon_emission1(photon_cal*200, 110);
#endif
        photon_emission2(photon_cal*200, 210); // normal photon mapping
        for (int y=0; y<HEIGHT; y++)
        {
            fprintf(stdout, "\nphoton estimation %4.2f%%", 100.0*y/(HEIGHT-1));
            for (int x=0; x<WIDTH; x++)
            {
                
                int i = x + (HEIGHT-y-1)*WIDTH;
                // calculate direction from eye to scene
                Vector eyeDir = cam.getEyeDir(x,y);
                Vector eyePos = cam.getEyePos(eyeDir);
#if MEDIUM
                d[i] =volume_radiance(Ray(eyePos,eyeDir.normal()),0, MAX(1,photon_cal));
#else
                d[i] = surface_radiance(Ray(eyePos,eyeDir.normal()),0,MAX(1,photon_cal/10)); // normal photon mapping
#endif
            }
        }
        transfer_to_pixels(d,pixel);
        glutPostRedisplay();
        photon_cal *=2;
    }
    
    return;
}


int main(int argc, char *argv[]) {
    int w=WIDTH, h=HEIGHT;
    if(argc>=2)
    {
        samples = MAX(atoi(argv[1]),5000);
    }
    else
    {
        samples = 5000;
    }

    Camera cam(w,h);
    cam.init();

    d = new Vector[w*h];
    pixel=new unsigned char[w*h*4];
    
    
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
    glutInitWindowSize(w, h);
    glutCreateWindow("Photon mapping");

    
#if MEDIUM
    photon_emission1(samples*100, 110);// volum photon mapping
#endif
    
#if !VISUAL
    photon_emission2(samples*100, 210); // normal photon mapping
#endif
    fprintf(stderr,"\r photon emission complete");
    
    for (int y=0; y<h; y++)
    {
        fprintf(stdout, "\nphoton estimation %4.2f%%", 100.0*y/(h-1));
        for (int x=0; x<w; x++)
        {

            int i = x + (h-y-1)*w;
            // calculate direction from eye to scene
            Vector eyeDir = cam.getEyeDir(x,y);
            Vector eyePos = cam.getEyePos(eyeDir);
#if MEDIUM && !VISUAL
            d[i] =volume_radiance(Ray(eyePos,eyeDir.normal()),0,samples);
#elif !VISUAL
            d[i] = surface_radiance(Ray(eyePos,eyeDir.normal()),0,samples); // normal photon mapping
#endif

        }
    }
    transfer_to_pixels(d,pixel);
    glutDisplayFunc(&my_display);

#if VISUAL
    glutIdleFunc(idle);
#endif

    glutMainLoop();
    delete [] d;
     
    
    
}
