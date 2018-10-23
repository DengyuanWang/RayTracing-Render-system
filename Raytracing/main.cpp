//
//  main.cpp
//  Raytracing
//
//  Created by 王登远 on 10/8/18.
//  Copyright © 2018 Dengyuan Wang. All rights reserved.
//


#include "Default_setting.h"
#include "Interface.hpp"
#include "World.hpp"

#define STB_IMAGE_IMPLEMENTATION //only place once in one .cpp file
#include "stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION //only place once in one .cpp files
#include "stb_image_write.h"



std::string gene_img(std::string );

int main(int argv,char *argc[]){
    std::string out_filename;
    std::string setting_filename = "directional_light.scn";
    out_filename = gene_img(setting_filename);
    
    printf("genereate img finished\n");
/*
    Interface interface;
    interface.create_window(640,480);
    // GUI
    
    interface.load_img(out_filename);

     
     
    while(true){//display loop
        if(interface.update_window()==true)//true for update_end
             break;
    }
    interface.close_window();
*/
    return EXIT_SUCCESS;
}
std::string gene_img(std::string setting_filename){
    int sx=3,sy=3;
    std::string rst = "";
    World world(setting_filename);
    printf("load setting of world begin\n");
    world.load_setting();
    printf("load setting of world finished\n");
    World_Setting W_Setting = world.get_worldSetting();
    
    //enlarge image for supersampling:
    world.Enlarge_resolution(sx, sy);
    printf("Generate ray begin\n");
    world.Generate_ray();
    printf("Generate ray finished\n");
    world.ray_tracing();
    //write world screen into bmp file;
    world.generate_img(sx,sy);
    world.save_img();
    rst = world.get_imagename();

    return rst;
}


/* Copyright (c) Mark J. Kilgard, 1994. */

/* This program is freely distributable without licensing fees
 and is provided without guarantee or warrantee expressed or
 implied. This program is -not- in the public domain. */
/*
#include <string.h>
#include <GLUT/glut.h>

void *font = GLUT_BITMAP_TIMES_ROMAN_24;
void *fonts[] =
{
    GLUT_BITMAP_9_BY_15,
    GLUT_BITMAP_TIMES_ROMAN_10,
    GLUT_BITMAP_TIMES_ROMAN_24
};
char defaultMessage[] = "GLUT means OpenGL.";
char *message = defaultMessage;

void
selectFont(int newfont)
{
    font = fonts[newfont];
    glutPostRedisplay();
}

void
selectMessage(int msg)
{
    switch (msg) {
        case 1:
            message = "abcdefghijklmnop";
            break;
        case 2:
            message = "ABCDEFGHIJKLMNOP";
            break;
    }
}

void
selectColor(int color)
{
    switch (color) {
        case 1:
            glColor3f(0.0, 1.0, 0.0);
            break;
        case 2:
            glColor3f(1.0, 0.0, 0.0);
            break;
        case 3:
            glColor3f(1.0, 1.0, 1.0);
            break;
    }
    glutPostRedisplay();
}

void
tick(void)
{
    glutPostRedisplay();
}

void
output(int x, int y, char *string)
{
    int len, i;
    
    glRasterPos2f(x, y);
    len = (int) strlen(string);
    for (i = 0; i < len; i++) {
        glutBitmapCharacter(font, string[i]);
    }
}

void
display(void)
{
    glClear(GL_COLOR_BUFFER_BIT);
    output(0, 24, "This is written in a GLUT bitmap font.");
    output(100, 100, message);
    output(50, 145, "(positioned in pixels with upper-left origin)");
    glutSwapBuffers();
}

void
reshape(int w, int h)
{
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, w, h, 0);
    glMatrixMode(GL_MODELVIEW);
}

int main(int argc, char **argv)
{
   // std::string setting_filename = "spheres1.scn";
   // gene_img(setting_filename);
    
    
    int i, msg_submenu, color_submenu;
    
    glutInit(&argc, argv);
    for (i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-mono")) {
            font = GLUT_BITMAP_9_BY_15;
        }
    }
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(500, 150);
    glutCreateWindow("Sphere");
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutIdleFunc(tick);
    msg_submenu = glutCreateMenu(selectMessage);
    glutAddMenuEntry("abc", 1);
    glutAddMenuEntry("ABC", 2);
    color_submenu = glutCreateMenu(selectColor);
    glutAddMenuEntry("Green", 1);
    glutAddMenuEntry("Red", 2);
    glutAddMenuEntry("White", 3);
    glutCreateMenu(selectFont);
    glutAddMenuEntry("9 by 15", 0);
    glutAddMenuEntry("Times Roman 10", 1);
    glutAddMenuEntry("Times Roman 24", 2);
    glutAddSubMenu("Messages", msg_submenu);
    glutAddSubMenu("Color", color_submenu);
    glutAttachMenu(GLUT_RIGHT_BUTTON);
    glutMainLoop();
    return 0;
}
*/
