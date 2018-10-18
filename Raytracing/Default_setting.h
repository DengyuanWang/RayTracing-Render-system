//
//  Default_setting.h
//  Raytracing
//
//  Created by 王登远 on 10/17/18.
//  Copyright © 2018 Dengyuan Wang. All rights reserved.
//

#ifndef Default_setting_h
#define Default_setting_h
#define NULL nullptr

#include <list>
typedef struct{
    //setting input filename
    string setting_file = " ";
    //output_image filename
    string output_image = "retraced.bmp";
    //camera px, py, pz, dx, dy, dz, ux, uy, uz, ha
    float camera[10] = {0, 0, 0, 0, 0, 1, 0, 1, 0, 45};
    //film_resolution width, height
    int film_resolution[2] = {480, 640};//width and height fliped
    //max_depth n
    float max_depth = 5;
    //background r, g, b
    float background[3] = {0, 0, 0};
    //material ar, ag, ab, dr, dg, db, sr, sg, sb, ns, tr, tg, tb, por
    float material[14] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 5, 0, 0, 0, 1};
    
}World_Setting;
typedef struct W_Objects{
    //sphere x, y, z, r
    enum obj_name{
        sphere = 0,
        plane = 1
    }obj_name;
    //material ar, ag, ab, dr, dg, db, sr, sg, sb, ns, tr, tg, tb, por
    float material[14] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 5, 0, 0, 0, 1};
    float *obj_parameters = nullptr;//ep: sphere[4] = {0, 0, 2, 1};
    W_Objects *next = nullptr;
}W_Objects,*W_Objects_pt;

typedef struct Rays{
    float Point[3];
    float Direction[3];
    float range = -1;
}Rays,*Ray_p;
typedef struct Screen{
    float position[3];
    float rgb[3];//range 0~255
}Screen,*Screen_p;
typedef struct R_values{
    W_Objects_pt obj;
    Lights lgt;
    float pos[3];
    float dir[3];
    float rgb[3];// range 0~1
}R_values;

#endif /* Default_setting_h */
