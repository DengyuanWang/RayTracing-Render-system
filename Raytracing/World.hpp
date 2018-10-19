//
//  World.hpp
//  Raytracing
//
//  Created by 王登远 on 10/9/18.
//  Copyright © 2018 Dengyuan Wang. All rights reserved.
//

#ifndef World_hpp
#define World_hpp



#include "image.h"
#include "pixel.h"
#include "Default_setting.h"
#include "Entities.hpp"
#include "Lights.hpp"


class World {
    Image *img=nullptr;
    World_Setting W_Settings;
    Ray_p Ray=nullptr;
    Screen_p Scn=nullptr;
    std::list <Entities> Ets;
    std::list <Lights> Lgts;
public:
    World(std::string);
//add wall and ground to the world;
    bool build_world();
    bool add_plane(float x,float y,float z,float dx,float dy,float dz);
//load setting and generate world
    bool load_setting();
    bool add_sphere(float x,float y,float z,float r);
    bool add_directional_light(float r,float g,float b,float x,float y,float z);
    bool add_point_light(float r,float g,float b,float x,float y,float z);
    bool add_spot_light(float r,float g,float b,float px,float py,float pz,float dx,float dy,float dz,float angle1,float angle2);
    bool add_ambient_light(float r,float g,float b);
//generate ray
    bool Generate_ray();
    int index_ray(int x_index,int y_index);
//ray tracing
    bool ray_tracing();
    bool check_intersect(R_values &R,Rays ray);
    bool calcu_color(R_values &R,Rays ray,int current_recursive_depth);
    void calcu_rgb_under_lights(float &R,float &G,float &B,R_values R_v,int lgt_index,bool delum_tag);
    void calcu_norm_vec(float &rx,float &ry,float &rz,float sx,float sy,float sz,float dx,float dy,float dz);
    
    bool generate_img(float sx,float sy);
    std::string get_imagename();
    bool Enlarge_resolution(float sx,float sy);
    bool save_img();
    World_Setting get_worldSetting();
};

#endif /* World_hpp */
