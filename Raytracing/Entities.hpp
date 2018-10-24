//
//  entities.hpp
//  Raytracing
//
//  Created by 王登远 on 10/17/18.
//  Copyright © 2018 Dengyuan Wang. All rights reserved.
//

#ifndef Entities_hpp
#define Entities_hpp

#include "Default_setting.h"

class Entities{
    
    //sphere x, y, z, r
    enum obj_name{
        sphere = 0,
        triangle = 1
    }entity_name;
    //material ar, ag, ab, dr, dg, db, sr, sg, sb, ns, tr, tg, tb, por
    float material[14] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 5, 0, 0, 0, 1};
    float entity_parameters[100];//ep: sphere[4] = {0, 0, 2, 1};
    int parameter_num = 0;//is based on name;
public:
    std::string get_entity_name();
    bool set_entity_name(std::string name);
    bool get_property_material(int index,float &value);
    bool set_property_material(int index,float &value);
    bool get_entity_parameter(int index,float &value);
    bool set_entity_parameter(int index,float &value);
    Rays get_normal_vec(float pos[3]);
    bool check_intersect(R_values &R_v,Rays ray);
    Rays get_reflection_ray(Rays ray,float hitpoint[3]);
    Rays get_refraction_ray(Rays ray,float hitpoint[3],float ior_in,float ior_out);
private:
    void calcu_norm_vec(float &rx,float &ry,float &rz,float sx,float sy,float sz,float dx,float dy,float dz){
        rx = dx-sx;ry = dy-sy;rz = dz-sz;
        float l_d = sqrt(rx*rx+ry*ry+rz*rz);
        rx = rx/l_d;ry = ry/l_d;rz = rz/l_d;
    }
    void calcu_cross_vec(float &rx,float &ry,float &rz,float ux,float uy,float uz,float vx,float vy,float vz){
        rx = uy*vz - uz*vy;//a2b3 - a3b2
        ry = uz*vx - ux*vz;//a3b1 - a1b3
        rz = ux*vy - uy*vx;//a1b2 - a2b1
        float l_d = sqrt(rx*rx+ry*ry+rz*rz);
        rx = rx/l_d;ry = ry/l_d;rz = rz/l_d;
    }
    bool triangleIntersection(float hitpoint[3],float &distance ,Rays ray);
};




#endif /* Entities_hpp */
