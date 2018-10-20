//
//  Lights.hpp
//  Raytracing
//
//  Created by 王登远 on 10/17/18.
//  Copyright © 2018 Dengyuan Wang. All rights reserved.
//

#ifndef Lights_hpp
#define Lights_hpp

#include "Default_setting.h"

#include "Entities.hpp"
class Lights{
    //directional_light r, g, b, x, y, z
    //float directional_light[6] = {1, 0, 0, 1, 0, 0};
    
    //point_light r, g, b, x, y, z
    //float point_light[6] = {1, 1, 1, 10, 10, 10};
    
    //spot_light r, g, b, px, py, pz, dx, dy, dz, angle1, angle2
    //float spot_light[11] ={1, 1, 1, 5, 5, 5, 1, 0, 0, 60, 120};
    
    //ambient_light r, g, b
    //float ambient_light[3] = {0, 0, 0};
    enum lgt_name{
        directional_light = 0,
        point_light = 1,
        spot_light = 2,
        ambient_light = 3
    }lgt_name;
    
    float lgt_parameters[100];//ep: r,g,b channel range 0~1 to 0~255;
    int parameter_num = 0;//is based on name;
public:
    std::string get_lgt_name();
    bool set_lgt_name(std::string name);
    bool get_lgt_parameter(int index,float value);
    bool set_lgt_parameter(int index,float value);
    bool calcu_color(float RGB[3],Entities *entity_it,float hit_point[3],float viewpoint[3]);
    Rays get_a_ray(float hit_point[3]);
private:
    void Lambertian_shading(float RBG[3],Entities *entity_it,float hit_point[3]);
    void specularity_shading_phong(float RGB[3], Entities *entity_it,float hit_point[3],float viewpoint[3]);
    void calcu_norm_vec(float &rx,float &ry,float &rz,float sx,float sy,float sz,float dx,float dy,float dz){
        rx = dx-sx;ry = dy-sy;rz = dz-sz;
        float l_d = sqrt(rx*rx+ry*ry+rz*rz);
        rx = rx/l_d;ry = ry/l_d;rz = rz/l_d;
        
    }
};
#endif /* Lights_hpp */
