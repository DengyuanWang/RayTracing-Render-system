//
//  Lights.hpp
//  Raytracing
//
//  Created by 王登远 on 10/17/18.
//  Copyright © 2018 Dengyuan Wang. All rights reserved.
//

#ifndef Lights_hpp
#define Lights_hpp

#include <stdio.h>
#include "Default_setting.h"
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
    std::string set_lgt_name();
    bool get_property_material(int index,float value);
    bool set_property_material(int index,float value);
    bool get_lgt_parameter(int index,float value);
    bool set_lgt_parameter(int index,float value);
    R_values get_normal_vec();
    bool check_intersect(R_values &R_v,Rays ray);
};
#endif /* Lights_hpp */
