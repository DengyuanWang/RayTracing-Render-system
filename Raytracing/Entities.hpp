//
//  entities.hpp
//  Raytracing
//
//  Created by 王登远 on 10/17/18.
//  Copyright © 2018 Dengyuan Wang. All rights reserved.
//

#ifndef Entities_hpp
#define Entities_hpp
#include <stdio.h>
#include "Default_setting.h"
class Entities{
    
    //sphere x, y, z, r
    enum obj_name{
        sphere = 0,
        plane = 1
    }obj_name;
    //material ar, ag, ab, dr, dg, db, sr, sg, sb, ns, tr, tg, tb, por
    float material[14] = {0, 0, 0, 1, 1, 1, 0, 0, 0, 5, 0, 0, 0, 1};
    float obj_parameters[100];//ep: sphere[4] = {0, 0, 2, 1};
    int parameter_num = 0;//is based on name;
public:
    std::string get_entity_name();
    std::string set_entity_name();
    bool get_property_material(int index,float value);
    bool set_property_material(int index,float value);
    bool get_entity_parameter(int index,float value);
    bool set_entity_parameter(int index,float value);
    R_values get_normal_vec();
    bool check_intersect(R_values &R_v,Rays ray);
};




#endif /* Entities_hpp */
