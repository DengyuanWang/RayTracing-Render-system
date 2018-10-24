//
//  Skybox.hpp
//  Raytracing
//
//  Created by 王登远 on 10/24/18.
//  Copyright © 2018 Dengyuan Wang. All rights reserved.
//

#ifndef Skybox_hpp
#define Skybox_hpp

#include <stdio.h>
#include "image.h"
#include "Lights.hpp"
class Skybox{
    
    float sky_size[2]={20,20};
public:
    std::list <Lights> Lgts;
    bool load_skybox(std::string filename);
    bool set_sky_size(float x,float y);
    
private:
    bool add_point_light(float r,float g,float b,float x,float y,float z){
        //point_light r, g, b, x, y, z
        //float point_light[6] = {255, 255, 255, 10, 10, 10};
        Lights tmp;
        tmp.set_lgt_name("point_light");
        tmp.set_lgt_parameter(0, r);tmp.set_lgt_parameter(1, g);tmp.set_lgt_parameter(2, b);
        tmp.set_lgt_parameter(3, x);tmp.set_lgt_parameter(4, y);tmp.set_lgt_parameter(5, z);
        Lgts.push_front(tmp);
        return true;
    }
};
#endif /* Skybox_hpp */
