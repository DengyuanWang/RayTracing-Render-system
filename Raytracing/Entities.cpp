//
//  entities.cpp
//  Raytracing
//
//  Created by 王登远 on 10/17/18.
//  Copyright © 2018 Dengyuan Wang. All rights reserved.
//

#include "Entities.hpp"

std::string Entities::get_entity_name(){
    switch (entity_name) {
        case obj_name::sphere:
            return "sphere";
        case obj_name::triangle:
            return "triangle";
        default:
            return "";
    }
}
bool Entities::set_entity_name(std::string name){
    if(name=="sphere"){
        entity_name =obj_name::sphere;
        parameter_num = 4;
        return true;
    }else if (name=="triangle"){
        entity_name =obj_name::triangle;
        parameter_num = 30;
        return true;
    }else{
        return false;
    }
}
bool Entities::get_property_material(int index,float &value){
    if(index>=14)
        return false;
    value = material[index];
    return true;
}
bool Entities::set_property_material(int index,float &value){
    if(index>=14)
        return false;
    material[index] = value;
    return true;
}
bool Entities::get_entity_parameter(int index,float &value){
    if(index>=parameter_num)
        return false;
    value = entity_parameters[index];
    return true;
}
bool Entities::set_entity_parameter(int index,float &value){
    if(index>=parameter_num)
        return false;
    entity_parameters[index] = value;
    return true;
}
Rays Entities::get_normal_vec(float pos[3]){
    Rays norm;
    switch (entity_name) {
        case obj_name::sphere:
            calcu_norm_vec(norm.Direction[0],norm.Direction[1],norm.Direction[2],
                           entity_parameters[0],entity_parameters[1],entity_parameters[2]
                                ,pos[0],pos[1],pos[2]);
            return norm;
        case obj_name::triangle:
            return norm;;
        default:
            return norm;
    }
}
Rays Entities::get_reflection_ray(R_values R_v){
    //R = ori_ray - 2*N*dot(ori_ray*N);
    Rays norm = get_normal_vec(R_v.pos);
    Rays ori_ray,reflect_ray;
    calcu_norm_vec(ori_ray.Direction[0],ori_ray.Direction[1],ori_ray.Direction[2],
                    R_v.viewpoint[0],R_v.viewpoint[1],R_v.viewpoint[2],
                            R_v.pos[0],R_v.pos[1],R_v.pos[2]);

    float dot = (ori_ray.Direction[0]*norm.Direction[0]+
                 ori_ray.Direction[1]*norm.Direction[1]+
                 ori_ray.Direction[2]*norm.Direction[2]);
    reflect_ray.Direction[0] = ori_ray.Direction[0]-2*dot*norm.Direction[0];
    reflect_ray.Direction[1] = ori_ray.Direction[1]-2*dot*norm.Direction[1];
    reflect_ray.Direction[2] = ori_ray.Direction[2]-2*dot*norm.Direction[2];
    return reflect_ray;
}


bool Entities::check_intersect(R_values &R_v,Rays ray){
     if(entity_name==obj_name::sphere){
             float x,x0,y,y0,z,z0,dx,dy,dz,r,t1,t2;
             x = ray.Point[0];y = ray.Point[1];z = ray.Point[2];
             dx = ray.Direction[0];dy = ray.Direction[1];dz = ray.Direction[2];
             x0 = entity_parameters[0];
             y0 = entity_parameters[1];
             z0 = entity_parameters[2];
             r  = entity_parameters[3];
             float a,b,c;
             a = dx*dx+dy*dy+dz*dz;
             b = 2.0*(dx*(x-x0)+dy*(y-y0)+dz*(z-z0));
             c = (x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0) - r*r;
             if(b*b-4.0*a*c>=0)//intersect
             {
                 t1 = (-b + sqrt(b*b-4.0*a*c))/(2.0*a);
                 t2 = (-b - sqrt(b*b-4.0*a*c))/(2.0*a);
                 if(fmax(t1,t2)<=0){}// no leagle one
                 else if(ray.range>=0&&fmin(t1,t2)>ray.range){}//check for segment. for ray the range will be -1, mean disable;
                 else if(ray.range>=0){
                     float t = fmin(fmin(ray.range,fmax(0,t1)),fmin(ray.range,fmax(0,t2)));
                     R_v.dis_in_t = t;
                     R_v.pos[0] = x+t*dx;R_v.pos[1] = y+t*dy;R_v.pos[2] = z+t*dz;
                     return true;
                 }else{
                     float t = fmin(fmax(0,t1),fmax(0,t2));
                     R_v.dis_in_t = t;
                     R_v.pos[0] = x+t*dx;R_v.pos[1] = y+t*dy;R_v.pos[2] = z+t*dz;
                     return true;
                 }
             }
         
     }
     return false;
}
