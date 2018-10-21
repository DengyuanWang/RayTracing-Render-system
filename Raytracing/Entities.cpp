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
Rays Entities::get_reflection_ray(Rays ray,float hitpoint[3]){//intput ray, which hit at hitpoint
    //R = ori_ray - 2*N*dot(ori_ray*N);
    Rays norm = get_normal_vec(hitpoint);
    Rays reflect_ray;
    //dot(inputray, norm)
    float dot = (ray.Direction[0]*norm.Direction[0]+
                 ray.Direction[1]*norm.Direction[1]+
                 ray.Direction[2]*norm.Direction[2]);
    reflect_ray.Direction[0] = ray.Direction[0]-2*dot*norm.Direction[0];
    reflect_ray.Direction[1] = ray.Direction[1]-2*dot*norm.Direction[1];
    reflect_ray.Direction[2] = ray.Direction[2]-2*dot*norm.Direction[2];
    reflect_ray.Point[0] = hitpoint[0]+0.001*reflect_ray.Direction[0];
    reflect_ray.Point[1] = hitpoint[1]+0.001*reflect_ray.Direction[1];
    reflect_ray.Point[2] = hitpoint[2]+0.001*reflect_ray.Direction[2];
    return reflect_ray;
}
Rays Entities::get_refraction_ray(Rays input_ray,float hitpoint[3],float ior_in,float ior_out){//intput ray, which hit at hitpoint
    //L2 = 1*sin(theta_refra)/sin(theta_input)
    //refraction = -norm*|1*cos(theta_refra)| + input_ray*(L2/|input_ray|) + norm*|L2*cos(theta_input)|
    Rays norm = get_normal_vec(hitpoint);
    Rays refraction;
    float theta_refra,theta_input;
    theta_input = acos(norm.Direction[0]*input_ray.Direction[0]+
                       norm.Direction[1]*input_ray.Direction[1]+
                       norm.Direction[2]*input_ray.Direction[2]);
    if(theta_input>M_PI/2)//means from out to in
    {
        //printf(">90\n");
        theta_input = acos( -norm.Direction[0]*input_ray.Direction[0]+
                            -norm.Direction[1]*input_ray.Direction[1]+
                            -norm.Direction[2]*input_ray.Direction[2]);
    }else if(theta_input<M_PI/2)//means from in to out
    {
        //printf("<90\n");
        norm.Direction[0] = -norm.Direction[0];
        norm.Direction[1] = -norm.Direction[1];
        norm.Direction[2] = -norm.Direction[2];
    }else{
        float len1 = pow(input_ray.Direction[0],2)+pow(input_ray.Direction[1],2)+pow(input_ray.Direction[2],2);
        float len2 = pow(norm.Direction[0],2)+pow(norm.Direction[1],2)+pow(norm.Direction[2],2);
        printf("error in refraction should not shoot it\n");
        exit(0);
    }
    
    float temp = sin(theta_input)*ior_in/ior_out;
    if(temp<1&&temp>-1)
        theta_refra = asin(sin(theta_input)*ior_in/ior_out);
    else{
        //printf("error in calcu refraction ray%f\n",abs(sin(theta_input)*ior_in/ior_out));
        refraction.range = -2;
        return refraction;
    }
    float L2 =1*sin(theta_refra)/sin(theta_input);
    float len =pow(input_ray.Direction[0],2)+pow(input_ray.Direction[1],2)+pow(input_ray.Direction[2],2);
    if(len-1>0.0001){
        printf("error in calcu refraction ray about input_ray length\n");
        exit(0);
    }
    refraction.Direction[0] = -norm.Direction[0]*         1*cos(theta_refra) +
                                input_ray.Direction[0]*   L2 +
                                norm.Direction[0]*        L2*cos(theta_input);
    refraction.Direction[1] = -norm.Direction[1]*         1*cos(theta_refra) +
                                input_ray.Direction[1]*   L2 +
                                norm.Direction[1]*        L2*cos(theta_input);
    refraction.Direction[2] = -norm.Direction[2]*         1*cos(theta_refra) +
                                input_ray.Direction[2]*   L2 +
                                norm.Direction[2]*        L2*cos(theta_input);
    float len_ray = pow(refraction.Direction[0],2)+pow(refraction.Direction[1],2)+pow(refraction.Direction[2],2);
    refraction.Direction[0]/=len_ray;refraction.Direction[1]/=len_ray;refraction.Direction[2]/=len_ray;
    refraction.Point[0] = hitpoint[0]+0.01*refraction.Direction[0];
    refraction.Point[1] = hitpoint[1]+0.01*refraction.Direction[1];
    refraction.Point[2] = hitpoint[2]+0.01*refraction.Direction[2];
    return refraction;
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
                 t1 = t1<0.001?0:t1;
                 t2 = t2<0.001?0:t2;
                 if(fmax(t1,t2)<=0){}// no leagle one
                 else if(t1>0&&t2>0){
                    if(ray.range>=0&&fmin(t1,t2)>ray.range){}//check for segment. for ray with range -1, means infinity;
                    else{
                         float t = fmin(t1,t2);//select nearer one
                         R_v.dis_in_t = t;
                         R_v.pos[0] = x+t*dx;R_v.pos[1] = y+t*dy;R_v.pos[2] = z+t*dz;
                         if(abs(t)<0.0001){
                             printf("error in check intersection\n");
                             exit(0);
                         }
                         return true;
                     }
                 }else if(t1>0){
                     if(ray.range>=0&&t1>ray.range){}//check for segment. for ray with range -1, means infinity;
                     else{
                         float t = t1;
                         R_v.dis_in_t = t;
                         R_v.pos[0] = x+t*dx;R_v.pos[1] = y+t*dy;R_v.pos[2] = z+t*dz;
                         if(abs(t)<0.001){
                             printf("error in check intersection\n");
                             exit(0);
                         }
                         return true;
                     }
                 }else{
                     if(ray.range>=0&&t2>ray.range){}//check for segment. for ray with range -1, means infinity;
                     else{
                         float t = t2;
                         R_v.dis_in_t = t;
                         R_v.pos[0] = x+t*dx;R_v.pos[1] = y+t*dy;R_v.pos[2] = z+t*dz;
                         if(abs(t)<0.001){
                             printf("error in check intersection\n");
                             exit(0);
                         }
                         return true;
                     }
                 }
             }
         
     }
     return false;
}
