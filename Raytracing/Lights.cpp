//
//  Lights.cpp
//  Raytracing
//
//  Created by 王登远 on 10/17/18.
//  Copyright © 2018 Dengyuan Wang. All rights reserved.
//

#include "Lights.hpp"
#include "Entities.hpp"
std::string Lights::get_lgt_name(){
    switch (lgt_name) {
        case spot_light:
            return "spot_light";
        case point_light:
            return "point_light";
        case ambient_light:
            return "ambient_light";
        case directional_light:
            return "directional_light";
        default:
            break;
    }
}
bool Lights::set_lgt_name(std::string name){
    if(name=="spot_light"){
        lgt_name = spot_light;
        parameter_num = 11;
    }else if (name=="point_light"){
        lgt_name = point_light;
        parameter_num = 6;
    }else if(name == "ambient_light"){
        lgt_name = ambient_light;
        parameter_num = 3;
    }else if (name == "directional_light"){
        lgt_name = directional_light;
        parameter_num = 6;
    }else{
        return false;
    }
    return true;
}
bool Lights::get_lgt_parameter(int index,float value){
    if(index>=parameter_num)
        return false;
    value = lgt_parameters[index];
    return true;
}
bool Lights::set_lgt_parameter(int index,float value){
    if(index>=parameter_num)
        return false;
    lgt_parameters[index] = value;
    return true;
}
bool Lights::calcu_color(float &R,float &G,float &B,R_values R_v){
    R_v.rgb[0] = 0;
    R_v.rgb[1] = 0;
    R_v.rgb[2] = 0;
    bool delum_tag=true;
    Entities *it = (Entities*)R_v.entity;
    switch (lgt_name) {
        case lgt_name::directional_light:
            delum_tag=false;
            break;
        case lgt_name::ambient_light:
            //printf("ambient r = %f,g = %f,b = %f\n",R,G,B);
            float ar,ag,ab;
            it->get_property_material(0, ar);
            it->get_property_material(1, ag);
            it->get_property_material(2, ab);
            //printf("ar = %f,ag = %f,ab = %f\n",ar,ag,ab);
            R = ar*lgt_parameters[0];
            G = ag*lgt_parameters[1];
            B = ab*lgt_parameters[2];
            return true;
            break;
        default:
            break;
    }
    //printf("others r = %f,g = %f,b = %f\n",R,G,B);
    if(lgt_name==lgt_name::directional_light)
         delum_tag=false;
    else delum_tag=true;
    Lambertian_shading(R_v);//add rgb to R_v
    specularity_shading_phong(R_v);//add rgb to R_v; error
    R = R_v.rgb[0];
    G = R_v.rgb[1];
    B = R_v.rgb[2];
    return true;
}
Rays Lights::get_a_ray(R_values R_v){//from a point to light
    Rays ray;
    float x = lgt_parameters[3],y = lgt_parameters[4],z = lgt_parameters[5];
    switch (lgt_name) {
        case lgt_name::directional_light:
            // normal direction
            // this ray is from light to point
            calcu_norm_vec(ray.Direction[0],ray.Direction[1],ray.Direction[2],
                               2*x,2*y,2*z,x,y,z);
            ray.Point[0] = R_v.pos[0]+0.001*ray.Direction[0];
            ray.Point[1] = R_v.pos[1]+0.001*ray.Direction[1];
            ray.Point[2] = R_v.pos[2]+0.001*ray.Direction[2];
            ray.range = -1;
            return ray;
        case lgt_name::ambient_light:
            ray.range = -2;
            return ray;
        default:
            // normal direction
            // this ray is from light to point

            calcu_norm_vec(ray.Direction[0],ray.Direction[1],ray.Direction[2],
                               R_v.pos[0],R_v.pos[1],R_v.pos[2],x,y,z);
            ray.Point[0] = R_v.pos[0]+0.001*ray.Direction[0];
            ray.Point[1] = R_v.pos[1]+0.001*ray.Direction[1];
            ray.Point[2] = R_v.pos[2]+0.001*ray.Direction[2];
            ray.range = (x-R_v.pos[0])/ray.Direction[0];
            return ray;
    }
    
}
void Lights::Lambertian_shading(R_values &R_v){
    bool delum_tag=true;
    if(lgt_name==lgt_name::directional_light)
        delum_tag = false;
    float dr,dg,db;
    Entities *it = (Entities*)R_v.entity;
    it->get_property_material(3, dr);it->get_property_material(4, dg);it->get_property_material(5, db);
    float n[3];
    Rays OP;
    OP = it->get_normal_vec(R_v.pos);
    n[0] = OP.Direction[0];n[1] = OP.Direction[1];n[2] = OP.Direction[2];
    Rays ray = get_a_ray(R_v);//get a ray from hit point to light
    float delta_r,delta_g,delta_b,theta;
    theta = 180*acos(n[0]*ray.Direction[0]+n[1]*ray.Direction[1]+n[2]*ray.Direction[2])/M_PI;
    //printf("theta = %f\n",theta);
    //printf("dr=%f dg=%f db=%f; lr=%f lg=%f lb=%f\n",dr,dg,db,L_header->lgt_parameters[0]
    // ,L_header->lgt_parameters[1],L_header->lgt_parameters[2]);
    float delum_factor=1;
    if(delum_tag){
        float x,y,z;
        x = lgt_parameters[3];y = lgt_parameters[4];z = lgt_parameters[5];
        delum_factor = delum_factor/(pow(R_v.pos[0]-x,2)+pow(R_v.pos[1]-y,2)+pow(R_v.pos[2]-z,2));
    }
    delta_r = dr*delum_factor*lgt_parameters[0]*fmax(0,cos(M_PI*theta/180));
    delta_g = dg*delum_factor*lgt_parameters[1]*fmax(0,cos(M_PI*theta/180));
    delta_b = db*delum_factor*lgt_parameters[2]*fmax(0,cos(M_PI*theta/180));
    R_v.rgb[0] += delta_r;
    R_v.rgb[1] += delta_g;
    R_v.rgb[2] += delta_b;
}
void Lights::specularity_shading_phong(R_values &R_v){
    float sr,sg,sb;
    Entities *it = (Entities*)R_v.entity;
    it->get_property_material(6, sr);it->get_property_material(7, sg);it->get_property_material(8, sb);
    bool delum_tag=true;
    if(lgt_name==lgt_name::directional_light)
        delum_tag = false;
    //phong shading specularity
    //Is_r = sr*I_r*dot(v,R)^n_
    //R = I - 2*N*dot(I*N)
    float v[3],n[3],R_[3];//view,normal,reflection
    float n_;//power factor
    it->get_property_material(9, n_);//get power factor
    Rays norm;
    norm = it->get_normal_vec(R_v.pos);//got norm
    n[0] = norm.Direction[0];n[1] = norm.Direction[1];n[2] = norm.Direction[2];
    
    Rays ray=get_a_ray(R_v);//from point to light
    ray.Direction[0] = -ray.Direction[0];ray.Direction[1] = -ray.Direction[1];ray.Direction[2] = -ray.Direction[2];//flip direction
    calcu_norm_vec(v[0],v[1],v[2],
                   R_v.pos[0],R_v.pos[1],R_v.pos[2]
                   ,R_v.viewpoint[0],R_v.viewpoint[1],R_v.viewpoint[2]);//a vector from pos to viewpoint
    float dot = (ray.Direction[0]*n[0]+ray.Direction[1]*n[1]+ray.Direction[2]*n[2]);
    R_[0] = ray.Direction[0]-2*dot*n[0];
    R_[1] = ray.Direction[1]-2*dot*n[1];
    R_[2] = ray.Direction[2]-2*dot*n[2];
    
    dot = pow(v[0]*R_[0]+v[1]*R_[1]+v[2]*R_[2],n_);

    
    float delum_factor=1;
    if(delum_tag){
        float x,y,z;
        x = lgt_parameters[3];y = lgt_parameters[4];z = lgt_parameters[5];
        delum_factor = delum_factor/(pow(R_v.pos[0]-x,2)+pow(R_v.pos[1]-y,2)+pow(R_v.pos[2]-z,2));
    }
    R_v.rgb[0] +=sr*delum_factor*lgt_parameters[0]*dot;
    R_v.rgb[1] +=sg*delum_factor*lgt_parameters[1]*dot;
    R_v.rgb[2] +=sb*delum_factor*lgt_parameters[2]*dot;
}
