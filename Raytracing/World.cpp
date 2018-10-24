//
//  World.cpp
//  Raytracing
//
//  Created by 王登远 on 10/9/18.
//  Copyright © 2018 Dengyuan Wang. All rights reserved.
//


#include "World.hpp"

//World class

World::World(std::string filename){
    
    W_Settings.setting_file = filename;//specify the setting file
    
}
bool World::add_skybox(){
    Environment.load_skybox("stpeters_cross.bmp");
    
    return true;
}
bool World::load_setting(){
    std::string line;
    
    std::string fileName = W_Settings.setting_file;
    
    // open the file containing the scene description
    std::ifstream input(fileName);
    
    // check for errors in opening the file
    if(input.fail()){
        std::cout << "Can't open file '" << fileName << "'" << std::endl;
        return false;
    }
    
    // determine the file size (this is optional -- feel free to delete the 6 lines below)
    std::streampos begin,end;
    begin = input.tellg();
    input.seekg(0, std::ios::end);
    end = input.tellg();
    std::cout << "File '" << fileName << "' is: " << (end-begin) << " bytes long.\n\n";
    input.seekg(0, std::ios::beg);
    
    
    //Loop through reading each line
    std::string command;
    while(input >> command) { //Read first word in the line (i.e., the command type)
        
        if (command[0] == '#'){
            getline(input, line); //skip rest of line
            std::cout << "Skipping comment: " << command  << line <<  std::endl;
            continue;
        }
        
        
        if (command == "sphere"){ //If the command is a sphere command
            float x,y,z,r;
            input >> x >> y >> z >> r;
            if(add_sphere(x, y, z, r)==true)
                printf("Sphere as position (%f,%f,%f) with radius %f\n",x,y,z,r);
            else printf("Sphere Add wrong as position (%f,%f,%f) with radius %f\n",x,y,z,r);
        }
        else if (command == "background"){ //If the command is a background command
            float r,g,b;
            input >> r >> g >> b;
            W_Settings.background[0] = r;
            W_Settings.background[1] = g;
            W_Settings.background[2] = b;
            printf("Background color of (%f,%f,%f)\n",r,g,b);
        }
        else if (command == "max_vertices"){ //If the command is a max_vertices command
            float max_vertices;
            input >> max_vertices;
            W_Settings.max_vertices= max_vertices;
            printf("Max_vertices of %f\n",max_vertices);
        }
        else if (command == "vertex"){ //If the command is a vertex command
            vertex_ vertex;
            input >> vertex.x>> vertex.y>> vertex.z;
            //push the vertex to the tail;
            if(Vertices_pool.size()<W_Settings.max_vertices){
                Vertices_pool.push_back(vertex);
                printf("vertex of (%f,%f,%f)\n",vertex.x,vertex.y,vertex.z);
            }
            else{
                printf("Vertices_pool Overflowed\n");
                exit(0);
            }
            
        }
        else if (command == "triangle"){
            int index1,index2,index3;
            input >> index1>> index2>>index3;
            add_triangle(index1, index2, index3);
            printf("add triangle %d %d %d\n",index1,index2,index3);
        }
        else if (command == "camera"){ //If the command is a camera command
            float px, py, pz, dx, dy, dz, ux, uy, uz, ha;
            input >>px >> py>> pz>> dx>> dy>> dz>> ux>> uy>> uz>> ha;
            float l_d = sqrt(dx*dx+dy*dy+dz*dz);//normalize
            float l_u = sqrt(ux*ux+uy*uy+uz*uz);//normalize
            W_Settings.camera[0] = px;W_Settings.camera[1] = py;W_Settings.camera[2] = pz;
            W_Settings.camera[3] = dx/l_d;W_Settings.camera[4] = dy/l_d;W_Settings.camera[5] = dz/l_d;
            W_Settings.camera[6] = ux/l_u;W_Settings.camera[7] = uy/l_u;W_Settings.camera[8] = uz/l_u;
            W_Settings.camera[9] = ha;
            printf("camera parameters of:\n(%f,%f,%f)\n(%f,%f,%f)\n(%f,%f,%f)\n(%f)\n",px, py, pz, dx, dy, dz, ux, uy, uz, ha);
        }
        else if (command == "film_resolution"){ //If the command is a film_resolution command
            float width, height;
            input >> width>> height;
            W_Settings.film_resolution[0] = width;
            W_Settings.film_resolution[1] = height;
            printf("file_resolution of (%f,%f)\n",width,height);
        }
        else if (command == "material"){ //If the command is a material command
            float ar, ag, ab, dr, dg, db, sr, sg, sb, ns, tr, tg, tb, por;
            input >> ar>> ag>> ab>> dr>> dg>> db>> sr>> sg>> sb>> ns>> tr>> tg>> tb>> por;
            W_Settings.material[0] = ar;W_Settings.material[1] = ag;W_Settings.material[2] = ab;
            W_Settings.material[3] = dr;W_Settings.material[4] = dg;W_Settings.material[5] = db;
            W_Settings.material[6] = sr;W_Settings.material[7] = sg;W_Settings.material[8] = sb;
            W_Settings.material[9] = ns;
            W_Settings.material[10] = tr;W_Settings.material[11] = tg;W_Settings.material[12] = tb;
            W_Settings.material[13] = por;
            printf("material parameters of:\n(%f,%f,%f)\n(%f,%f,%f)\n(%f,%f,%f)\n(%f,%f,%f)\n(%f,%f)\n",ar, ag, ab, dr, dg, db, sr, sg, sb, ns, tr, tg, tb, por);
        }
        else if (command == "directional_light"){ //If the command is a directional_light command
            float r, g, b, x, y, z;
            input >> r>> g>> b>> x>> y>> z;
            if(add_directional_light(r,g,b,x,y,z)==true)
                printf("directional_light parameters of:\n(%f,%f,%f)\n(%f,%f,%f)\n",r, g, b, x, y, z);
            else printf("directional_light add error parameters of:\n(%f,%f,%f)\n(%f,%f,%f)\n",r, g, b, x, y, z);
        }
        else if (command == "point_light"){ //If the command is a point_light command
            float r, g, b, x, y, z;
            input >> r>> g>> b>> x>> y>> z;
            if(add_point_light(r,g,b,x,y,z)==true)
                printf("point_light parameters of:\n(%f,%f,%f)\n(%f,%f,%f)\n",r, g, b, x, y, z);
            else printf("point_light add error parameters of:\n(%f,%f,%f)\n(%f,%f,%f)\n",r, g, b, x, y, z);
        }
        else if (command == "spot_light"){ //If the command is a spot_light command
            float r, g, b, px, py, pz, dx, dy, dz, angle1, angle2;
            input >> r>> g>> b>> px>> py>> pz>> dx>> dy>> dz>> angle1>> angle2;
            if(add_spot_light(r, g, b, px, py, pz, dx, dy, dz, angle1, angle2)==true)
                printf("spot_light parameters of:\n(%f,%f,%f)\n(%f,%f,%f)\n(%f,%f,%f)\n(%f,%f)\n",r, g, b, px, py, pz, dx, dy, dz, angle1, angle2);
            else printf("spot_light add error parameters of:\n(%f,%f,%f)\n(%f,%f,%f)\n(%f,%f,%f)\n(%f,%f)\n",r, g, b, px, py, pz, dx, dy, dz, angle1, angle2);
        }
        else if (command == "ambient_light"){ //If the command is a ambient_light command
            float r, g, b;
            input >> r>> g>> b;
            if(add_ambient_light(r,g,b)==true)
                printf("ambient_light parameters of:\n(%f,%f,%f)\n",r, g, b);
            else printf("ambient_light add error parameters of:\n(%f,%f,%f)\n",r, g, b);
        }
        else if (command == "output_image"){ //If the command is an output_image command
            std::string outFile;
            input >> outFile;
            W_Settings.output_image = outFile;
            printf("Render to file named: %s\n", outFile.c_str());
        }
        else if (command == "max_depth"){ //If the command is an max_depth command
            float max_depth;
            input >> max_depth;
            W_Settings.max_depth = max_depth;
            printf("max_depth of: %f\n", max_depth);
        }
        else if (command == "skybox"){ //If the command is an max_depth command
            W_Settings.skybox_tag = true;
            add_skybox();
            printf("Sky box true\n");
        }
        else {
            getline(input, line); //skip rest of line
            std::cout << "WARNING. Do not know command: " << command << std::endl;
        }
    }
    if(W_Settings.skybox_tag)
    {
        int i=0;
        for(std::list<Lights>::iterator it=Environment.Lgts.begin();it!=Environment.Lgts.end();++it)
        {
            Lights t = *it;
            Lgts.push_front(t);
            ++i;
        }
        printf("i = %d\n",i);
    }
    return true;
}
bool World::add_plane(float x,float y,float z,float dx,float dy,float dz){//one point and one n could decide a plane
    
    return false;
}
bool World::add_sphere(float x,float y,float z,float r){
    Entities tmp;
    tmp.set_entity_name("sphere");
    tmp.set_entity_parameter(0, x);
    tmp.set_entity_parameter(1, y);
    tmp.set_entity_parameter(2, z);
    tmp.set_entity_parameter(3, r);
    for(int i=0;i<14;i++)
        tmp.set_property_material(i,  W_Settings.material[i]);
    Ets.push_front(tmp);
    return true;
}
bool World::add_triangle(int index1,int index2,int index3){
    float vertex1[3],vertex2[3],vertex3[3];
    int i=1,j=0;;
    for(std::list<vertex_>::iterator it=Vertices_pool.begin();it!=Vertices_pool.end();++it,++i){
        if(i==index1){
            vertex1[0] = it->x;vertex1[1] = it->y;vertex1[2] = it->z;++j;
        }else if(i==index2){
            vertex2[0] = it->x;vertex2[1] = it->y;vertex2[2] = it->z;++j;
        }else if (i==index3){
            vertex3[0] = it->x;vertex3[1] = it->y;vertex3[2] = it->z;++j;
        }else{
            
        }
        if(j==3)
            break;
    }
    Entities tmp;
    tmp.set_entity_name("triangle");
    tmp.set_entity_parameter(0, vertex1[0]);tmp.set_entity_parameter(1, vertex1[1]);tmp.set_entity_parameter(2, vertex1[2]);
    tmp.set_entity_parameter(3, vertex2[0]);tmp.set_entity_parameter(4, vertex2[1]);tmp.set_entity_parameter(5, vertex2[2]);
    tmp.set_entity_parameter(6, vertex3[0]);tmp.set_entity_parameter(7, vertex3[1]);tmp.set_entity_parameter(8, vertex3[2]);
    for(int i=0;i<14;i++)
        tmp.set_property_material(i,  W_Settings.material[i]);
    Ets.push_front(tmp);
    return true;
}
bool World::add_directional_light(float r,float g,float b,float x,float y,float z){
    //point_light r, g, b, x, y, z
    //float point_light[6] = {255, 255, 255, 10, 10, 10};
    Lights tmp;
    tmp.set_lgt_name("directional_light");
    tmp.set_lgt_parameter(0, r);tmp.set_lgt_parameter(1, g);tmp.set_lgt_parameter(2, b);
    tmp.set_lgt_parameter(3, x);tmp.set_lgt_parameter(4, y);tmp.set_lgt_parameter(5, z);
    Lgts.push_front(tmp);
    return true;
}
bool World::add_point_light(float r,float g,float b,float x,float y,float z){
    //point_light r, g, b, x, y, z
    //float point_light[6] = {255, 255, 255, 10, 10, 10};
    Lights tmp;
    tmp.set_lgt_name("point_light");
    tmp.set_lgt_parameter(0, r);tmp.set_lgt_parameter(1, g);tmp.set_lgt_parameter(2, b);
    tmp.set_lgt_parameter(3, x);tmp.set_lgt_parameter(4, y);tmp.set_lgt_parameter(5, z);
    Lgts.push_front(tmp);
    return true;
}
bool World::add_spot_light(float r,float g,float b,float px,float py,float pz,float dx,float dy,float dz,float angle1,float angle2){
    //point_light r, g, b, x, y, z
    //float point_light[6] = {255, 255, 255, 10, 10, 10};
    Lights tmp;
    tmp.set_lgt_name("spot_light");
    tmp.set_lgt_parameter(0, r);tmp.set_lgt_parameter(1, g);tmp.set_lgt_parameter(2, b);
    tmp.set_lgt_parameter(3, px);tmp.set_lgt_parameter(4, py);tmp.set_lgt_parameter(5, pz);
    tmp.set_lgt_parameter(6, dx);tmp.set_lgt_parameter(7, dy);tmp.set_lgt_parameter(8, dz);
    tmp.set_lgt_parameter(9, angle1);tmp.set_lgt_parameter(10, angle2);
    Lgts.push_front(tmp);
    return true;
}
bool World::add_ambient_light(float r,float g,float b){
    //point_light r, g, b
    Lights tmp;
    tmp.set_lgt_name("ambient_light");
    tmp.set_lgt_parameter(0, r);tmp.set_lgt_parameter(1, g);tmp.set_lgt_parameter(2, b);
    Lgts.push_front(tmp);
    return true;
}
bool World::ray_tracing(){
    printf("ray tracing begin:_________________________\n");
    int width = W_Settings.film_resolution[0];
    int height = W_Settings.film_resolution[1];
    int percentage = -1;
    #pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP
    for(int i=0;i<width;i++)
        for(int j=0;j<height;j++)
        {
            if(percentage<(100*i*j)/(width*height)){
                percentage =(100*i*j)/(width*height);
                printf("processing: %d\n",percentage);
            }
            R_values R;
            if(check_intersect(R, Ray[index_ray(i, j)])==true){//R include objname, first intersect pos
                R.viewpoint[0] = W_Settings.camera[0];//set viewpoint to camera pos
                R.viewpoint[1] = W_Settings.camera[1];
                R.viewpoint[2] = W_Settings.camera[2];
                
                //input entity_it,viewpoint,hitpoint
                calcu_color(Scn[index_ray(i, j)].rgb, R,Ray[index_ray(i, j)],W_Settings.max_depth);
                Scn[index_ray(i, j)].rgb[0] *= 255;
                Scn[index_ray(i, j)].rgb[1] *= 255;
                Scn[index_ray(i, j)].rgb[2] *= 255;
            }else{
                Scn[index_ray(i, j)].rgb[0] = W_Settings.background[0]*255;
                Scn[index_ray(i, j)].rgb[1] = W_Settings.background[1]*255;
                Scn[index_ray(i, j)].rgb[2] = W_Settings.background[2]*255;
            }
        }

    printf("Raytracing finished\n");
    return true;
}
bool World::check_intersect(R_values &R, Rays ray){//ray collision check
    R.dis_in_t = -1;
    R_values T_r;
    bool collision_tag = false;
    float t=-1;
    for(std::list<Entities>::iterator it=Ets.begin();it!=Ets.end();++it){
        if(it->check_intersect(T_r, ray)){//return pos, distance
            if(t>T_r.dis_in_t||collision_tag ==false){//if nearer
                R=T_r;
                R.entity = (void*)&(*it);
                t = T_r.dis_in_t;
                //if(it->get_entity_name()=="triangle")
                //{
                  //  printf("collision with triangle\n");
                //}
            }
            collision_tag = true;
        }
    }
    return collision_tag;
}
bool World::calcu_color(float RGB[3],R_values R_v,Rays input_ray,int current_recursive_depth){
    if(pow(input_ray.Direction[0],2)+pow(input_ray.Direction[1],2)+pow(input_ray.Direction[2],2)-1>0.0001){
        float len = pow(input_ray.Direction[0],2)+pow(input_ray.Direction[1],2)+pow(input_ray.Direction[2],2);
        printf("error input ray len:%f\n",len);
    }
    
    //return: rgb,   input: entity pointer, viewpoint, hitpoint,
    //get input values from R_v;
    float hit_point[3],view_point[3];
    Entities *entity_it;
    entity_it = (Entities*)R_v.entity;
    hit_point[0] = R_v.pos[0];hit_point[1] = R_v.pos[1];hit_point[2] = R_v.pos[2];
    view_point[0] = R_v.viewpoint[0];view_point[1] = R_v.viewpoint[1];view_point[2] = R_v.viewpoint[2];
    float RGB_refl[3]={0,0,0},RGB_refr[3]={0,0,0},RGB_other[3]={0,0,0};
    R_values R_v1,R_v2,R_v3;
    R_v1 = R_v;R_v2 = R_v;R_v3 = R_v;
    
    calcu_refl_color(RGB_refl,  R_v1, input_ray, current_recursive_depth);
    calcu_refra_color(RGB_refr, R_v2, input_ray, current_recursive_depth);
    calcu_other_color(RGB_other,R_v3, input_ray);
    
    //return rgb
    //RGB_refl[0] = 0;RGB_refl[1]=0;RGB_refl[2]=0;
    RGB[0] = RGB_other[0]+RGB_refl[0]+RGB_refr[0];
    RGB[1] = RGB_other[1]+RGB_refl[1]+RGB_refr[1];
    RGB[2] = RGB_other[2]+RGB_refl[2]+RGB_refr[2];
    RGB[0] = RGB[0]<0?0:(RGB[0]>1?1:RGB[0]);
    RGB[1] = RGB[1]<0?0:(RGB[1]>1?1:RGB[1]);
    RGB[2] = RGB[2]<0?0:(RGB[2]>1?1:RGB[2]);
    return true;

}
bool World::calcu_refl_color(float RGB[3], R_values R_v,Rays input_ray,int current_recursive_depth){
    if(pow(input_ray.Direction[0],2)+pow(input_ray.Direction[1],2)+pow(input_ray.Direction[2],2)-1>0.0001){
        float len = pow(input_ray.Direction[0],2)+pow(input_ray.Direction[1],2)+pow(input_ray.Direction[2],2);
        printf("error input ray len:%f\n",len);
    }
    //return: rgb,   input: entity pointer, viewpoint, hitpoint,
    //get input values from R_v;
    float hit_point[3],view_point[3];
    Entities *entity_it;
    entity_it = (Entities*)R_v.entity;
    hit_point[0] = R_v.pos[0];hit_point[1] = R_v.pos[1];hit_point[2] = R_v.pos[2];
    view_point[0] = R_v.viewpoint[0];view_point[1] = R_v.viewpoint[1];view_point[2] = R_v.viewpoint[2];
    RGB[0] = 0;RGB[1] = 0;RGB[2] = 0;
    float ks_r,ks_g,ks_b;
    entity_it->get_property_material(6, ks_r);
    entity_it->get_property_material(7, ks_g);
    entity_it->get_property_material(8, ks_b);
    if(current_recursive_depth==0){//R_v will be reset in this time
        //end recursive until aim depth been reached or current ray can't hit any obj
        RGB[0] = W_Settings.background[0];RGB[1] = W_Settings.background[1];RGB[2] = W_Settings.background[2];
        //printf("rgb[0]:%f in calcu_refl color\n",RGB[0]);
        return true;
    }
    //recursively calculate reflection;
    R_values R_v_recursive;
    Rays ray_reflection;//shoot out from hit point
    ray_reflection = entity_it->get_reflection_ray(input_ray, hit_point);//input: viewpoint and hit point,output:ray
    float RGB_recursive[3]={0,0,0};
    float distance=0;
    
    if(check_intersect(R_v_recursive, ray_reflection)==true){//error in reflection
        //check wether the reflection ray could hit something or not
        R_v_recursive.viewpoint[0] = R_v.pos[0];//set new view point as current hit point
        R_v_recursive.viewpoint[1] = R_v.pos[1];
        R_v_recursive.viewpoint[2] = R_v.pos[2];
        //entitiy pointer and hitpoint are already saved in R_v_recursive
        distance = R_v_recursive.dis_in_t;//presave the distance, since it would be removed in calcu_color().
        calcu_color(RGB_recursive , R_v_recursive, ray_reflection, current_recursive_depth-1);
        //init rgb
        //RGB_recursive[0] = 255;
        //RGB_recursive[1] = 255;
        //RGB_recursive[2] = 255;
    }
    else{
        RGB_recursive[0] = W_Settings.background[0];
        RGB_recursive[1] = W_Settings.background[1];
        RGB_recursive[2] = W_Settings.background[2];
    }
    float delum_factor = 1;//deluminate factor for reflection
    RGB[0] = ks_r*delum_factor*RGB_recursive[0];
    RGB[1] = ks_g*delum_factor*RGB_recursive[1];
    RGB[2] = ks_b*delum_factor*RGB_recursive[2];
   
    RGB[0] = RGB[0]<0?0:(RGB[0]>1?1:RGB[0]);
    RGB[1] = RGB[1]<0?0:(RGB[1]>1?1:RGB[1]);
    RGB[2] = RGB[2]<0?0:(RGB[2]>1?1:RGB[2]);
    if(isnan(RGB[0])||isnan(RGB[1])||isnan(RGB[2])){
        printf("reflection: rgb: %f %f %f\n",RGB[0],RGB[1],RGB[2]);
    }
    
    return true;
}
bool World::calcu_refra_color(float RGB[3],R_values R_v,Rays input_ray,int current_recursive_depth){
    if(pow(input_ray.Direction[0],2)+pow(input_ray.Direction[1],2)+pow(input_ray.Direction[2],2)-1>0.0001){
        float len = pow(input_ray.Direction[0],2)+pow(input_ray.Direction[1],2)+pow(input_ray.Direction[2],2);
        printf("error input ray len:%f\n",len);
    }
    //return: rgb,   input: entity pointer, viewpoint, hitpoint,
    //get input values from R_v;
    float hit_point[3],view_point[3];
    Entities *entity_it;
    entity_it = (Entities*)R_v.entity;
    hit_point[0] = R_v.pos[0];hit_point[1] = R_v.pos[1];hit_point[2] = R_v.pos[2];
    view_point[0] = R_v.viewpoint[0];view_point[1] = R_v.viewpoint[1];view_point[2] = R_v.viewpoint[2];
    RGB[0] = 0;RGB[1] = 0;RGB[2] = 0;
    float kt_r,kt_g,kt_b;
    entity_it->get_property_material(10, kt_r);
    entity_it->get_property_material(11, kt_g);
    entity_it->get_property_material(12, kt_b);
    if(current_recursive_depth==0){//R_v will be reset in this time
        //end recursive until aim depth been reached or current ray can't hit any obj
        RGB[0] = W_Settings.background[0];RGB[1] = W_Settings.background[1];RGB[2] = W_Settings.background[2];
        return true;
    }
    //recursively calculate refraction;
    R_values R_v_recursive;
    Rays ray_refraction;//shoot out from hit point
    float ior_in,ior_out;
    //decide which one is input ior which one is output ior based on norm and input ray angle
    Rays norm = entity_it->get_normal_vec(hit_point);
    float temp = norm.Direction[0]*input_ray.Direction[0]+
                    norm.Direction[1]*input_ray.Direction[1]+
                    norm.Direction[2]*input_ray.Direction[2];
    float theta = acos(temp);
    if(theta>M_PI/2)//means from out to in
    {
        ior_in = 1;
        entity_it->get_property_material(13, ior_out);
    }else if(theta<M_PI/2)//means from in to out
    {
        entity_it->get_property_material(13, ior_in);
        ior_out = 1;
    }else{
        float len1 = pow(input_ray.Direction[0],2)+pow(input_ray.Direction[1],2)+pow(input_ray.Direction[2],2);
        float len2 = pow(norm.Direction[0],2)+pow(norm.Direction[1],2)+pow(norm.Direction[2],2);
        ior_in = 0;ior_out=0;
        printf("error in refraction should not shoot it\n");
        RGB[0] = 0;RGB[1] = 0;RGB[2] = 0;
        return true;
    } 
    ray_refraction = entity_it->get_refraction_ray(input_ray, hit_point, ior_in, ior_out);//input: viewpoint and hit point,output:ray
    float RGB_recursive[3]={0,0,0};
    float distance=0;
//
    //
    //
    //error in check intersection, wrong point output and direction should not be zero
    //
    //
    if(round(ray_refraction.range)!=-2&&check_intersect(R_v_recursive, ray_refraction)==true){//error in check_intersection
        //check wether the reflection ray could hit something or not
        R_v_recursive.viewpoint[0] = R_v.pos[0];//set new view point as current hit point
        R_v_recursive.viewpoint[1] = R_v.pos[1];
        R_v_recursive.viewpoint[2] = R_v.pos[2];
        //entitiy pointer and hitpoint are already saved in R_v_recursive
        distance = R_v_recursive.dis_in_t;//presave the distance, since it would be removed in calcu_color().
        calcu_color(RGB_recursive , R_v_recursive, ray_refraction, current_recursive_depth-1);
    }
    else{
        RGB_recursive[0] = W_Settings.background[0];
        RGB_recursive[1] = W_Settings.background[1];
        RGB_recursive[2] = W_Settings.background[2];
    }
    float delum_factor = 1;//deluminate factor for reflection
    RGB[0] = kt_r*delum_factor*RGB_recursive[0];
    RGB[1] = kt_g*delum_factor*RGB_recursive[1];
    RGB[2] = kt_b*delum_factor*RGB_recursive[2];

    RGB[0] = RGB[0]<0?0:(RGB[0]>1?1:RGB[0]);
    RGB[1] = RGB[1]<0?0:(RGB[1]>1?1:RGB[1]);
    RGB[2] = RGB[2]<0?0:(RGB[2]>1?1:RGB[2]);
    if(isnan(RGB[0])||isnan(RGB[1])||isnan(RGB[2])){
        printf("refraction: rgb: %f %f %f\n",RGB[0],RGB[1],RGB[2]);
    }
    //printf("rgb[0]: %fin calcu_refr color\n",RGB[0]);
    return true;
}
bool World::calcu_other_color(float RGB[3], R_values R_v,Rays input_ray){
    if(pow(input_ray.Direction[0],2)+pow(input_ray.Direction[1],2)+pow(input_ray.Direction[2],2)-1>0.0001){
        float len = pow(input_ray.Direction[0],2)+pow(input_ray.Direction[1],2)+pow(input_ray.Direction[2],2);
        printf("error input ray len\n");
    }
    //get input values from R_v;
    RGB[0] = 0;RGB[1] = 0;RGB[2] = 0;
    float hit_point[3],view_point[3];
    Entities *entity_it;
    entity_it = (Entities*)R_v.entity;
    hit_point[0] = R_v.pos[0];hit_point[1] = R_v.pos[1];hit_point[2] = R_v.pos[2];
    view_point[0] = R_v.viewpoint[0];view_point[1] = R_v.viewpoint[1];view_point[2] = R_v.viewpoint[2];
    //SET other variants
    //shading
    std::list<Lights>::iterator it;
    for(it=Lgts.begin();it!=Lgts.end();++it){
        std::string name = it->get_lgt_name();
        Rays ray = it->get_a_ray(hit_point);//ray to light
        R_values C_R;
        bool visible_tag=true;
        if(round(ray.range)!=-2){// if it's not ambient light,then check visibility
            visible_tag = false;
            if(check_intersect(C_R, ray)==true)//C_R include obj,first intersect pos, min distance
            {    //float dis = sqrt(pow(C_R.pos[0]-R_v.pos[0],2)+pow(C_R.pos[1]-R_v.pos[1],2)+pow(C_R.pos[2]-R_v.pos[2],2));
                if(C_R.dis_in_t<0.0001&&(void*)C_R.entity==(void*)entity_it)//collision of point itself
                    visible_tag=true;
            }
            else//no collision between light and obj
                visible_tag=true;
        }
        if(visible_tag){
            float tmp_RGB[3]={0,0,0};
            //printf("color viewpoint: %f %f %f\n",view_point[0],view_point[1],view_point[2]);
            it->calcu_color(tmp_RGB,entity_it,hit_point,view_point);
            
            
            RGB[0] += tmp_RGB[0]<0?0:tmp_RGB[0];
            RGB[1] += tmp_RGB[1]<0?0:tmp_RGB[1];
            RGB[2] += tmp_RGB[2]<0?0:tmp_RGB[2];
            if(isnan(RGB[0])||isnan(RGB[1])||isnan(RGB[2])){
                printf("other: rgb: %f %f %f\n",tmp_RGB[0],tmp_RGB[1],tmp_RGB[2]);
            }
        }
    }
    RGB[0] = RGB[0]<0?0:(RGB[0]>1?1:RGB[0]);
    RGB[1] = RGB[1]<0?0:(RGB[1]>1?1:RGB[1]);
    RGB[2] = RGB[2]<0?0:(RGB[2]>1?1:RGB[2]);
    
    return true;
}
void World::calcu_norm_vec(float &rx,float &ry,float &rz,float sx,float sy,float sz,float dx,float dy,float dz){
    rx = dx-sx;ry = dy-sy;rz = dz-sz;
    float l_d = sqrt(rx*rx+ry*ry+rz*rz);
    rx = rx/l_d;ry = ry/l_d;rz = rz/l_d;
}
bool World::generate_img(float sx,float sy){
    // notice: the width and height are flipped during computing so, here we flipped it back
    img =  new Image(W_Settings.film_resolution[1]/sy,W_Settings.film_resolution[0]/sx);
    //Scn = [width][height];
    for(int i=0;i<img->width;i++)
        for(int j = 0;j<img->height;j++)
        {
            Pixel p;
            if(Scn==nullptr)
                p.SetClamp(255, 255, 255);
            else{
                float r=0,g=0,b=0;
                for(int x=0;x<sx;x++)
                    for(int y=0;y<sy;y++)
                    {
                        r+=Scn[index_ray(j*sy+y, i*sx+x)].rgb[0];
                        g+=Scn[index_ray(j*sy+y, i*sx+x)].rgb[1];
                        b+=Scn[index_ray(j*sy+y, i*sx+x)].rgb[2];
                    }
                p.SetClamp(r/(sx*sy), g/(sx*sy),b/(sx*sy));
            }
            img->GetPixel(i, img->height-j-1) = p;
        }
    return true;
}
bool World::Enlarge_resolution(float sx,float sy){
    W_Settings.film_resolution[0]*=sx;W_Settings.film_resolution[1]*=sy;
    return true;
}
bool World::save_img(){
    char *out_name = new char(W_Settings.output_image.length()+1);
    strcpy(out_name,W_Settings.output_image.c_str());
    img->Write(out_name);
    return true;
}
std::string World::get_imagename(){
    return W_Settings.output_image;
}
World_Setting World::get_worldSetting(){
    return W_Settings;
}
bool World::Generate_ray(){
    bool cross_3(float *d,float *u,float *v);
    int width,height;
    width = W_Settings.film_resolution[0];
    height = W_Settings.film_resolution[1];
    Ray = new Rays [width * height];
    
    float z_up[3] = {W_Settings.camera[6],W_Settings.camera[7],W_Settings.camera[8]};//got axis z/up
    float x_drt[3] = {W_Settings.camera[3],W_Settings.camera[4],W_Settings.camera[5]};//got axis x/direction
    float y_left[3];//got axis y/left
    cross_3(y_left, z_up ,x_drt);
    float c_image[3];
    float dst = width/(2*tan(M_PI*W_Settings.camera[9]/180));
    
    for(int i=0;i<3;i++)
        c_image[i] = W_Settings.camera[i] + dst*x_drt[i];
    float Offset[3];
    for(int i=0;i<3;i++)
        Offset[i] = -(width/2.0)*z_up[i]-(height/2.0)*y_left[i];
    Scn = new Screen[width*height];
    for(int i=0;i<width;i++)
        for(int j=0;j<height;j++)
        {
            for(int k=0;k<3;k++){// calculate pixels position and ray
                Scn[index_ray(i, j)].position[k] = c_image[k]+Offset[k]+i*z_up[k]+j*y_left[k];
                Scn[index_ray(i, j)].rgb[k] = 0;
                // to make sure every thing afrond of the eye/camera could be seen
                // wee chose eye to be the start point of the ray instead points on screen
                // later, to meet the resolution requirement we will scale the image into that resolution
                Ray[index_ray(i, j)].Point[k] = W_Settings.camera[k];
                Ray[index_ray(i, j)].Direction[k] = Scn[index_ray(i, j)].position[k] - W_Settings.camera[k];
            }
            //normalize rays direction
            float l = sqrt(Ray[index_ray(i, j)].Direction[0]*Ray[index_ray(i, j)].Direction[0]+
                Ray[index_ray(i, j)].Direction[1]*Ray[index_ray(i, j)].Direction[1]+
                    Ray[index_ray(i, j)].Direction[2]*Ray[index_ray(i, j)].Direction[2]);
            for(int k=0;k<3;k++)
                Ray[index_ray(i, j)].Direction[k] = Ray[index_ray(i, j)].Direction[k]/l;
            //
        }
    return true;
}
bool cross_3(float *d,float *u,float *v){
    //cross(u, v)=(u2v3 - u3v2)e1 + (u3v1 - u1v3)e2 + (u1v2 - u2v1)e3;
    try {
        d[0] = u[2-1]*v[3-1] - u[3-1]*v[2-1];
        d[1] = u[3-1]*v[1-1] - u[1-1]*v[3-1];
        d[2] = u[1-1]*v[2-1] - u[2-1]*v[1-1];
    } catch (const std::exception& e) {
        printf("error in cross_3\n");
        return false;
    }

    return true;
}
int World::index_ray(int x_index,int y_index){//x: 0~width,y: 0~height
    return x_index*W_Settings.film_resolution[1]+y_index;
}
