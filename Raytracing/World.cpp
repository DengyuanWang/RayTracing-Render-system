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
    if(!build_world())
        printf("build world error\n");
}
bool World::build_world(){
    //add two wall
    
    //add ground
    
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
            printf("max_deptj of: %f\n", max_depth);
        }
        else {
            getline(input, line); //skip rest of line
            std::cout << "WARNING. Do not know command: " << command << std::endl;
        }
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
    /*
    W_Objects_pt O_pt=nullptr,pt_head=Objs;
    if (pt_head==nullptr)//not W_Objects yet
    {
        Objs = new W_Objects;
        Objs->obj_name = W_Objects::sphere;
        for(int i=0;i<14;i++)
            Objs->material[i] = W_Settings.material[i];
        Objs->obj_parameters = new float[4];
        Objs->obj_parameters[0] = x;
        Objs->obj_parameters[1] = y;
        Objs->obj_parameters[2] = z;
        Objs->obj_parameters[3] = r;
        Objs->next = nullptr;
        return true;
    }else{//at least one object
        O_pt = pt_head->next;
        while(O_pt!=nullptr){
            pt_head = O_pt;
            O_pt = O_pt->next;
        }
        O_pt = nullptr;
        O_pt = new W_Objects;
        O_pt->obj_name = W_Objects::sphere;
        for(int i=0;i<sizeof(W_Settings.material);i++)
            O_pt->material[i] = W_Settings.material[i];
        O_pt->obj_parameters = nullptr;
        O_pt->obj_parameters = new float[4];
        O_pt->obj_parameters[0] = x;
        O_pt->obj_parameters[1] = y;
        O_pt->obj_parameters[2] = z;
        O_pt->obj_parameters[3] = r;
        O_pt->next = nullptr;
        pt_head->next = O_pt;
        return true;
    }*/
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
    for(int i=0;i<width;i++)
        for(int j=0;j<height;j++)
        {
            if(percentage<(100*i*j)/(width*height)){
                percentage =(100*i*j)/(width*height);
                printf("processing: %d\n",percentage);
            }
            R_values R;
            if(check_intersect(R, Ray[index_ray(i, j)])==true){//R include objname, first intersect pos
                calcu_color(R,Ray[index_ray(i, j)],W_Settings.max_depth);
                Scn[index_ray(i, j)].rgb[0] = R.rgb[0]*255;
                Scn[index_ray(i, j)].rgb[1] = R.rgb[1]*255;
                Scn[index_ray(i, j)].rgb[2] = R.rgb[2]*255;
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

    R_values T_r;
    bool collision_tag = false;
    float t=-1;
    for(std::list<Entities>::iterator it=Ets.begin();it!=Ets.end();++it){
        collision_tag = collision_tag||it->check_intersect(T_r, ray);
        if(t>T_r.dis_in_t){
            R=T_r;
            R.entity = &(*it);
            t = T_r.dis_in_t;
        }
    }
    return collision_tag;
}
bool World::calcu_color(R_values &R_v,Rays ray,int current_recursive_depth){
    Pixel p;
    R_values In_v=R_v;
    Entities *R_it = (Entities*)R_v.entity;
    
    float R = 0,G = 0,B = 0;
    if(current_recursive_depth==0||check_intersect(In_v, ray)==false){//R_v will be reset in this time;
        //end recursive until aim depth or current ray can't hit any obj
        R_v.rgb[0] = R>1?1:R;R_v.rgb[1] = G>1?1:G;R_v.rgb[2] = B>1?1:B;
        return true;
    }else{// recursively calculate reflection;
        R_values R_v_recursive;
        
        Rays ray_reflection = R_it->get_reflection_ray(R_v);//input: viewpoint and hit point,output:ray
        calcu_color(R_v_recursive, ray_reflection, current_recursive_depth-1);
        //init rgb
        float kr_r,kr_g,kr_b;
        Entities *IN_it = (Entities*)In_v.entity;
        IN_it->get_property_material(6, kr_r);
        IN_it->get_property_material(7, kr_g);
        IN_it->get_property_material(8, kr_b);
        In_v.rgb[0] = kr_r*R_v_recursive.rgb[0];
        In_v.rgb[1] = kr_g*R_v_recursive.rgb[1];
        In_v.rgb[2] = kr_b*R_v_recursive.rgb[2];
    }
    //shading
    std::list<Lights>::iterator it;
    for(it=Lgts.begin();it!=Lgts.end();++it){
        Rays ray = it->get_a_ray(In_v);//ray to light
        R_values C_R;
        bool visible_tag=false;
        if(check_intersect(C_R, ray)==true){//C_R include obj,first intersect pos, min distance
            //float dis = sqrt(pow(C_R.pos[0]-R_v.pos[0],2)+pow(C_R.pos[1]-R_v.pos[1],2)+pow(C_R.pos[2]-R_v.pos[2],2));
            if(C_R.dis_in_t<0.0001&&C_R.entity==In_v.entity){//collision of point itself
                visible_tag=true;
            }
        }
        else//no collision between light and obj
            visible_tag=true;
        if(visible_tag){
            In_v.viewpoint[0] = W_Settings.camera[0];In_v.viewpoint[1] = W_Settings.camera[1];In_v.viewpoint[2] = W_Settings.camera[2];
            it->calcu_color(R, G, B, In_v);
        }
    }
    R += In_v.rgb[0];
    G += In_v.rgb[1];
    B += In_v.rgb[2];
    R_v.rgb[0] = R>1?1:R;R_v.rgb[1] = G>1?1:G;R_v.rgb[2] = B>1?1:B;
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
