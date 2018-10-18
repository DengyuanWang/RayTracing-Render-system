//
//  World.cpp
//  Raytracing
//
//  Created by 王登远 on 10/9/18.
//  Copyright © 2018 Dengyuan Wang. All rights reserved.
//

#include "World.hpp"

//World class

World::World(string filename){
    
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
    string line;
    
    string fileName = W_Settings.setting_file;
    
    // open the file containing the scene description
    ifstream input(fileName);
    
    // check for errors in opening the file
    if(input.fail()){
        cout << "Can't open file '" << fileName << "'" << endl;
        return false;
    }
    
    // determine the file size (this is optional -- feel free to delete the 6 lines below)
    streampos begin,end;
    begin = input.tellg();
    input.seekg(0, ios::end);
    end = input.tellg();
    cout << "File '" << fileName << "' is: " << (end-begin) << " bytes long.\n\n";
    input.seekg(0, ios::beg);
    
    
    //Loop through reading each line
    string command;
    while(input >> command) { //Read first word in the line (i.e., the command type)
        
        if (command[0] == '#'){
            getline(input, line); //skip rest of line
            cout << "Skipping comment: " << command  << line <<  endl;
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
            string outFile;
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
            cout << "WARNING. Do not know command: " << command << endl;
        }
    }
    
    return true;
}
bool World::add_plane(float x,float y,float z,float dx,float dy,float dz){//one point and one n could decide a plane
    W_Objects_pt O_pt=nullptr,pt_head=Objs;
    if (pt_head==nullptr)//not W_Objects yet
    {
        Objs = new W_Objects;
        Objs->obj_name = W_Objects::plane;
        for(int i=0;i<14;i++)
            Objs->material[i] = W_Settings.material[i];
        Objs->obj_parameters = new float[6];
        Objs->obj_parameters[0] = x;
        Objs->obj_parameters[1] = y;
        Objs->obj_parameters[2] = z;
        Objs->obj_parameters[3] = dx;
        Objs->obj_parameters[4] = dy;
        Objs->obj_parameters[5] = dz;
        Objs->next = nullptr;
        return true;
    }else{//at least one object
        O_pt = pt_head->next;
        while(O_pt!=nullptr){
            pt_head = O_pt;
            O_pt = O_pt->next;
        }
        O_pt = new W_Objects;
        O_pt->obj_name = W_Objects::plane;
        for(int i=0;i<sizeof(W_Settings.material);i++)
            O_pt->material[i] = W_Settings.material[i];
        O_pt->obj_parameters = new float[6];
        O_pt->obj_parameters[0] = x;
        O_pt->obj_parameters[1] = y;
        O_pt->obj_parameters[2] = z;
        O_pt->obj_parameters[3] = dx;
        O_pt->obj_parameters[4] = dy;
        O_pt->obj_parameters[5] = dz;
        O_pt->next = nullptr;
        pt_head->next = O_pt;
        return true;
    }
}
bool World::add_sphere(float x,float y,float z,float r){
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
        
    }
}
bool World::add_directional_light(float r,float g,float b,float x,float y,float z){
    //point_light r, g, b, x, y, z
    //float point_light[6] = {255, 255, 255, 10, 10, 10};
    Lights_pt head = Lgts,l_pt;
    if(Lgts==nullptr){
        Lgts = new Lights;
        Lgts->lgt_name = Lights::directional_light;
        Lgts->lgt_parameters = new float[6];
        Lgts->lgt_parameters[0] = r;Lgts->lgt_parameters[1] = g;Lgts->lgt_parameters[2] = b;
        Lgts->lgt_parameters[3] = x;Lgts->lgt_parameters[4] = y;Lgts->lgt_parameters[5] = z;
        Lgts->next = nullptr;
        return true;
    }else{
        while(head->next!=nullptr){
            head = head->next;
        }
        l_pt = new Lights;
        l_pt->lgt_name = Lights::directional_light;
        l_pt->lgt_parameters = new float[6];
        l_pt->lgt_parameters[0] = r;l_pt->lgt_parameters[1] = g;l_pt->lgt_parameters[2] = b;
        l_pt->lgt_parameters[3] = x;l_pt->lgt_parameters[4] = y;l_pt->lgt_parameters[5] = z;
        l_pt->next = nullptr;
        head->next = l_pt;
        return true;
    }
    return false;
}
bool World::add_point_light(float r,float g,float b,float x,float y,float z){
    //point_light r, g, b, x, y, z
    //float point_light[6] = {255, 255, 255, 10, 10, 10};
    Lights_pt head = Lgts,l_pt;
    if(Lgts==nullptr){
        Lgts = new Lights;
        Lgts->lgt_name = Lights::point_light;
        Lgts->lgt_parameters = new float[6];
        Lgts->lgt_parameters[0] = r;Lgts->lgt_parameters[1] = g;Lgts->lgt_parameters[2] = b;
        Lgts->lgt_parameters[3] = x;Lgts->lgt_parameters[4] = y;Lgts->lgt_parameters[5] = z;
        Lgts->next = nullptr;
        return true;
    }else{
        while(head->next!=nullptr){
            head = head->next;
        }
        l_pt = new Lights;
        l_pt->lgt_name = Lights::point_light;
        l_pt->lgt_parameters = new float[6];
        l_pt->lgt_parameters[0] = r;l_pt->lgt_parameters[1] = g;l_pt->lgt_parameters[2] = b;
        l_pt->lgt_parameters[3] = x;l_pt->lgt_parameters[4] = y;l_pt->lgt_parameters[5] = z;
        l_pt->next = nullptr;
        head->next = l_pt;
        return true;
    }
    return false;
}
bool World::add_spot_light(float r,float g,float b,float px,float py,float pz,float dx,float dy,float dz,float angle1,float angle2){
    //point_light r, g, b, x, y, z
    //float point_light[6] = {255, 255, 255, 10, 10, 10};
    Lights_pt head = Lgts,l_pt;
    if(Lgts==nullptr){
        Lgts = new Lights;
        Lgts->lgt_name = Lights::spot_light;
        Lgts->lgt_parameters = new float[11];
        Lgts->lgt_parameters[0] = r;Lgts->lgt_parameters[1] = g;Lgts->lgt_parameters[2] = b;
        Lgts->lgt_parameters[3] = px;Lgts->lgt_parameters[4] = py;Lgts->lgt_parameters[5] = pz;
        Lgts->lgt_parameters[6] = dx;Lgts->lgt_parameters[7] = dy;Lgts->lgt_parameters[8] = dz;
        Lgts->lgt_parameters[9] = angle1;Lgts->lgt_parameters[10] = angle2;
        Lgts->next = nullptr;
        return true;
    }else{
        while(head->next!=nullptr){
            head = head->next;
        }
        l_pt = new Lights;
        l_pt->lgt_name = Lights::spot_light;
        l_pt->lgt_parameters = new float[11];
        l_pt->lgt_parameters[0] = r;l_pt->lgt_parameters[1] = g;l_pt->lgt_parameters[2] = b;
        l_pt->lgt_parameters[3] = px;l_pt->lgt_parameters[4] = py;l_pt->lgt_parameters[5] = pz;
        l_pt->lgt_parameters[6] = dx;l_pt->lgt_parameters[7] = dy;l_pt->lgt_parameters[8] = dz;
        l_pt->lgt_parameters[9] = angle1;l_pt->lgt_parameters[10] = angle2;
        l_pt->next = nullptr;
        head->next = l_pt;
        return true;
    }
    return false;
}
bool World::add_ambient_light(float r,float g,float b){
    //point_light r, g, b
    Lights_pt head = Lgts,l_pt;
    if(Lgts==nullptr){
        Lgts = new Lights;
        Lgts->lgt_name = Lights::ambient_light;
        Lgts->lgt_parameters = new float[3];
        Lgts->lgt_parameters[0] = r;Lgts->lgt_parameters[1] = g;Lgts->lgt_parameters[2] = b;
        Lgts->next = nullptr;
        return true;
    }else{
        while(head->next!=nullptr){
            //printf("head->lgt_name: %d",head->lgt_name);
            head = head->next;
        }
        l_pt = new Lights;
        l_pt->lgt_name = Lights::ambient_light;
        l_pt->lgt_parameters = new float[3];
        l_pt->lgt_parameters[0] = r;l_pt->lgt_parameters[1] = g;l_pt->lgt_parameters[2] = b;
        l_pt->next = nullptr;
        head->next = l_pt;
        return true;
    }
    return false;
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
            
            /* typedef struct R_values{
                    W_Objects_pt obj;
                    Lights lgt;
                    float pos[3];
                    float rgb[3];// range 0~1
                }R_values;*/
            R_values R;
            if(check_intersect(R, Ray[index_ray(i, j)])==true){
                calcu_color(R,W_Settings.max_depth);
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
    W_Objects_pt header=Objs;
    float t = 0;
    bool collision_tag = false;
    while(header!=nullptr){
        if(header->obj_name==W_Objects::sphere){
            float x,x0,y,y0,z,z0,dx,dy,dz,r,t1,t2;
            x = ray.Point[0];y = ray.Point[1];z = ray.Point[2];
            dx = ray.Direction[0];dy = ray.Direction[1];dz = ray.Direction[2];
            x0 = header->obj_parameters[0];
            y0 = header->obj_parameters[1];
            z0 = header->obj_parameters[2];
            r  = header->obj_parameters[3];
            float a,b,c;
            a = dx*dx+dy*dy+dz*dz;
            b = 2.0*(dx*(x-x0)+dy*(y-y0)+dz*(z-z0));
            c = (x-x0)*(x-x0) + (y-y0)*(y-y0) + (z-z0)*(z-z0) - r*r;
            if(b*b-4.0*a*c>=0)//intersect
            {
                t1 = (-b + sqrt(b*b-4.0*a*c))/(2.0*a);
                t2 = (-b - sqrt(b*b-4.0*a*c))/(2.0*a);
                if(fmax(t1,t2)<=0){// no leagle one
                    
                }
                else if(ray.range>=0&&fmin(t1,t2)>ray.range)//check for segment. for ray the range will be -1, mean disable;
                {
                    
                }
                else if(t1<=0){//only t2
                    if(collision_tag==false||t>=t2){
                        t = t2;
                        R.obj = header;
                        R.pos[0] = x+t*dx;R.pos[1] = y+t*dy;R.pos[2] = z+t*dz;
                        collision_tag = true;
                    }
                }
                else if(t2<=0){//only t1
                    if(collision_tag==false||t>=t1){
                        t = t1;
                        R.obj = header;
                        R.pos[0] = x+t*dx;R.pos[1] = y+t*dy;R.pos[2] = z+t*dz;
                        collision_tag = true;
                    }
                }
                else{//intersect two points
                    if(collision_tag==false){//no former collision
                        if(t1<t2){//save nearest one
                            t = t1;
                            R.obj = header;
                            R.pos[0] = x+t*dx;R.pos[1] = y+t*dy;R.pos[2] = z+t*dz;
                        }else{
                            t = t2;
                            R.obj = header;
                            R.pos[0] = x+t*dx;R.pos[1] = y+t*dy;R.pos[2] = z+t*dz;
                        }
                        collision_tag = true;
                    }else{//have former collision
                        if(t1<=t2&&t1<=t){//save the smallest one between these two and former nearest one
                            t = t1;
                            R.obj = header;
                            R.pos[0] = x+t*dx;R.pos[1] = y+t*dy;R.pos[2] = z+t*dz;
                        }else if(t2<=t&&t2<=t1){
                            t = t2;
                            R.obj = header;
                            R.pos[0] = x+t*dx;R.pos[1] = y+t*dy;R.pos[2] = z+t*dz;
                        }
                    }

                }
                
            }
        }else{
            printf("error obj name:%d in check_intersect\n",header->obj_name);
            return false;
        }
        header = header->next;
    }
    return collision_tag;
}
bool World::calcu_color(R_values &R_v,int current_recursive_depth){
    Pixel p;
    float R = 0,G = 0,B = 0;
    if(current_recursive_depth==0){
        R_v.rgb[0] = R>1?1:R;R_v.rgb[1] = G>1?1:G;R_v.rgb[2] = B>1?1:B;
        return true;
    }
    //shading
    Lights_pt L_header = Lgts;
    while(L_header!=nullptr){
        if(L_header->lgt_name==Lights::point_light){
            calcu_rgb_under_lights(R,G,B,R_v,L_header,true);
        }
        else if(L_header->lgt_name==Lights::directional_light){// add the directional_light to the object
            calcu_rgb_under_lights(R,G,B,R_v,L_header,false);
        }
        else if(L_header->lgt_name==Lights::ambient_light){// add the ambient_light to the object
            float ar = R_v.obj->material[0];float ag = R_v.obj->material[1];float ab = R_v.obj->material[2];
            R+=ar*L_header->lgt_parameters[0];
            G+=ag*L_header->lgt_parameters[1];
            B+=ab*L_header->lgt_parameters[2];
        }
        L_header = L_header->next;
    }
    //reflection:
    R_values R_v2;
    R_v2 = R_v;
    Rays ray_reflect;
    ray_reflect.Point[0] = R_v.pos[0];ray_reflect.Point[1] = R_v.pos[1];ray_reflect.Point[2] = R_v.pos[2];
    //check_intersect(R, Ray[index_ray(i, j)]
    calcu_color(R_v2,current_recursive_depth-1);//value saved in R_v
    //Refrection:
    
    R_v.rgb[0] += R>1?1:R;R_v.rgb[1] += G>1?1:G;R_v.rgb[2] += B>1?1:B;
    return true;
}
void World::calcu_rgb_under_lights(float &R,float &G,float &B,R_values R_v,Lights_pt L_header,bool delum_tag){
    float dr = R_v.obj->material[3];float dg = R_v.obj->material[4];float db = R_v.obj->material[5];
    float sr = R_v.obj->material[6];float sg = R_v.obj->material[7];float sb = R_v.obj->material[8];
    float x = L_header->lgt_parameters[3],y = L_header->lgt_parameters[4],z = L_header->lgt_parameters[5];
    R_values C_R;
    Rays ray;
    bool visible_tag = false;
    // normal direction
    // this ray is from light to point
    if(delum_tag)
        calcu_norm_vec(ray.Direction[0],ray.Direction[1],ray.Direction[2],
                   R_v.pos[0],R_v.pos[1],R_v.pos[2],x,y,z);
    else
        calcu_norm_vec(ray.Direction[0],ray.Direction[1],ray.Direction[2],
                       2*x,2*y,2*z,x,y,z);
    ray.Point[0] = R_v.pos[0]+0.01*ray.Direction[0];
    ray.Point[1] = R_v.pos[1]+0.01*ray.Direction[1];
    ray.Point[2] = R_v.pos[2]+0.01*ray.Direction[2];
    if(delum_tag)
        ray.range = (x-R_v.pos[0])/ray.Direction[0];
    else
        ray.range = -1;
    if(check_intersect(C_R, ray)==true){
        float dis = sqrt(pow(C_R.pos[0]-R_v.pos[0],2)+pow(C_R.pos[1]-R_v.pos[1],2)+pow(C_R.pos[2]-R_v.pos[2],2));
        if(dis<0.0001){//collision of itself
            //G = 255;
            visible_tag=true;
        }
    }
    else//no collision between light and obj
        visible_tag=true;
    if(visible_tag)// the object can see at least one light
    {
        float delum_factor;
        if(delum_tag)
            delum_factor = 1/(pow(R_v.pos[0]-x,2)+pow(R_v.pos[1]-y,2)+pow(R_v.pos[2]-z,2));
        else
            delum_factor = 1;
        
        //Lambertian shading
        
        float n[3];
        calcu_norm_vec(n[0],n[1],n[2],
                       R_v.obj->obj_parameters[0],R_v.obj->obj_parameters[1],R_v.obj->obj_parameters[2]
                       ,R_v.pos[0],R_v.pos[1],R_v.pos[2]);
        float delta_r,delta_g,delta_b,theta;
        theta = 180*acos(n[0]*ray.Direction[0]+n[1]*ray.Direction[1]+n[2]*ray.Direction[2])/M_PI;
        //printf("theta = %f\n",theta);
        //printf("dr=%f dg=%f db=%f; lr=%f lg=%f lb=%f\n",dr,dg,db,L_header->lgt_parameters[0]
        // ,L_header->lgt_parameters[1],L_header->lgt_parameters[2]);
        delta_r = dr*delum_factor*L_header->lgt_parameters[0]*fmax(0,cos(M_PI*theta/180));
        delta_g = dg*delum_factor*L_header->lgt_parameters[1]*fmax(0,cos(M_PI*theta/180));
        delta_b = db*delum_factor*L_header->lgt_parameters[2]*fmax(0,cos(M_PI*theta/180));
        R += delta_r;
        G += delta_g;
        B += delta_b;
        //phong shading specularity
        //Is_r = dr*I_r*dot(v,R)^n_
        //R = I - 2*N*dot(I*N)
        float v[3];
        float R_[3];
        float n_ = R_v.obj->material[9];
        calcu_norm_vec(v[0],v[1],v[2],
                       R_v.obj->obj_parameters[0],R_v.obj->obj_parameters[1],R_v.obj->obj_parameters[2]
                       ,W_Settings.camera[0],W_Settings.camera[1],W_Settings.camera[2]);
        float dot = (ray.Direction[0]*n[0]+ray.Direction[1]*n[1]+ray.Direction[2]*n[2]);
        R_[0] = ray.Direction[0]-2*dot*n[0];
        R_[1] = ray.Direction[1]-2*dot*n[1];
        R_[2] = ray.Direction[2]-2*dot*n[2];
        dot = pow(v[0]*R_[0]+v[1]*R_[1]+v[2]*R_[2],n_);
        
        R+=sr*delum_factor*L_header->lgt_parameters[0]*dot;
        G+=sg*delum_factor*L_header->lgt_parameters[1]*dot;
        B+=sb*delum_factor*L_header->lgt_parameters[2]*dot;
    }
    
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
string World::get_imagename(){
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
World::~World(){
    delete img;
    img=nullptr;
//    World_Setting W_Settings;
    delete Ray;
    Ray=nullptr;
    delete Scn;
    Scn=nullptr;
    if(Objs!=nullptr){
        W_Objects_pt tail,head=Objs;
        do {
            tail = head->next;
            head->next = nullptr;
            delete head->obj_parameters;
            delete head;
            head = tail;
        }while(tail!=nullptr);
        Objs = nullptr;
    }
    if(Lgts!=nullptr){
        Lights_pt tail,head=Lgts;
        do {
            tail = head->next;
            head->next = nullptr;
            delete head->lgt_parameters;
            delete head;
            head = tail;
        }while(tail!=nullptr);
        Lgts = nullptr;
    }
}
