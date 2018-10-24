//
//  Skybox.cpp
//  Raytracing
//
//  Created by 王登远 on 10/24/18.
//  Copyright © 2018 Dengyuan Wang. All rights reserved.
//

#include "Skybox.hpp"
bool Skybox::load_skybox(std::string filename)
{
    char *cstr = new char[filename.length() + 1];
    strcpy(cstr, filename.c_str());
    Image img(cstr);
    Image *Sky_up=nullptr,*Sky_down=nullptr,*Sky_front=nullptr,*Sky_back=nullptr,*Sky_left=nullptr,*Sky_right=nullptr;
    
    Sky_up = new Image(img.width/3,img.width/3);
    Sky_down = new Image(img.width/3,img.width/3);
    Sky_left = new Image(img.width/3,img.width/3);
    Sky_right = new Image(img.width/3,img.width/3);
    Sky_front = new Image(img.width/3,img.width/3);
    Sky_back = new Image(img.width/3,img.width/3);
    float upx=0,upy=0,downx=0,downy=0,leftx=0,lefty=0,rightx=0,righty=0,frontx=0,fronty=0,backx=0,backy=0;
    delete[] cstr;
    if(img.width<img.height)
    {
        for(int i=0;i<img.width;i++)
        {
            for(int j=0;j<img.height;j++)
            {
                
                int x = 3*i/img.width;//0,1,2
                int y = 4*j/img.height;//0,1,2,3
                if(x==1&&y==0)//right
                {
                    Pixel p = img.GetPixel(i, j);
                    if(rightx<img.width/3&&righty<img.width/3)
                        Sky_right->GetPixel(rightx, righty) = p;
                    if(++righty==Sky_right->height)
                    {
                        rightx++;righty=0;
                    }
                    
                }else if (x==0&&y==1)//front
                {
                    Pixel p = img.GetPixel(i, j);
                    if(frontx<img.width/3&&fronty<img.width/3)
                        Sky_front->GetPixel(frontx, fronty) = p;
                    if(++fronty==Sky_front->height)
                    {
                        frontx++;fronty=0;
                    }
                    
                }else if (x==1&&y==1)//down
                {
                    Pixel p = img.GetPixel(i, j);
                    if(downx<img.width/3&&downy<img.width/3)
                        Sky_down->GetPixel(downx, downy) = p;
                    if(++downy==Sky_down->height)
                    {
                        downx++;downy=0;
                    }
                    
                }else if (x==2&&y==1)//back
                {
                    Pixel p = img.GetPixel(i, j);
                    if(backx<img.width/3&&backy<img.width/3)
                        Sky_back->GetPixel(backx, backy) = p;
                    if(++backy==Sky_back->height)
                    {
                        backx++;backy=0;
                    }
                    
                }else if (x==1&&y==2)//left
                {
                    Pixel p = img.GetPixel(i, j);
                    if(leftx<img.width/3&&lefty<img.width/3)
                        Sky_left->GetPixel(leftx, lefty) = p;
                    if(++lefty==Sky_left->height)
                    {
                        leftx++;lefty=0;
                    }
                    
                }else if (x==1&&y==3)//Up
                {
                    Pixel p = img.GetPixel(i, j);
                    if(upx<img.width/3&&upy<img.width/3)
                        Sky_up->GetPixel(upx, upy) = p;
                    if(++upy==Sky_up->height)
                    {
                        upx++;upy=0;
                    }
                    
                }
            }
        }
        Sky_up = Sky_up->Scale(float(sky_size[0])/(img.width/3), float(sky_size[0])/(img.width/3));
        Sky_down = Sky_down->Scale(float(sky_size[0])/(img.width/3), float(sky_size[0])/(img.width/3));
        Sky_left = Sky_left->Scale(float(sky_size[0])/(img.width/3), float(sky_size[0])/(img.width/3));
        Sky_right = Sky_right->Scale(float(sky_size[0])/(img.width/3), float(sky_size[0])/(img.width/3));
        Sky_front = Sky_front->Scale(float(sky_size[0])/(img.width/3), float(sky_size[0])/(img.width/3));
        Sky_back = Sky_back->Scale(float(sky_size[0])/(img.width/3), float(sky_size[0])/(img.width/3));
        for(int i=0;i<sky_size[0];i++)
            for(int j=0;j<sky_size[1];j++)
            {
                Pixel p;
                float x,y,z;
                p = Sky_right->GetPixel(i, j);
                x = -sky_size[0]/2;y = sky_size[0]/2-j;z = sky_size[0]/2-i;
                add_point_light(p.r,p.g,p.b,x,y,z);
                
                p = Sky_left->GetPixel(i, j);
                x = sky_size[0]/2;y = sky_size[0]/2-j;z = i - sky_size[0]/2;
                add_point_light(p.r,p.g,p.b,x,y,z);
                
                p = Sky_front->GetPixel(i, j);
                x = sky_size[0]/2 - i;y = sky_size[0]/2-j;z = sky_size[0]/2;
                add_point_light(p.r,p.g,p.b,x,y,z);
                
                p = Sky_back->GetPixel(i, j);
                x = i-sky_size[0]/2;y = sky_size[0]/2-j;z = -sky_size[0]/2;
                add_point_light(p.r,p.g,p.b,x,y,z);
                
                p = Sky_up->GetPixel(i, j);
                x = sky_size[0]/2-j;y = sky_size[0]/2;z = sky_size[0]/2-i;
                add_point_light(p.r,p.g,p.b,x,y,z);
                p = Sky_down->GetPixel(i, j);
                x = j-sky_size[0]/2;y = -sky_size[0]/2;z = sky_size[0]/2-i;
                add_point_light(p.r,p.g,p.b,x,y,z);
                
            }
        delete Sky_up;
        delete Sky_down;
        delete Sky_left;
        delete Sky_right;
        delete Sky_front;
        delete Sky_back;
        return true;
    }
    return false;
}
bool Skybox::set_sky_size(float x,float y)
{
    x = x<y?x:y;
    sky_size[0] = x;
    sky_size[1] = x;
    return true;
}


