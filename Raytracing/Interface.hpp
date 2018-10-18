//
//  Interface.hpp
//  Raytracing
//
//  Created by 王登远 on 10/9/18.
//  Copyright © 2018 Dengyuan Wang. All rights reserved.
//

#ifndef Interface_hpp
#define Interface_hpp

#ifndef Default_setting_h
#define Default_setting_h
#include "Default_setting.h"

#endif
#include "SDL2/SDL.h"
#include <stdio.h>
#include <iostream>

class Interface {
    SDL_Surface *imageSurface = nullptr;
    SDL_Surface *windowSurface = nullptr;
    SDL_Window *window=nullptr;
    SDL_Event window_event;
    int WIDTH =800,HEIGHT = 600;
public:
    Interface();
    bool create_window(int width,int height);
    
    bool load_img(std::string filename);
    
    bool update_window();
    
    bool close_window();
    
    ~Interface();
};



#endif /* Interface_hpp */
