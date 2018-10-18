//
//  Interface.cpp
//  Raytracing
//
//  Created by 王登远 on 10/9/18.
//  Copyright © 2018 Dengyuan Wang. All rights reserved.
//

#include "Interface.hpp"
// Interface class

Interface::Interface(){

}

bool Interface::create_window(int width,int height){
    WIDTH = width;
    HEIGHT = height;
    //Create a window
    if(SDL_Init(SDL_INIT_EVERYTHING)<0){//init SDL
        std::cout<<"SDL could not initialization! SDL ERROR: "<<SDL_GetError()<<std::endl;
        exit(EXIT_FAILURE);
    }
    window = SDL_CreateWindow("Hello World!", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, WIDTH, HEIGHT, SDL_WINDOW_OPENGL);
    windowSurface = SDL_GetWindowSurface(window);//get window surface from window
    if( window == nullptr){
        std::cout<<"SDL could not create window: "<<SDL_GetError()<<std::endl;
        return EXIT_FAILURE;
    }
    
    
    return false;
}
bool Interface::load_img(std::string filename){
    //load Image
    imageSurface = SDL_LoadBMP(filename.c_str());
    if(imageSurface == nullptr){
        std::cout<<"SDL could not load image: "<<SDL_GetError()<<std::endl;
        return EXIT_FAILURE;
    }
    return true;
}
bool Interface::update_window(){//true for update_end
    if(SDL_PollEvent(&window_event)){
        if(SDL_QUIT==window_event.type){
            return true;
        }
    }
    //transfer image into windows surface
    SDL_BlitSurface(imageSurface, nullptr, windowSurface, nullptr);
    //update display
    SDL_UpdateWindowSurface(window);
    return false;
}
bool Interface::close_window(){
    // clear all content
    SDL_FreeSurface(imageSurface);
    SDL_FreeSurface(windowSurface);
    imageSurface = nullptr;
    windowSurface = nullptr;
    SDL_DestroyWindow(window);
    window = nullptr;
    SDL_Quit();
    return false;
}
Interface::~Interface(){
    
}
