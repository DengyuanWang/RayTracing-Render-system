#include "image.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <algorithm>
#include <iostream>
/**
 * Image
 **/
Image::Image (int width_, int height_){

    assert(width_ > 0);
    assert(height_ > 0);

    width           = width_;
    height          = height_;
    num_pixels      = width * height;
    sampling_method = IMAGE_SAMPLING_POINT;
    
    data.raw = new uint8_t[num_pixels*4];
	int b = 0; //which byte to write to
	for (int j = 0; j < height; j++){
		for (int i = 0; i < width; i++){
			data.raw[b++] = 0;
			data.raw[b++] = 0;
			data.raw[b++] = 0;
			data.raw[b++] = 0;
		}
	}

    assert(data.raw != NULL);
}

Image::Image (const Image& src){
	
	width           = src.width;
    height          = src.height;
    num_pixels      = width * height;
    sampling_method = IMAGE_SAMPLING_POINT;
    
    data.raw = new uint8_t[num_pixels*4];
    
    //memcpy(data.raw, src.data.raw, num_pixels);
    *data.raw = *src.data.raw;
}

Image::Image (char* fname){
	int numComponents; //(e.g., Y, YA, RGB, or RGBA)
	data.raw = stbi_load(fname, &width, &height, &numComponents, 4);
	
	if (data.raw == NULL){
		printf("Error loading image: %s", fname);
		exit(-1);
	}
	

	num_pixels = width * height;
	sampling_method = IMAGE_SAMPLING_POINT;
	
}

Image::Image (char* fname,char* type){
    int encode_channel(int channel);
    FILE *in_file = fopen(fname, "r");
    if (! in_file ) // equivalent to saying if ( in_file == NULL )
    {
        printf("oops, file can't be read\n");
        exit(-1);
    }
    fscanf(in_file, "%d %d ", & width,&height);
    Image* img2 = new Image(width,height);
    int r,g,b;
    char pixel;
    Pixel p;
    // attempt to read the next line and store
    // the value in the "number" variable
    for (int i=0;i<width;i++)
        for(int j=0;j<height;j++)
        {
            fscanf(in_file, "%c", & pixel );
            r = pixel/25;
            g = (pixel-r*25)/5;
            b = pixel-r*25-g*5;
            p.SetClamp(encode_channel(r),
                       encode_channel(g),
                       encode_channel(b));
            img2->GetPixel(i, j) = p;
        }
    sampling_method = IMAGE_SAMPLING_POINT;
    img2->sampling_method = 1;
    Image *img = img2->Scale(2, 2);
    width*=2;
    height*=2;
    num_pixels = width * height;
    data.raw = new uint8_t[num_pixels*4];
    
    for (int i=0;i<width;i++)
        for(int j=0;j<height;j++)
        {
            GetPixel(i, j) = img->GetPixel(i, j);
        }

    fclose(in_file);
}



Image::~Image (){
    delete data.raw;
    data.raw = NULL;
}

void Image::Write(char* fname){
	
	int lastc = strlen(fname);

	switch (fname[lastc-1]){
	   case 'g': //jpeg (or jpg) or png
	     if (fname[lastc-2] == 'p' || fname[lastc-2] == 'e') //jpeg or jpg
	        stbi_write_jpg(fname, width, height, 4, data.raw, 95);  //95% jpeg quality
	     else //png
	        stbi_write_png(fname, width, height, 4, data.raw, width*4);
	     break;
	   case 'a': //tga (targa)
	     stbi_write_tga(fname, width, height, 4, data.raw);
	     break;
	   case 'p': //bmp
	   default:
	     stbi_write_bmp(fname, width, height, 4, data.raw);
	}
}
void Image::Lossy_compress_Write(char* fname)
{
    FILE *out_file = fopen(fname, "w");
    char **pixel_format;
    pixel_format = new char*[width];
    int decode_channel(int rgb);
    
    fprintf(out_file, "%d %d ", width/2,height/2);

    Image *img = Scale(0.5, 0.5);
    img->OrderedDither(2);
    for(int i=0;i<width/2;i++)
    {
        pixel_format[i] = new char[height];
        for(int j=0;j<height/2;j++)
        {
            Pixel p = img->GetPixel(i, j);
            
            pixel_format[i][j] = 25*decode_channel(p.r)
                                   +5*decode_channel(p.g)
                                    +1*decode_channel(p.b);
            fprintf(out_file, "%c", pixel_format[i][j]);
        }
    }
    
    fclose(out_file);
}
int encode_channel(int channel)
{
    int nbits = 2;
    int *color_list = new int[pow(2,nbits)];
    for(int i=0;i<sizeof(color_list);i++)
        color_list[i] = i*pow(2.0,8-nbits)+pow(2.0,8-nbits-1);
    channel =color_list[channel];
    delete [] color_list;
    return channel;
}
int decode_channel(int rgb)
{
    int nbits = 2;
    int *color_list = new int[pow(2,nbits)];
    for(int i=0;i<sizeof(color_list);i++)
        color_list[i] = i*pow(2.0,8-nbits)+pow(2.0,8-nbits-1);
    int i = 0;
    while(rgb!=color_list[i]){
        i++;
    }
    delete [] color_list;
    return i;
}
void Image::AddNoise (double factor)
{
	/* WORK HERE */
    int x,y;
    for (x = 0 ; x < Width() ; x++)
    {
        for (y = 0 ; y < Height() ; y++)
        {
            Pixel p = GetPixel(x, y);
            Pixel rand_p = PixelRandom();
            GetPixel(x,y) = PixelLerp (p, rand_p, factor);
        }
    }
}

void Image::Brighten (double factor)
{
	int x,y;
	for (x = 0 ; x < Width() ; x++)
	{
		for (y = 0 ; y < Height() ; y++)
		{
			Pixel p = GetPixel(x, y);
			Pixel scaled_p = p*factor;
			GetPixel(x,y) = scaled_p;
		}
	}
}


void Image::ChangeContrast (double factor)
{
    int x,y;
    Component max_l=0,min_l=255;
    float sum_l=0;
    double sum_c=0;
    for (x = 0 ; x < width ; x++)
    {
        for (y = 0 ; y < height ; y++)
        {
            Pixel p = GetPixel(x, y);
            if (max_l<p.Luminance())
                max_l =p.Luminance();
            if (min_l>p.Luminance())
                min_l=p.Luminance();
            sum_l+=float(p.Luminance());
        }
    }
    for (x = 0 ; x < width ; x++)
    {
        for (y = 0 ; y < height ; y++)
        {
            Pixel p = GetPixel(x, y);
            
            sum_c += pow(float(p.Luminance())-sum_l,2);
        }
    }
    float contrast =  (float(max_l)-float(min_l))/float(max_l)+float(min_l);
    //contrast = sqrt(sum_c/((width+1)*(height+1)));
    float factor_c = factor*(259.0 * (contrast + 255.0)) / (255.0 * ((259.0 - contrast)));
    for (x = 0 ; x < Width() ; x++)
    {
        for (y = 0 ; y < Height() ; y++)
        {
            Pixel p = GetPixel(x, y);
            p.SetClamp(128+(p.r-128)*factor_c, 128+(p.g-128)*factor_c, 128+(p.b-128)*factor_c);
            GetPixel(x, y) = p;
        }
    }

}


void Image::ChangeSaturation(double factor)
{
	/* WORK HERE */
    float h = 0, r;          // hue shift (in degrees)
    float s = factor, g;        // saturation multiplier (scalar)
    float v = 1, b;      // value multiplier (scalar)
    float vsu;
    float vsw;
    for (int x = 0 ; x < width ; x++)
    {
        for (int y = 0 ; y < height ; y++)
        {
            //Using YIQ color space as a intermedia
            /*
             T_YIQ = [.299 .587 .114
                        .596  -.274 -.321
                        .211  -.523 .311  ]
             T_RGB = [1 .956 .621
                      1 −.272−.647
                      1 -1.107 1.705]
             T_hue = [1 0 0
                     0 U -W
                     0 W U  ] where U = cos(theta),W = sin(theta),theta = H*pi/180
             Then :
                           [ V  0    0
             T_HSV = T_RGB*  0 VSU -VSW  * T_YIQ
                             0 VSW VSU  ]
             RGB2 = T_HSV*RGB1;
             */
            Pixel p = GetPixel(x, y);
            vsu = v*s*cos(h*M_PI/180);
            vsw = v*s*sin(h*M_PI/180);
            r = (.299*v + .701*vsu + .168*vsw)*p.r
               +(.587*v - .587*vsu + .330*vsw)*p.g
               +(.114*v - .114*vsu - .497*vsw)*p.b;
            g = (.299*v - .299*vsu - .328*vsw)*p.r
               +(.587*v + .413*vsu + .035*vsw)*p.g
               +(.114*v - .114*vsu + .292*vsw)*p.b;
            b = (.299*v - .300*vsu + 1.25*vsw)*p.r
               +(.587*v - .588*vsu - 1.05*vsw)*p.g
               +(.114*v + .886*vsu - .203*vsw)*p.b;
            p.SetClamp(r, g, b);
            GetPixel(x, y) = p;
        }
    }
}

Image* Image::Crop(int x, int y, int w, int h)
{
	/* WORK HERE */
    int start_x = fmax(0,x), end_x = fmin(width,x+w);
    int start_y = fmax(0,y), end_y = fmin(height,y+h);
    if (end_x<=start_x||end_y<=start_y)
        return NULL;
    Image* img2 = new Image(end_x-start_x,end_y-start_y);
    for(int i=0;i<img2->width;i++)
    {
        for(int j=0;j<img2->height;j++)
        {
            Pixel p = GetPixel(start_x+i, start_y+j);
            img2->SetPixel(i, j, p);
        }
    }
    return img2;
}


void Image::ExtractChannel(int channel)
{
    // 1 for r, 2 for g, 3 for b
	/* WORK HERE */
    int x,y;
    for (x = 0 ; x < Width() ; x++)
    {
        for (y = 0 ; y < Height() ; y++)
        {
            Pixel p = GetPixel(x, y);
            switch (channel) {
                case 1:
                    p.Set(p.r, 0, 0);
                    break;
                case 2:
                    p.Set(0, p.g, 0);
                    break;
                case 3:
                    p.Set(0, 0, p.b);
                    break;
                default:
                    break;
            }
            
            GetPixel(x,y) = p;
        }
    }
}


void Image::Quantize (int nbits)
{
	/* WORK HERE */
    int nr,ng,nb;
    int *color_list = new int[pow(2,nbits)];
    for(int i=0;i<sizeof(color_list);i++)
        color_list[i] = i*pow(2.0,8-nbits)+pow(2.0,8-nbits-1);
    int min_dist = 255,index_min = 0;
    
    for(int i=0;i<width;i++)
        for(int j = 0;j<height;j++)
            {
                Pixel p = GetPixel(i, j);
                min_dist = 255;index_min = -1;
                for(int k=0;k<sizeof(color_list);k++)
                    if(min_dist>abs(color_list[k]-p.r))
                    {
                        min_dist=abs(color_list[k]-p.r);
                        index_min = k;
                    }
                nr = color_list[index_min];
                min_dist = 255;index_min = -1;
                for(int k=0;k<sizeof(color_list);k++)
                    if(min_dist>abs(color_list[k]-p.g))
                    {
                        min_dist=abs(color_list[k]-p.g);
                        index_min = k;
                    }
                ng = color_list[index_min];
                min_dist = 255;index_min = -1;
                for(int k=0;k<sizeof(color_list);k++)
                    if(min_dist>abs(color_list[k]-p.b))
                    {
                        min_dist=abs(color_list[k]-p.b);
                        index_min = k;
                    }
                nb = color_list[index_min];
                /*nr =pow(2.0,8-nbits)*int((p.r+1) /pow(2.0,8-nbits))+pow(2.0,8-nbits-1);
                ng =pow(2.0,8-nbits)*int((p.g+1) /pow(2.0,8-nbits))+pow(2.0,8-nbits-1);
                nb =pow(2.0,8-nbits)*int((p.b+1) /pow(2.0,8-nbits))+pow(2.0,8-nbits-1);*/
                //printf("%d %d %dVS%d %d %d\n",p.r,p.g,p.b,nr,ng,nb);
                p.SetClamp(nr, ng, nb);
                GetPixel(i, j) = p;
            }
    

}

void Image::RandomDither (int nbits)
{
	/* WORK HERE */
    int rand_range = pow(2.0,8-nbits)/2-1;
    int nr,ng,nb;
    int *color_list = new int[pow(2,nbits)];
    for(int i=0;i<sizeof(color_list);i++)
        color_list[i] = i*pow(2.0,8-nbits)+pow(2.0,8-nbits-1);
    int min_dist = 255,index_min = 0;
    
    for(int i=0;i<width;i++)
        for(int j = 0;j<height;j++)
        {
            Pixel p = GetPixel(i, j);
            p.SetClamp(p.r+pow(-1.0,(rand()%1))*(rand()%rand_range),
                       p.g+pow(-1.0,(rand()%1))*(rand()%rand_range),
                       p.b+pow(-1.0,(rand()%1))*(rand()%rand_range));
            min_dist = 255;index_min = -1;
            for(int k=0;k<sizeof(color_list);k++)
                if(min_dist>abs(color_list[k]-p.r))
                {
                    min_dist=abs(color_list[k]-p.r);
                    index_min = k;
                }
            nr = color_list[index_min];
            min_dist = 255;index_min = -1;
            for(int k=0;k<sizeof(color_list);k++)
                if(min_dist>abs(color_list[k]-p.g))
                {
                    min_dist=abs(color_list[k]-p.g);
                    index_min = k;
                }
            ng = color_list[index_min];
            min_dist = 255;index_min = -1;
            for(int k=0;k<sizeof(color_list);k++)
                if(min_dist>abs(color_list[k]-p.b))
                {
                    min_dist=abs(color_list[k]-p.b);
                    index_min = k;
                }
            nb = color_list[index_min];
            /*nr =pow(2.0,8-nbits)*int((p.r+1) /pow(2.0,8-nbits))+pow(2.0,8-nbits-1);
             ng =pow(2.0,8-nbits)*int((p.g+1) /pow(2.0,8-nbits))+pow(2.0,8-nbits-1);
             nb =pow(2.0,8-nbits)*int((p.b+1) /pow(2.0,8-nbits))+pow(2.0,8-nbits-1);*/
            //printf("%d %d %dVS%d %d %d\n",p.r,p.g,p.b,nr,ng,nb);
            p.SetClamp(nr, ng, nb);
            GetPixel(i, j) = p;
        }

}


static int Bayer4[4][4] =
{
    {15,  7, 13,  5},
    { 3, 11,  1,  9},
    {12,  4, 14,  6},
    { 0,  8,  2, 10}
};


void Image::OrderedDither(int nbits)
{
    float ** calcu_M_n(int n);
    int matrix_size = fmax(4.0,pow(2.0,8-nbits)/2);
    float **M_n = calcu_M_n(matrix_size);
    
    /* WORK HERE */
    int nr,ng,nb;
    int *color_list = new int[pow(2,nbits)];
    for(int i=0;i<sizeof(color_list);i++)
        color_list[i] = i*pow(2.0,8-nbits)+pow(2.0,8-nbits-1);
    int min_dist = 255,index_min = 0;
    
    for(int i=0;i<width;i++)
        for(int j = 0;j<height;j++)
        {
            Pixel p = GetPixel(i, j);
            
            int ind_x = i%matrix_size;
            int ind_y = j%matrix_size;
            double tr,tg,tb;
            tr = p.r+256.0/pow(2.0,nbits)*(M_n[ind_x][ind_y]-0.5);
            tg = p.g+256.0/pow(2.0,nbits)*(M_n[ind_x][ind_y]-0.5);
            tb = p.b+256.0/pow(2.0,nbits)*(M_n[ind_x][ind_y]-0.5);
            
            p.SetClamp(tr,tg,tb);
            min_dist = 255;index_min = -1;
            for(int k=0;k<sizeof(color_list);k++)
                if(min_dist>abs(color_list[k]-p.r))
                {
                    min_dist=abs(color_list[k]-p.r);
                    index_min = k;
                }
            nr = color_list[index_min];
            min_dist = 255;index_min = -1;
            for(int k=0;k<sizeof(color_list);k++)
                if(min_dist>abs(color_list[k]-p.g))
                {
                    min_dist=abs(color_list[k]-p.g);
                    index_min = k;
                }
            ng = color_list[index_min];
            min_dist = 255;index_min = -1;
            for(int k=0;k<sizeof(color_list);k++)
                if(min_dist>abs(color_list[k]-p.b))
                {
                    min_dist=abs(color_list[k]-p.b);
                    index_min = k;
                }
            nb = color_list[index_min];
            /*nr =pow(2.0,8-nbits)*int((p.r+1) /pow(2.0,8-nbits))+pow(2.0,8-nbits-1);
             ng =pow(2.0,8-nbits)*int((p.g+1) /pow(2.0,8-nbits))+pow(2.0,8-nbits-1);
             nb =pow(2.0,8-nbits)*int((p.b+1) /pow(2.0,8-nbits))+pow(2.0,8-nbits-1);*/
            //printf("%d %d %dVS%d %d %d\n",p.r,p.g,p.b,nr,ng,nb);
            p.SetClamp(nr, ng, nb);
            GetPixel(i, j) = p;
        }

}
float ** calcu_M_n(int matrix_size)
{//return a n*n matirx
    float **M_nd2,**M_n;
    M_n = new float*[matrix_size];
    for(int i=0;i<matrix_size;i++)
        M_n[i] = new float[matrix_size];
    if (matrix_size==2)
    {
        M_n[0][0] = 0/4;M_n[0][1] = 2/4.0;
        M_n[1][0] = 3/4.0;M_n[1][1] = 1/4.0;
        return M_n;
    }
    M_nd2 = calcu_M_n(matrix_size/2);
    for(int i=0;i<matrix_size;i++)
        for(int j=0;j<matrix_size;j++)
        {
            if (i<matrix_size/2)
            {
                if(j<matrix_size/2){
                    M_n[i][j] =  M_nd2[i][j];
                }else{
                    M_n[i][j] =  M_nd2[i][j-matrix_size/2]+2.0/(matrix_size*matrix_size);
                }
            }else{
                if(j<matrix_size/2){
                    M_n[i][j] = M_nd2[i-matrix_size/2][j]+3.0/(matrix_size*matrix_size);
                }else{
                    M_n[i][j] = M_nd2[i-matrix_size/2][j-matrix_size/2]+1.0/(matrix_size*matrix_size);
                }
            }
        }
    delete M_nd2;
    return M_n;
}
/* Error-diffusion parameters */
const double
    ALPHA = 7.0 / 16.0,
    BETA  = 3.0 / 16.0,
    GAMMA = 5.0 / 16.0,
    DELTA = 1.0 / 16.0;

void Image::FloydSteinbergDither(int nbits)
{
	/* WORK HERE */
    Pixel oldpixel,newpixel;
    int r,g,b,e_r,e_g,e_b;
    int closest_value(int value,int nbits);
    for (int x=0;x<width;x++)
        for (int y=0;y<height;y++)
        {
            oldpixel  = GetPixel(x, y);
            r  = closest_value(oldpixel.r,nbits);
            g  = closest_value(oldpixel.g,nbits);
            b  = closest_value(oldpixel.b,nbits);
            newpixel.SetClamp(r, g, b);
            GetPixel(x, y) = newpixel;
            //quant_error  := oldpixel - newpixel
            e_r = oldpixel.r-newpixel.r;
            e_g = oldpixel.g-newpixel.g;
            e_b = oldpixel.b-newpixel.b;
            if (x+1<width)
            {//pixel[x + 1][y    ] := pixel[x + 1][y    ] + quant_error * 7 / 16
                oldpixel  = GetPixel(x+1, y);
                oldpixel.SetClamp(oldpixel.r+e_r*7.0/16, oldpixel.g+e_g*7.0/16, oldpixel.b+e_b*7.0/16);
                GetPixel(x+1, y) = oldpixel;
            }
            if (x-1>=0&&y+1<height)
            {//pixel[x - 1][y + 1] := pixel[x - 1][y + 1] + quant_error * 3 / 16
                oldpixel  = GetPixel(x-1, y+1);
                oldpixel.SetClamp(oldpixel.r+e_r*3.0/16, oldpixel.g+e_g*3.0/16, oldpixel.b+e_b*3.0/16);
                GetPixel(x-1, y+1) = oldpixel;
            }
            if (y+1<height)
            {//pixel[x    ][y + 1] := pixel[x    ][y + 1] + quant_error * 5 / 16
                oldpixel  = GetPixel(x, y+1);
                oldpixel.SetClamp(oldpixel.r+e_r*5.0/16, oldpixel.g+e_g*5.0/16, oldpixel.b+e_b*5.0/16);
                GetPixel(x, y+1) = oldpixel;
            }
            if (x+1<width&&y+1<height)
            {//pixel[x + 1][y + 1] := pixel[x + 1][y + 1] + quant_error * 1 / 16
                oldpixel  = GetPixel(x+1, y+1);
                oldpixel.SetClamp(oldpixel.r+e_r*1.0/16, oldpixel.g+e_g*1.0/16, oldpixel.b+e_b*1.0/16);
                GetPixel(x+1, y+1) = oldpixel;
            }
            
        }
}
int closest_value(int value,int nbits){
    return round(value / (pow(2.0,8-nbits)-1))*pow(2.0,8-nbits)+pow(2.0,8-nbits-1);
}
void Image::Blur(int n)
{
	/* WORK HERE */
    if (n%2==0)
    {
        printf("input kernel size is even;Formolized to odd by add 1\n");
    }
    double sigma = 1.0;
    double *kernel,sumk=0;
    kernel = new double[n];
    for (int i=0;i<n;i++)
    {
        kernel[i] =  (1/(sigma*sqrt(2*M_PI)))  *  pow(M_E,-.5*pow((i-(n-1)/2)/sigma,2));
        //printf("%f \n",kernel[i]);
        sumk+=kernel[i];
    }
    //for (int i=0;i<n;i++)
      //  kernel[i] =kernel[i] / sumk;
    Image* img2 = new Image(width,height);
    Pixel p;
    float r,g,b;
    for(int x=0;x<width;x++)
    {
        for (int y=0;y<height;y++)
        {
            r=0;g=0;b=0;
            for(int i=0;i<n;i++)
            {
                if (x+(i-(n-1)/2)>=0&&x+(i-(n-1)/2)<width)
                    p = GetPixel(x+(i-(n-1)/2), y);
                else
                    p = GetPixel(x-(i-(n-1)/2), y);
                r+=p.r* kernel[i];g+=p.g* kernel[i];b+=p.b* kernel[i];
            }
            p.SetClamp(r, g, b);
            img2->SetPixel(x, y, p);
        }
    }
    for(int x=0;x<width;x++)
    {
        for (int y=0;y<height;y++)
        {
            r=0;g=0;b=0;
            for(int i=0;i<n;i++)
            {
                if (y+(i-(n-1)/2)>=0&&y+(i-(n-1)/2)<height)
                    p = img2->GetPixel(x, y+(i-(n-1)/2));
                else
                    p = img2->GetPixel(x, y-(i-(n-1)/2));
                r+=p.r* kernel[i];g+=p.g* kernel[i];b+=p.b* kernel[i];
            }
            p.SetClamp(r, g, b);
            GetPixel(x, y) = p;
        }
    }
   /* for(int x=0;x<width;x++)
            for (int y=0;y<height;y++)
                GetPixel(x, y) = img->GetPixel(x, y);*/
    delete img2;
    delete [] kernel;
}

void Image::Sharpen(int n)
{
	/* WORK HERE */
    Image img2 = Image(*this);
    Pixel p1,p2,p3;
    int r,g,b;
    float alpha = 2;
    for(int x=0;x<width;x++)
        for(int y=0;y<height;y++)
            img2.SetPixel(x, y, GetPixel(x, y));
    this->Blur(n);
    for(int x=0;x<width;x++)
        for(int y=0;y<height;y++)
        {
            p1 = img2.GetPixel(x, y);
            p2 = this->GetPixel(x, y);
            r = ((1+alpha)*p1.r)-alpha*p2.r;
            g = ((1+alpha)*p1.g)-alpha*p2.g;
            b = ((1+alpha)*p1.b)-alpha*p2.b;
            p3.SetClamp(r, g, b);
            this->SetPixel(x, y, p3);
            
        }
}

void Image::EdgeDetect()
{
	/* WORK HERE */
    Pixel p;
    float r,g,b;
    int n = 2;
    int index_x,index_y;
    int kernel_1[3][3] = {{ 1, 0,-1},
                          { 0, 0, 0},
                            {-1, 0, 1}};
    int kernel_2[3][3] = {{ 0,-1, 0},
                          {-1, 4,-1},
                            {0 ,-1, 0}};
    int kernel_3[3][3] = {{-1,-1,-1},
                          {-1,8,-1},
                            {-1,-1,-1}};
    int kernel[3][3];
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            kernel[i][j] = kernel_3[i][j];
    
    Image* img2 = new Image(width,height);
    for(int x=0;x<width;x++)
        for(int y=0;y<height;y++)
        {
            r=0;g=0;b=0;
            for(int i=0;i<3;i++)
                for(int j=0;j<3;j++)
                {
                    index_x =x+(i-(n-1)/2);
                    if(index_x<0||index_x>=width)
                    {
                        index_x =x-(i-(n-1)/2);
                    }
                    index_y = y+(j-(n-1)/2);
                    if(index_y<0||index_y>=height)
                    {
                        index_y =y-(j-(n-1)/2);
                    }
                   p = GetPixel(index_x, index_y);
                   r+=p.r* kernel[i][j];g+=p.g* kernel[i][j];b+=p.b* kernel[i][j];
                }
            
            p.SetClamp(r, g, b);
            img2->SetPixel(x, y, p);
        }
    for(int x=0;x<width;x++)
        for(int y=0;y<height;y++)
        {
            p = img2->GetPixel(x, y);
            p.SetClamp(p.Luminance(), p.Luminance(), p.Luminance());
            GetPixel(x, y) =p;
        }
    
    delete img2;
}
void Image::Charcoal(){
    Pixel p,p2,p3;
    Image *gray_img =new Image(width,height);
    Image *blur_inverse_img =new Image(width,height);
    Image *edge_img =new Image(width,height);
    for(int x=0;x<width;x++)
        for(int y=0;y<height;y++)
        {
            p = GetPixel(x, y);
            edge_img->GetPixel(x, y) = p;
            p2.SetClamp(255-p.Luminance(), 255-p.Luminance(), 255-p.Luminance());
            p.SetClamp(p.Luminance(), p.Luminance(), p.Luminance());
            gray_img->GetPixel(x, y) =p;
            blur_inverse_img->GetPixel(x, y) = p2;
        }
    edge_img->EdgeDetect();
    edge_img->Blur(3);
    blur_inverse_img->Blur(5);
    for(int x=0;x<width;x++)
        for(int y=0;y<height;y++)
        {
            p = gray_img->GetPixel(x, y);
            p2 = blur_inverse_img->GetPixel(x, y);
            p3 = edge_img->GetPixel(x, y);
            double factor = 0.3;
            int luminance = (1-factor)*((p.Luminance()/255.0)/(1-p2.Luminance()/255.0)*255.0)+ factor*(255-p3.Luminance());
            p.SetClamp(luminance, luminance, luminance);
            GetPixel(x, y) = p;
        }
}
Image* Image::Scale(double sx, double sy)
{
	/* WORK HERE */
    Image* img2 = new Image(int(sx*width),int(sy*height));
    Pixel p;
    float index_x,index_y;
    for (int i=0;i<int(sx*width);i++)
        for(int j=0;j<int(sy*height);j++)
        {
            index_x = i/sx;
            index_y = j/sy;
            p = Sample(index_x, index_y);
            img2->SetPixel(i, j, p);
        }
	return img2;
}

Image* Image::Rotate(double angle)
{
	/* WORK HERE */
    angle = angle*M_PI/180;
    double R_T[4][4]={{cos(-angle), -sin(-angle), 0, 0},
                    {sin(-angle), cos(-angle), 0, 0},
                    {0,             0,          1, 0},
                    {0,             0,          0, 1}};
    int square_len = sqrt(width*width+height*height);
    Image* img2 = new Image(square_len,square_len);
    Pixel p;
    float index_x,index_y;
    for (int i=0;i<square_len;i++)
        for(int j=0;j<square_len;j++)
        {
            index_x = (i-square_len/2)*R_T[0][0]+(j-square_len/2)*R_T[0][1]+width/2;
            index_y = (i-square_len/2)*R_T[1][0]+(j-square_len/2)*R_T[1][1]+height/2;
            p = Sample(index_x, index_y);
            img2->SetPixel(i, j, p);
        }
    return img2;

}

void Image::Fun()
{
	/* WORK HERE */
    double angle = 0,r,angle_origin;
    double calcu_angle(int x,int y);
    
    double r_range = (sqrt(pow(width,2)+pow(height,2))/2)/2;
    Image* img2 = new Image(width,height);
    Pixel p;
    
    float index_x,index_y;
    for(int i=0;i<width;i++)
        for(int j=0;j<height;j++)
        {
            angle = calcu_angle(i-width/2, j-height/2);
            //printf("%d %d %f\n",i-width/2, j-height/2,calcu_angle(i-width/2, j-height/2)/M_PI*180);
            if(i-width/2==0&&j-height/2==0)
            {
                index_x = i;
                index_y = j;
            }else{
                r = sqrt((i-width/2)*(i-width/2)+(j-height/2)*(j-height/2));
                
                //r = pow(r/r_range,1.5)*r_range/fmin(1.0/abs(sin(angle)),1.0/abs(cos(angle)));
                angle_origin = atan(pow(r/r_range,1.0/2))/(r/r_range);
                index_x = (width/fmin(width,height))*angle_origin*(i-width/2)+width/2;
                index_y = (height/fmin(width,height))*angle_origin*(j-height/2)+height/2;
            }
            p = Sample(index_x, index_y);
            img2->SetPixel(i, j, p);
        }
    
    for(int i=0;i<width;i++)
        for(int j=0;j<height;j++)
        {
            p = img2->GetPixel(i, j);
            GetPixel(i, j) = p;
        }
    delete img2;
}
double calcu_angle(int r_x,int r_y)
{
    // rotation should be in relative coordinate
    double g_angle;
    //Compute the new angle, g_angle, based on the mouse positions
    double r;
    r = sqrt(r_x*r_x+r_y*r_y);
    
    g_angle = acos(r_x/r);
    if(r_y<0)
        g_angle = -g_angle;
    while(g_angle>2*M_PI)
        g_angle = g_angle-2*M_PI;
    while(g_angle<0)
        g_angle = g_angle+2*M_PI;

    return g_angle;
}
/**
 * Image Sample
 **/
void Image::SetSamplingMethod(int method)
{
    assert((method >= 0) && (method < IMAGE_N_SAMPLING_METHODS));
    sampling_method = method;
}


Pixel Image::Sample (double u, double v){
    /* WORK HERE */
    if (u<0||v<0||u>=width||v>=height)
        return Pixel();
    if (sampling_method==0)//point sampling
    {
        int x = round(u);
        int y = round(v);
        if (x<width&&y<height&&x>=0&&y>=0)
            return GetPixel(x, y);
        else
            return Pixel();
    }
    else if(sampling_method==1){//bilinear sampling
        
        int x_f = floor(u),y_f = floor(v);
        int x_c = ceil(u),y_c = ceil(v);
        double a = u-x_f,b = v-y_f;
        x_f = x_f<0?0:x_f;
        x_f = x_f>=width?width-1:x_f;
        x_c = x_c<0?0:x_c;
        x_c = x_c>=width?width-1:x_c;
        y_f = y_f<0?0:y_f;
        y_f = y_f>=height?height-1:y_f;
        y_c = y_c<0?0:y_c;
        y_c = y_c>=height?height-1:y_c;
        
        
        
        
        Pixel p1 = GetPixel(x_f, y_f);
        Pixel p2 = GetPixel(x_f, y_c);
        Pixel p3 = GetPixel(x_c, y_f);
        Pixel p4 = GetPixel(x_c, y_c);
        Pixel p;
        p.SetClamp((1-a)*(1-b)*p1.r + a*(1-b)*p3.r + a*b*p4.r + (1-a)*b*p2.r
                   , (1-a)*(1-b)*p1.g + a*(1-b)*p3.g + a*b*p4.g + (1-a)*b*p2.g
                   , (1-a)*(1-b)*p1.b + a*(1-b)*p3.b + a*b*p4.b + (1-a)*b*p2.b);
        
        return p;
    }
    else if(sampling_method==2){//Gaussian sampling
        int n=5;
        if (n%2==0)
        {
            printf("input kernel size is even;Formolized to odd by add 1\n");
        }
        double sigma = 1.0;
        double **kernel;
        kernel = new double *[n];
        for(int i=0;i<n;i++)
        {
            kernel[i] = new double [n];
            for(int j=0;j<n;j++)
            {
                kernel[i][j] =  (1/(pow(sigma,2.0)*2*M_PI))
                *  pow(M_E,-.5*(pow((i-(n-1)/2),2)+pow((j-(n-1)/2),2))/pow(sigma,2.0));
                //printf("%f \n",kernel[i]);
            }
        }
        
        //for (int i=0;i<n;i++)
        //  kernel[i] =kernel[i] / sumk;
        Pixel p;
        float r,g,b;
        int x = u;
        int y = v;
        r=0;g=0;b=0;
        for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
        {
            int ind_x,ind_y;
            if (x+(i-(n-1)/2)>=0&&x+(i-(n-1)/2)<width)
                ind_x = x+(i-(n-1)/2);
            else
                ind_x = x-(i-(n-1)/2);
            if (y+(i-(n-1)/2)>=0&&y+(i-(n-1)/2)<height)
                ind_y = y+(i-(n-1)/2);
            else
                ind_y = y-(i-(n-1)/2);
            p = GetPixel(ind_x, ind_y);
            r+=p.r* kernel[i][j];g+=p.g* kernel[i][j];b+=p.b* kernel[i][j];
        }
        p.SetClamp(r, g, b);
        return p;
    }
    else if(sampling_method==3){//uniform supersampling
        
        return Pixel();
        }
    else
        return Pixel();
}

