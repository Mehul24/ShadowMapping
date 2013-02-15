/*   CS580 HW   */
#include    "stdafx.h"  
#include	"Gz.h"
#include	"disp.h"

int GzNewFrameBuffer(char** framebuffer, int width, int height)
{
/* create a framebuffer:
 -- allocate memory for framebuffer : (sizeof)GzPixel x width x height
 -- pass back pointer 
*/
	int size = 0;
	size = sizeof(GzPixel) * width * height;
	*framebuffer = new char[size];
	return GZ_SUCCESS;
}

int GzNewDisplay(GzDisplay	**display, GzDisplayClass dispClass, int xRes, int yRes)
{

/* create a display:
  -- allocate memory for indicated class and resolution
  -- pass back pointer to GzDisplay object in display
*/
	*display = new GzDisplay() ;
	((*display)->xres) = xRes;
	(*display)->yres = yRes;
	(*display)->dispClass = dispClass;
	(*display)->fbuf = new GzPixel[sizeof(GzPixel) * xRes * yRes];
	return GZ_SUCCESS;
}


int GzFreeDisplay(GzDisplay	*display)
{
	/* clean up, free memory */	
	free(display);
	return GZ_SUCCESS;
}


int GzGetDisplayParams(GzDisplay *display, int *xRes, int *yRes, GzDisplayClass	*dispClass)
{
/* pass back values for an open display */
	*xRes = display->xres;
	*yRes = display->yres;
	*dispClass = display->dispClass;
	return GZ_SUCCESS;
}


int GzInitDisplay(GzDisplay	*display)
{
/* set everything to some default values - start a new frame */

	int i,j;
	for(i=0; i <display->xres;i++)
	{
		for(j=0;j<display->yres ;j++)
		{
			
			display->fbuf[ARRAY(i,j)].red = 0x1FF;
			display->fbuf[ARRAY(i,j)].blue = 0x1FF;
			display->fbuf[ARRAY(i,j)].green = 0x1FF;
			display->fbuf[ARRAY(i,j)].alpha = 1;
			display->fbuf[ARRAY(i,j)].z = MAXINT;
		}
	}

	return GZ_SUCCESS;
}


int GzPutDisplay(GzDisplay *display, int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
/* write pixel values into the display */
	int index;
	if(i>256)
		i = 256;
	if(i<0)
		i = 0;
	if(j>256)
		j = 256;
	if(j<0)
		j = 0;
	if(r>4095)
		r = 4095;
	if(r<0)
		r = 0;
	if(g>4095)
		g = 4095;
	if(g<0)
		g = 0;
	if(b>4095)
		b = 4095;
	if(b<0)
		b = 0;
	
	index= ARRAY(i,j);
	display->fbuf[index].red = r;
	display->fbuf[index].blue = b;
	display->fbuf[index].green = g;
	display->fbuf[index].alpha = a;
	display->fbuf[index].z = z;
//	printf("(%d,%d)\n", i, j);
	return GZ_SUCCESS;
}


int GzGetDisplay(GzDisplay *display, int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth *z)
{
	/* pass back pixel value in the display */
	/* check display class to see what vars are valid */
	int index;
	if(i>256)
		i = 256;
	if(i<0)
		i = 0;
	if(j>256)
		j = 256;
	if(j<0)
		j = 0;
	if(*r>4095)
		*r = 4095;
	if(*r<0)
		*r = 0;
	if(*g>4095)
		*g = 4095;
	if(*g<0)
		*g = 0;
	if(*b>4095)
		*b = 4095;
	if(*b<0)
		*b = 0;

	index= ARRAY(i,j);
	*r = display->fbuf[index].red;
	*b = display->fbuf[index].blue;
	*g = display->fbuf[index].green;
	*a = display->fbuf[index].alpha;
	*z = display->fbuf[index].z;

	return GZ_SUCCESS;
}


int GzFlushDisplay2File(FILE* outfile, GzDisplay *display)
{

	/* write pixels to ppm file based on display class -- "P6 %d %d 255\r" */
	int size,i,j,k;
	char rgb[3];
	fprintf(outfile,"P6 %d %d 255\n",256,256);
	for(j=0;j<display->yres;j++)
	{
		for(i=0;i<display->xres;i++)
		{
			k=ARRAY(i,j);
			rgb[0]=(char)(display->fbuf[k].red>>4) & 0xff;
			rgb[1]=(char)(display->fbuf[k].green>>4) & 0xff;
			rgb[2]=(char)(display->fbuf[k].blue>>4) & 0xff;
			fwrite(rgb,sizeof(rgb),1,outfile);
		}
	}

	return GZ_SUCCESS;
}

int GzFlushDisplay2FrameBuffer(char* framebuffer, GzDisplay *display)
{

	/* write pixels to framebuffer: 
		- Put the pixels into the frame buffer
		- Caution: store the pixel to the frame buffer as the order of blue, green, and red 
		- Not red, green, and blue !!!
	*/
	
	int size,i,j,k,cnt=0;
	for(j=0;j<display->yres;j++)
	{
		for(i=0;i<display->xres;i++)
		{
			k=ARRAY(i,j);
			framebuffer[cnt++]=(display->fbuf[k].blue>>4) & 0xff;
			framebuffer[cnt++]=(display->fbuf[k].green>>4) & 0xff;
			framebuffer[cnt++]=(display->fbuf[k].red>>4) & 0xff;
		}
	}

	return GZ_SUCCESS;
}