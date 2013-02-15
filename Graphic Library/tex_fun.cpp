/* Texture functions for cs580 GzLib	*/
#include    "stdafx.h" 
#include	"stdio.h"
#include	"Gz.h"
#include "rend.h"
#include<math.h>

#define	FINDCOLOR(x,y)	(x+(y*xs))


GzColor	*image;
int xs, ys;
int reset = 1;


/* Image texture function */
int tex_fun(float u, float v, GzColor color)
{
  unsigned char		pixel[3];
  unsigned char     dummy;
  char  		foo[8];
  int   		i, j;
  FILE			*fd;

  if (reset) {          /* open and load texture file */
    fd = fopen ("texture", "rb");
    if (fd == NULL) {
      fprintf (stderr, "texture file not found\n");
      exit(-1);
    }
    fscanf (fd, "%s %d %d %c", foo, &xs, &ys, &dummy);
    image = (GzColor*)malloc(sizeof(GzColor)*(xs+1)*(ys+1));
    if (image == NULL) {
      fprintf (stderr, "malloc for texture image failed\n");
      exit(-1);
    }

    for (i = 0; i < xs*ys; i++) {	/* create array of GzColor values */
      fread(pixel, sizeof(pixel), 1, fd);
      image[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
      image[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
      image[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
      }

    reset = 0;          /* init is done */
	fclose(fd);
  }

/* bounds-test u,v to make sure nothing will overflow image array bounds */
/* determine texture cell corner values and perform bilinear interpolation */
/* set color to interpolated GzColor value and return */


  //s = floor y ; t = floor X ; ABCD
 // Color(p) = s t C + (1-s) t D + s (1-t) B + (1-s) (1-t) A  

  if(u<0)
	  u=0;
  if(v<0)
	  v=0;
  if(u>1)
	  u=1;
  if(v>1)
	  v=1;

  float s , t;
  GzColor temp;
  int ax,ay,bx,by,cx,cy,dx,dy;

  u = u*(xs - 1);
  v = v*(ys - 1);

  s = -floor(v) +v;
  t = -floor(u) + u;

  ax = floor(u);
  ay = floor(v);

  bx = floor(u);
  by = ceil(v);

  cx = ceil(u);
  cy = ceil(v);

  dx = ceil(u);
  dy = floor(v);

  for(int i = 0 ;i<3;i++)
	  color[i] =  s*t*(image[FINDCOLOR(cx,cy)][i]) + (1-s)*t*(image[FINDCOLOR(dx,dy)][i]) + s*(1-t)*(image[FINDCOLOR(bx,by)][i]) + (1-s)*(1-t)*(image[FINDCOLOR(ax,ay)][i]);
 }


/* Procedural texture function */
int ptex_fun(float u, float v, GzColor color)
{

 if(u<0)
	  u=0;
  if(v<0)
	  v=0;
  if(u>1)
	  u=1;
  if(v>1)
	  v=1;

  memcpy(color,"0x03F",sizeof(GzColor));
	
  int scaleu = u/0.1f;
  int scalev = v/0.1f;

  if((scaleu %2 ==0 && scalev %2 ==0) || (scaleu %2 ==1 && scalev %2 ==1))
  {
	  color[0] = 1;
	  color[1] = 0 ;
	  color[2] = 1;
  }

  	return 0;

}

