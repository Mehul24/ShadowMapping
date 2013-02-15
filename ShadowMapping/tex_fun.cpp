/* Texture functions for cs580 GzLib	*/
#include    "stdafx.h" 
#include    "stdio.h"
#include    "Gz.h"
#include    "math.h"

GzColor*    image;
int         xs, ys;
int         reset   = 1;

/* Image texture function */
// TODO: Document
int tex_fun(float u, float v, GzColor color)
{
    unsigned char pixel[3];
    unsigned char dummy;
    char          foo[8];
    int           i, j;
    FILE*         fd;

    // Open and load texture file
    if (reset)
    {
        fd = fopen ("texture", "rb");
        if(fd == NULL)
        {
            fprintf (stderr, "texture file not found\n");
            exit(-1);
        }

        fscanf (fd, "%s %d %d %c", foo, &xs, &ys, &dummy);
        image = (GzColor*)malloc(sizeof(GzColor)*(xs+1)*(ys+1));

        if (image == NULL)
        {
            fprintf (stderr, "malloc for texture image failed\n");
            exit(-1);
        }

        // Create an array of GzColor values
        for (i = 0; i < xs*ys; i++)
        {
            fread(pixel, sizeof(pixel), 1, fd);
            image[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
            image[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
            image[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
        }

        // Initialization is done
        reset = 0;
        fclose(fd);
    }

    memset(color, 0x00, sizeof(GzColor));

    u = u < 0.0f ? 0.0f : u;
    u = u > 1.0f ? 1.0f : u;
    v = v < 0.0f ? 0.0f : v;
    v = v > 1.0f ? 1.0f : v;

    float uTranslate = u * (xs - 1);
    float vTranslate = v * (ys - 1);

    unsigned int uLeft = 0, uRight     = 0;
    unsigned int vTop  = 0, vBottom    = 0;

    float du = 0, dv = 0;

    uLeft   = floor(uTranslate);
    uRight  = ceil(uTranslate);
    vTop    = floor(vTranslate);
    vBottom = ceil(vTranslate);

    du = uTranslate - uLeft;
    dv = vTranslate - vTop;

    unsigned int topLeftColorIndex      = (xs) * vTop + uLeft;
    unsigned int topRightColorIndex     = (xs) * vTop + uRight;
    unsigned int bottomLeftColorIndex   = (xs) * vBottom + uLeft;
    unsigned int bottomRightColorIndex  = (xs) * vBottom + uRight;

    for(unsigned int colorIndex = 0; colorIndex < 3; ++colorIndex)
    {
        float topLeftColor      = image[topLeftColorIndex][colorIndex];
        float topRightColor     = image[topRightColorIndex][colorIndex];
        float bottomLeftColor   = image[bottomLeftColorIndex][colorIndex];
        float bottomRightColor  = image[bottomRightColorIndex][colorIndex];

        color[colorIndex] =
            (1 - du) * (1 - dv) * topLeftColor + du * (1 - dv) * topRightColor +
            (1 - du) * dv * bottomLeftColor + du * dv * bottomRightColor;
    }

    return GZ_SUCCESS;
}

/* Procedural texture function */
//TODO: Document
int ptex_fun(float u, float v, GzColor color)
{
    memset(color, 0x00, sizeof(GzColor));

    // Boundary condition
    u = u < 0.0f ? 0.0f : u;
    u = u > 1.0f ? 1.0f : u;
    v = v < 0.0f ? 0.0f : v;
    v = v > 1.0f ? 1.0f : v;

    // This value determines the number of boxes. Increasing it reduces the number of boxes
    const float scaleFactor = 0.2f;

    unsigned int uValue = (unsigned int) (u / scaleFactor);
    unsigned int vValue = (unsigned int) (v / scaleFactor);

    if(uValue % 2 == 0 && vValue % 2 == 0 || uValue % 2 == 1 && vValue % 2 == 1)
    {
        color[0] = 1.0f;
        color[1] = 1.0f;
    }

    return GZ_SUCCESS;
}
