/*************************************************************************************************
 * CS580 Homework 1                                                                              *
 *                                                                                               *
 * Copyright (C) 2012. Karthik Nallapeta Subramanya and University of Southern California.       *
 *                                                                                               *
 * Contact:                                                                                      *
 * karthik.nallapeta@gmail.com                                                                   *
 * knsubram@usc.edu                                                                              *
 *************************************************************************************************/

#include "stdafx.h"
#include "Gz.h"
#include "disp.h"

/**
 * Creates a frameBuffer for the display by allocating (sizeof)GzPixel * width * height bytes.
 * <P>
 * @param frameBuffer Output parameter representing the frameBuffer of display
 * @param width The width of the frameBuffer
 * @param height The height of the frameBuffer
 * @return Operation success value
 */
int GzNewFrameBuffer(char** frameBuffer, int width, int height)
{
    *frameBuffer = (char *) malloc (sizeof(GzPixel) * width * height);

    if(*frameBuffer == NULL)
    {
        return GZ_FAILURE;
    }

    return GZ_SUCCESS;
}

/**
 * Creates and initializes a display object and allocates memory for indicated class and resolution.
 * <P>
 * @param display Output parameter representing the display object
 * @param displayClass
 * @param xResolution Resolution of the x axis
 * @param yResolution Resolution of the y axis
 * @return Operation success value
 */
int GzNewDisplay(GzDisplay **display, GzDisplayClass displayClass, int xResolution, int yResolution)
{
    if(xResolution < 0)
    {
        xResolution = 0;
    }

    if(xResolution >= MAXXRES)
    {
        xResolution = MAXXRES;
    }

    if(yResolution < 0)
    {
        yResolution = 0;
    }

    if(yResolution >= MAXYRES)
    {
        yResolution = MAXYRES;
    }

    *display = (GzDisplay*) malloc (sizeof(GzDisplay));

    if(*display == NULL)
    {
        return GZ_FAILURE;
    }

    (*display)->fbuf = (GzPixel*) malloc (sizeof(GzPixel) * xResolution * yResolution);

    if((*display)->fbuf == NULL)
    {
        free(*display);
        *display = NULL;

        return GZ_FAILURE;
    }

    (*display)->xres        = xResolution - 1;
    (*display)->yres        = yResolution - 1;
    (*display)->open        = 1;
    (*display)->dispClass   = displayClass;

    return GZ_SUCCESS;
}

/**
 * Frees memory of the display object and its children.
 * <P>
 * @param display The display object
 * @return Operation success value
 */
int GzFreeDisplay(GzDisplay *display)
{
    if(display != NULL)
    {
        free(display->fbuf);
        display->fbuf = NULL;
        free(display);
        display = NULL;
    }

    return GZ_SUCCESS;
}

/**
 * Pass back values for an open display.
 * <P>
 * @param display The display object
 * @param xResolution X axis resolution output parameter
 * @param yResolution Y axis resolution output parameter
 * @param displayClass
 * @return Operation success value
 */
int GzGetDisplayParams(GzDisplay *display, int *xResolution, int *yResolution, GzDisplayClass *displayClass)
{
    if(display->open != 1)
    {
        return GZ_FAILURE;
    }

    *xResolution    = display->xres;
    *yResolution    = display->yres;
    *displayClass   = display->dispClass;

    return GZ_SUCCESS;
}

/**
 * Set the display to default values; start a new frame.
 * <P>
 * @param display The display object
 * @return Operation success value
 */
int GzInitDisplay(GzDisplay *display)
{
    if(display->open != 1)
    {
        return GZ_FAILURE;
    }

    int numberOfPixels = (display->xres + 1) * (display->yres + 1);

    for(unsigned int index = 0; index < numberOfPixels; ++index)
    {
        //Initialize to the color 0x807060
        display->fbuf[index].red    = 0x800;
        display->fbuf[index].green  = 0x700;
        display->fbuf[index].blue   = 0x600;
        display->fbuf[index].alpha  = 1;

        memset(&display->fbuf[index].z, 0xFF, sizeof(display->fbuf[index].z));
        display->fbuf[index].z = MAXINT;
    }

    return GZ_SUCCESS;
}

/**
 * Write pixel values into the display.
 * <P>
 * @param display The display object
 * @param xPosition The X position of the pixel
 * @param yPosition The Y position of the pixel
 * @param r The red intensity of the pixel
 * @param g The green intensity of the pixel
 * @param b The blue intensity of the pixel
 * @param a The alpha intensity of the pixel
 * @param z The Z location of the pixel
 * @return Operation success value
 */
int GzPutDisplay(GzDisplay *display, int xPosition, int yPosition, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
    if(display->open != 1)
    {
        return GZ_FAILURE;
    }

    if(xPosition < 0 || yPosition < 0 || xPosition > display->xres || yPosition > display->yres)
    {
        return GZ_SUCCESS;
    }

    if(r > 4095)
    {
        r = 4095;
    }

    if(g > 4095)
    {
        g = 4095;
    }

    if(b > 4095)
    {
        b = 4095;
    }

    const unsigned int  xResolution     = display->xres + 1;
    unsigned int        absoluteIndex   = yPosition * xResolution + xPosition;

    display->fbuf[absoluteIndex].red    = r;
    display->fbuf[absoluteIndex].green  = g;
    display->fbuf[absoluteIndex].blue   = b;
    display->fbuf[absoluteIndex].alpha  = a;
    display->fbuf[absoluteIndex].z      = z;

    return GZ_SUCCESS;
}

/**
 * Gets back pixel value in the display for given coordinates.
 * <P>
 * @param display The display object
 * @param xPosition The X position of the pixel
 * @param yPosition The Y position of the pixel
 * @param r The red intensity of the pixel: Output parameter
 * @param g The green intensity of the pixel: Output parameter
 * @param b The blue intensity of the pixel: Output parameter
 * @param a The alpha intensity of the pixel: Output parameter
 * @param z The Z location of the pixel: Output parameter
 * @return Operation success value. GZ_FAILURE is returned when the xPosition and yPosition are not in the bounds
 */
int GzGetDisplay(GzDisplay *display, int xPosition, int yPosition, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth *z)
{
	/* pass back pixel value in the display */
	/* check display class to see what vars are valid */
    if(display->open != 1)
    {
        return GZ_FAILURE;
    }

    if(xPosition < 0 || xPosition > display->xres || yPosition < 0 || yPosition > display->yres)
    {
        return GZ_FAILURE;
    }

    const unsigned int  xResolution     = display->xres + 1;
    unsigned int        absoluteIndex   = yPosition * xResolution + xPosition;

    *r = display->fbuf[absoluteIndex].red;
    *g = display->fbuf[absoluteIndex].green;
    *b = display->fbuf[absoluteIndex].blue;
    *a = display->fbuf[absoluteIndex].alpha;
    *z = display->fbuf[absoluteIndex].z;

    return GZ_SUCCESS;
}

/**
 * Write pixels to a file in ppm format based on the contents of display class
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * PPM File Format: "P6 %d %d 255\n" followed by RGB pixel data. 8 bits for each of the RGB components.
 *
 * @param outFile The output file descriptor where image is stored
 * @param display The display object
 * @return Operation success value
 */
int GzFlushDisplay2File(FILE* outFile, GzDisplay *display)
{
    if(display->open != 1)
    {
        return GZ_FAILURE;
    }

    char string[100];
    sprintf(string, "P6 %d %d 255\n", display->xres+1, display->yres+1);

    fwrite(string, sizeof(char), strlen(string), outFile);

    for(int y=0; y<=display->yres; ++y)
    {
        for(int x=0; x <= display->xres; ++x)
        {
            int index = (display->xres + 1) * y + x;

            /*
             * Every intensity is stored as a 16bit short value. Of this, the lower 12 bits are valid.
             * Right shift it by 4 bits to save it in the file.
             */
            char colors[3];
            colors[0] = display->fbuf[index].red >> 4;
            colors[1] = display->fbuf[index].green >> 4;
            colors[2] = display->fbuf[index].blue >> 4;

            fwrite(colors, sizeof(char), 3, outFile);
        }
    }

    return GZ_SUCCESS;
}

/**
 * Write pixels in the display object to the frameBuffer.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * The pixels are stored in the frame buffer in the order blue, green, and red.
 *
 * @param frameBuffer The frameBuffer to populate
 * @param display The display object
 * @return Operation success value
 */
int GzFlushDisplay2FrameBuffer(char* frameBuffer, GzDisplay *display)
{
    if(display->open != 1)
    {
        return GZ_FAILURE;
    }

    const unsigned int xResolution      = display->xres + 1;
    const unsigned int numberOfColors   = 3;

    for(int yIndex = 0; yIndex <= display->yres; ++yIndex)
    {
        for(int xIndex = 0; xIndex <= display->xres; ++xIndex)
        {
            unsigned int absoluteIndex      = yIndex * xResolution + xIndex;
            unsigned int frameBufferIndex   = 3 * absoluteIndex;

            /*
             * Every intensity is stored as a 16bit short value. Of this, the lower 12 bits are valid.
             * Right shift it by 4 bits to save it in the buffer required format.
             */
            frameBuffer[frameBufferIndex + 0] = display->fbuf[absoluteIndex].blue >> 4;
            frameBuffer[frameBufferIndex + 1] = display->fbuf[absoluteIndex].green >> 4;
            frameBuffer[frameBufferIndex + 2] = display->fbuf[absoluteIndex].red >> 4;
        }
    }

    return GZ_SUCCESS;
}
