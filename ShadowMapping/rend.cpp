#include "stdafx.h"
#include "stdio.h"
#include "math.h"
#include "Gz.h"
#include "rend.h"
#include "gzMath.h"

/**
 * Converts a property from screen space to image space.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * None.
 *
 * @param image Output property in image space
 * @param screen Input property in screen space
 * @param vzprime A value calculated based on screen z at the point where the property is present
 */
void GzConvertScreenToImage(GzCoord image, GzCoord screen, float vzprime)
{
    for(unsigned int index = 0; index < 3; ++index)
    {
        image[index] = screen[index] * (vzprime + 1.0f);
    }
}

/**
 * Converts a property from image space to screen space.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * None.
 *
 * @param screen Output property in screen space
 * @param image Input property in image space
 * @param vzprime A value calculated based on screen z at the point where the property is present
 */
void GzConvertImageToScreen(GzCoord screen, GzCoord image, float vzprime)
{
    for(unsigned int index = 0; index < 3; ++index)
    {
        if(vzprime + 1.0f == 0.0f)
        {
            vzprime = 0.0f;
        }

        screen[index] = image[index] / (vzprime + 1.0f);
    }
}

/**
 * Computes the Xpi and Xiw for given camera parameters.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * None.
 *
 * @param camera The camera parameters.
 */
void MakeCameraMatrices(GzCamera *camera)
{
    if(camera == NULL)
    {
        return;
    }

    GzMakeIdentityMatrix(camera->Xpi);
    camera->Xpi[3][2] = tan( (camera->FOV / 2) * PI / 180.0);

    GzCoord IX = {0};
    GzCoord IY = {0};
    GzCoord IZ = {0};

    // Z = cl / ||cl||
    CreateVector(IZ, camera->lookat, camera->position);
    NormalizeVector(IZ);

    // Y, Up Prime
    CopyVector(IY, IZ);
    float magnitude = VectorDotProduct(camera->worldup, IZ);
    MultiplyScalarToVector(magnitude, IY);
    CreateVector(IY, camera->worldup, IY);
    NormalizeVector(IY);

    VectorCrossProduct(IX, IY, IZ);

    float translationX = - VectorDotProduct(IX, camera->position);
    float translationY = - VectorDotProduct(IY, camera->position);
    float translationZ = - VectorDotProduct(IZ, camera->position);

    camera->Xiw[0][0] = IX[0];
    camera->Xiw[0][1] = IX[1];
    camera->Xiw[0][2] = IX[2];
    camera->Xiw[0][3] = translationX;

    camera->Xiw[1][0] = IY[0];
    camera->Xiw[1][1] = IY[1];
    camera->Xiw[1][2] = IY[2];
    camera->Xiw[1][3] = translationY;

    camera->Xiw[2][0] = IZ[0];
    camera->Xiw[2][1] = IZ[1];
    camera->Xiw[2][2] = IZ[2];
    camera->Xiw[2][3] = translationZ;

    camera->Xiw[3][0] = 0.0f;
    camera->Xiw[3][1] = 0.0f;
    camera->Xiw[3][2] = 0.0f;
    camera->Xiw[3][3] = 1.0f;

    // Inverse camera transform, i.e. image to world space
    camera->Xwi[0][0] = IX[0];
    camera->Xwi[1][0] = IX[1];
    camera->Xwi[2][0] = IX[2];
    camera->Xwi[0][3] = camera->position[0];

    camera->Xwi[0][1] = IY[0];
    camera->Xwi[1][1] = IY[1];
    camera->Xwi[2][1] = IY[2];
    camera->Xwi[1][3] = camera->position[1];

    camera->Xwi[0][2] = IZ[0];
    camera->Xwi[1][2] = IZ[1];
    camera->Xwi[2][2] = IZ[2];
    camera->Xwi[2][3] = camera->position[2];

    camera->Xwi[3][0] = 0.0f;
    camera->Xwi[3][1] = 0.0f;
    camera->Xwi[3][2] = 0.0f;
    camera->Xwi[3][3] = 1.0f;
}

//TODO
void MakeLightMatrices(GzLight *light)
{
    if(light == NULL)
    {
        return;
    }

    GzMakeIdentityMatrix(light->Xlw);
    light->Xlw[3][2] = tan( (light->FOV / 2) * PI / 180.0);

    GzCoord IX = {0};
    GzCoord IY = {0};
    GzCoord IZ = {0};

    // Z = cl / ||cl||
    CreateVector(IZ, light->lookat, light->position);
    NormalizeVector(IZ);

    // Y, Up Prime
    CopyVector(IY, IZ);
    float magnitude = VectorDotProduct(light->worldup, IZ);
    MultiplyScalarToVector(magnitude, IY);
    CreateVector(IY, light->worldup, IY);
    NormalizeVector(IY);

    VectorCrossProduct(IX, IY, IZ);

    float translationX = - VectorDotProduct(IX, light->position);
    float translationY = - VectorDotProduct(IY, light->position);
    float translationZ = - VectorDotProduct(IZ, light->position);

    light->Xlw[0][0] = IX[0];
    light->Xlw[0][1] = IX[1];
    light->Xlw[0][2] = IX[2];
    light->Xlw[0][3] = translationX;

    light->Xlw[1][0] = IY[0];
    light->Xlw[1][1] = IY[1];
    light->Xlw[1][2] = IY[2];
    light->Xlw[1][3] = translationY;

    light->Xlw[2][0] = IZ[0];
    light->Xlw[2][1] = IZ[1];
    light->Xlw[2][2] = IZ[2];
    light->Xlw[2][3] = translationZ;

    light->Xlw[3][0] = 0.0f;
    light->Xlw[3][1] = 0.0f;
    light->Xlw[3][2] = 0.0f;
    light->Xlw[3][3] = 1.0f;

    //compute Xpl matrix
    GzMakeIdentityMatrix(light->Xpl);
    light->Xpl[3][2] = tan( (light->FOV / 2) * PI / 180.0);
}

/**
 * Dynamically create a GzRender object
 * KEEP CLOSED UNTIL BEGINRENDER INITS ARE DONE
 * Span interpolator needs pointer to display for pixel writes
 * Check for legal class GZ_Z_BUFFER_RENDER
 * Compute the Xsp
 * Initialize the default camera
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * None.
 *
 * @param frameBuffer Output parameter representing the frameBuffer of display
 * @param width The width of the frameBuffer
 * @param height The height of the frameBuffer
 * @return Operation success value
 */
int GzNewRender(GzRender** render, GzRenderClass renderClass, GzDisplay* display)
{
    *render = (GzRender*) malloc (sizeof(GzRender));
    memset(*render, 0x00, sizeof(GzRender));

    if(renderClass != GZ_Z_BUFFER_RENDER)
    {
        (*render)->open = 0;
        return GZ_FAILURE;
    }

    (*render)->renderClass  = renderClass;
    (*render)->display      = display;
    (*render)->open         = 1;
    (*render)->matlevel     = -1;
    (*render)->matlevelsXwm = -1;

    // Setup default camera
    (*render)->camera.FOV = DEFAULT_FOV;

    for(unsigned int index = 0; index < 3; ++index)
    {
        (*render)->camera.lookat[index] = 0;
    }

    (*render)->camera.position[0] = DEFAULT_IM_X;
    (*render)->camera.position[1] = DEFAULT_IM_Y;
    (*render)->camera.position[2] = DEFAULT_IM_Z;

    (*render)->camera.worldup[0] = 0;
    (*render)->camera.worldup[1] = 1;
    (*render)->camera.worldup[2] = 0;

    MakeCameraMatrices( &((*render)->camera) );

    // Initiailize the interpolation mode to flat shading
    (*render)->interp_mode = GZ_NONE;

    // Initialize lighting parameters
    (*render)->numlights = -1;

    for(unsigned int index = 0; index < 3; ++index)
    {
        // Intialize drawing color to black
        (*render)->flatcolor[index] = 0.0f;

        // Intialize material properties
        (*render)->Ka[index] = 0.0f;
        (*render)->Kd[index] = 0.0f;
        (*render)->Ks[index] = 0.0f;

        // Intialize ambient light properties
        (*render)->ambientlight.color[index] = 0.0f;
        (*render)->ambientlight.direction[index] = 0.0f;
    }

    // Intialize specular constant
    (*render)->spec = 0;

    (*render)->aliasX = 0.0f;
    (*render)->aliasY = 0.0f;

    GzMakeIdentityMatrix((*render)->Xlbp);
    (*render)->Xlbp[0][0] = LIGHT_RES_X/2;
    (*render)->Xlbp[1][1] = -LIGHT_RES_Y/2;
    (*render)->Xlbp[0][3]= LIGHT_RES_X/2;
    (*render)->Xlbp[1][3]= LIGHT_RES_Y/2;

    float d = 1 / tan( (40 / 2) * PI / 180.0);
    (*render)->Xlbp[2][2] = INT_MAX / d;

    return GZ_SUCCESS;
}

/**
 * Free all renderer resources.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * The calling function has to set GzRender to NULL.
 *
 * @param render The render that is freed.
 * @return A GZ_SUCCESS value
 */
int GzFreeRender(GzRender *render)
{
    free(render);
    return GZ_SUCCESS;
}

/**
 * Setup for the start of each frame. Initialize the frame buffer.
 * Compute Xiw and projection xform Xpi from the camera definition.
 * Push on Xpi and Xiw on the stack.
 * Now the stack contains Xsw and app can push model Xforms.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * None.
 *
 * @param render The render object which contains the display object to be initialized.
 * @return GZ_SUCCESS or GZ_FAILURE value
 */
int GzBeginRender(GzRender *render)
{
    float d = 1 / tan( (render->camera.FOV / 2) * PI / 180.0);

    // Setup Xsp
    GzMakeIdentityMatrix(render->Xsp);

    render->Xsp[0][0] = (render->display->xres + 1) / 2;
    render->Xsp[1][1] = -(render->display->yres + 1) / 2;
    render->Xsp[0][3] = (render->display->xres + 1) / 2;
    render->Xsp[1][3] = (render->display->yres + 1) / 2;
    render->Xsp[2][2] = INT_MAX / d;

    GzPushMatrix(render, render->Xsp);

    GzPushMatrix(render, render->camera.Xpi);
    GzPushMatrix(render, render->camera.Xiw);

    if(render == NULL || render-> open != 1)
    {
        return GZ_FAILURE;
    }

    return GzInitDisplay(render->display);
}

/**
 * 
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * None.
 *
 * @param render The render object which contains the display object to be initialized.
 * @return GZ_SUCCESS or GZ_FAILURE value
 */
int GzPutCamera(GzRender *render, GzCamera *camera)
{
    render->camera.FOV = camera->FOV;

    for(unsigned int index = 0; index < 3; ++index)
    {
        render->camera.lookat[index]    = camera->lookat[index];
        render->camera.position[index]  = camera->position[index];
        render->camera.worldup[index]   = camera->worldup[index];
    }

    MakeCameraMatrices(&(render->camera));

    return GZ_SUCCESS;
}

/**
 * Adds a transformation matrix to the top of the stack.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * Multiplies the current top matrix with the new modelling matrix and pushes the result to the top of the stack.
 *
 * @param render The render object for which the Xform stack is updated.
 * @return GZ_FAILURE in case of overflow else GZ_SUCCESS.
 */
int GzPushMatrix(GzRender *render, GzMatrix matrix)
{
    return GzPushMatrix(render, matrix, false);
}

/**
 * Adds a transformation matrix to the top of the stack.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * Multiplies the current top matrix with the new modelling matrix and pushes the result to the top of the stack.
 *
 * @param render The render object for which the Xform stack is updated.
 * @param matrix Matrix multiplied and added on top of stack
 * @param isModelMatrix Switch so that the Xwm matrix is populated correctly
 * @return GZ_FAILURE in case of overflow else GZ_SUCCESS.
 */
int GzPushMatrix(GzRender *render, GzMatrix matrix, bool isModelMatrix)
{
    // Check for overflow condition
    if(render->matlevel >= 99)
    {
        return GZ_FAILURE;
    }

    ++(render->matlevel);

    if(render->matlevel == 0)
    {
        for(unsigned int column = 0; column < 4; ++column)
        {
            for(unsigned int row = 0; row < 4; ++row)
            {
                render->Ximage[render->matlevel][column][row] = matrix[column][row];
            }
        }

        // Push an identity matrix on the normal stack instead of screen projection
        GzMakeIdentityMatrix(render->Xnorm[render->matlevel]);
    }
    else
    {
        Matrix4x4MultiplyBy4x4(render->Ximage[render->matlevel], render->Ximage[render->matlevel - 1], matrix);

        // Push an identity matrix on the normal stack in case of Xpi
        if(render->matlevel == 1)
        {
            GzMakeIdentityMatrix(render->Xnorm[render->matlevel]);
        }
        else
        {
            GzMatrix normalFriendlyMatrix;
            MatrixCopy(normalFriendlyMatrix, matrix);

            // Remove translation components
            for(unsigned int index = 0; index < 3; ++index)
            {
                normalFriendlyMatrix[3][index] = 0;
                normalFriendlyMatrix[index][3] = 0;
            }

            // Compute inverse scalar of matrix
            float a = normalFriendlyMatrix[0][0];
            float b = normalFriendlyMatrix[0][1];
            float c = normalFriendlyMatrix[0][2];
            float scale = sqrt(a*a + b*b + c*c);
            float inverseScale = scale == 0.0f ? 1 : 1 / sqrt(a*a + b*b + c*c);

            MatrixMultiplyByScalar(normalFriendlyMatrix, inverseScale);

            // Assuming that the normalFriendlyMatrix is now normalized
            Matrix4x4MultiplyBy4x4(render->Xnorm[render->matlevel], render->Xnorm[render->matlevel - 1], normalFriendlyMatrix);
        }
    }

    if(isModelMatrix == true)
    {
        ++render->matlevelsXwm;
        if(render->matlevelsXwm == 0)
        {
            for(unsigned int column = 0; column < 4; ++column)
            {
                for(unsigned int row = 0; row < 4; ++row)
                {
                    render->Xwm[render->matlevelsXwm][column][row] = matrix[column][row];
                }
            }
        }
        else
        {
            Matrix4x4MultiplyBy4x4(render->Xwm[render->matlevelsXwm], render->Xwm[render->matlevelsXwm - 1], matrix);
        }
    }

    return GZ_SUCCESS;
}

/**
 * Removes a transformation matrix from the top of the stack.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * Pops the last element in the stack.
 *
 * @param render The render object for which the Xform stack is updated.
 * @return GZ_FAILURE in case of underflow else GZ_SUCCESS.
 */
int GzPopMatrix(GzRender *render)
{
    // Check for underflow condition
    if(render->matlevel == -1)
    {
        return GZ_FAILURE;
    }

    --(render->matlevel);

    return GZ_SUCCESS;
}

/**
 * Set renderer attribute states (e.g. GZ_RGB_COLOR default color).
 * Later set shaders, interpolators, texture maps and lights.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * None.
 *
 * @param render The render object whose attributes are to be set
 * @param numAttributes The number of attributes that have been passed
 * @param nameList The array of all tokens passed
 * @param valueList The list of values that have been passed
 * @return A GZ_SUCCESS value
 */
int GzPutAttribute(GzRender* render, int numAttributes, GzToken* nameList, GzPointer* valueList)
{
    if(render == NULL || render-> open != 1)
    {
        return GZ_FAILURE;
    }

    for(unsigned int index = 0; index < numAttributes; ++index)
    {
        switch(nameList[index])
        {
            case GZ_RGB_COLOR:
            {
                float* color = (float*) valueList[index];
                for(unsigned int colorIndex = 0; colorIndex < 3; ++colorIndex)
                {
                    render->flatcolor[colorIndex] = color[colorIndex];
                }

                break;
            }

            case GZ_AMBIENT_LIGHT:
            {
                GzLight* light = (GzLight*) valueList[index];
                render->ambientlight = *light;
                break;
            }

            case GZ_DIRECTIONAL_LIGHT:
            {
                GzLight *light = (GzLight*) valueList[index];

                if(render->numlights >= MAX_LIGHTS - 1)
                {
                    break;
                }

                ++render->numlights;
                memcpy(&render->lights[render->numlights], light, sizeof(GzLight));

                if(render->lights[render->numlights].isCastShadows == true)
                {
                    const unsigned int lightRes = LIGHT_RES_X * LIGHT_RES_Y;
                    render->lights[render->numlights].shadowZ = new float[lightRes];

                    // Set all z buffer positions to inifinity. Since we don't know how to set it to infinity, use INT_MAX for now.
                    for(unsigned int index = 0; index < lightRes; ++index)
                    {
                        render->lights[render->numlights].shadowZ[index] = INT_MAX;
                    }

                    //TODO: Compute the matrices
                    MakeLightMatrices(&render->lights[render->numlights]);
                }

                break;
            }

            case GZ_AMBIENT_COEFFICIENT:
            {
                GzColor* ambientColor = (GzColor*) valueList[index];
                memcpy(render->Ka, *ambientColor, sizeof(GzColor));
                break;
            }

            case GZ_DIFFUSE_COEFFICIENT:
            {
                GzColor* diffuseColor = (GzColor*) valueList[index];
                memcpy(render->Kd, *diffuseColor, sizeof(GzColor));
                break;
            }

            case GZ_SPECULAR_COEFFICIENT:
            {
                GzColor* specularColor = (GzColor*) valueList[index];
                memcpy(render->Ks, *specularColor, sizeof(GzColor));
                break;
            }

            case GZ_DISTRIBUTION_COEFFICIENT:
            {
                float* specularPower = (float*) valueList[index];
                render->spec = *specularPower;
                break;
            }

            case GZ_INTERPOLATE:
            {
                int* interpolationMethod = (int*) valueList[index];
                render->interp_mode = *interpolationMethod;
                break;
            }

            case GZ_TEXTURE_MAP:
            {
                render->tex_fun = (GzTexture) valueList[index];
                break;
            }

            case GZ_AASHIFTX:
            {
                float *aliasX = (float *) valueList[index];
                render->aliasX = *aliasX;
                break;
            }

            case GZ_AASHIFTY:
            {
                float *aliasY = (float *) valueList[index];
                render->aliasY = *aliasY;
                break;
            }
        }
    }

    return GZ_SUCCESS;
}

/**
 * Pass in a triangle description with tokens and values corresponding to
 *      GZ_NULL_TOKEN:  Do nothing, no values
 *      GZ_POSITION:    3 vertx positions in model space
 *      GZ_NORMAL:      3 normal positions in model space
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * None.
 *
 * @param render The render to which the triangle is rendered
 * @param numParts The number of tokens passed
 * @param nameList The list of tokens passed
 * @param valueList The list values corresponding to the tokens
 * @return A success or failure of the render operation
 */
int GzPutTriangle(GzRender* render, int numParts, GzToken* nameList, GzPointer* valueList)
{
    if(render == NULL || render-> open != 1)
    {
        return GZ_FAILURE;
    }

    GzCoord         vertexList[3]   = {0};
    GzCoord         normalList[3]   = {0};
    GzTextureIndex  uvList[3]       = {0};

    for(unsigned int index = 0; index < numParts; ++index)
    {
        switch(nameList[index])
        {
            case GZ_POSITION:
            {
                // The vertexList is of size 3 i.e. GzCoord array of size 3.
                // Each GzCoord has components x, y and z.
                GzCoord *tempVertexList = (GzCoord*) valueList[index];
                memcpy(vertexList, tempVertexList, sizeof(GzCoord[3]));
                break;
            }

            case GZ_NORMAL:
            {
                GzCoord *tempNormalList = (GzCoord*) valueList[index];
                memcpy(normalList, tempNormalList, sizeof(GzCoord[3]));
            }

            case GZ_TEXTURE_INDEX:
            {
                GzTextureIndex *tempUvList = (GzTextureIndex*) valueList[index];
                memcpy(uvList, tempUvList, sizeof(GzTextureIndex[3]));
            }
        }
    }

    drawTriangle(render, vertexList[0], vertexList[1], vertexList[2], normalList, uvList);

    return GZ_SUCCESS;
}

/**
 * Convert float color to GzIntensity short
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * Not a part of the API.
 *
 * @param Float of the color
 * @return GzIntensity of the color
 */
short ctoi(float color)
{
    return (short) ((int)(color * ((1 << 12) - 1)));
}

/**
 * The function determines the position of a point in an edge and classify the edge as a Top or Left edge.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * The method assumes the coordinate system as Left handed coordinate system but changing the direction
 * assumptions with respect to the reference coordinate system will also work.
 *
 * 1) Slope of the edge is determined and the edge is classified as positive or negative slope
 * 2) The edge points are sorted by their x coordinates
 * 3) The direction of the vertex is computed with respect to the edge using the equation of the
 *      Implicit Form of A Line -> (1)
 *          -> If the slope of the edge is positive, and the equation (1) gives a positive value,
 *              the point lies right of the common edge, i.e. the edge is a left edge
 *          -> If the slope of the edge is negative, and the equation (1) gives a positive value,
 *              the point lies bottom of the common edge , i.e. the edge is a top edge
 *          -> If the slope is infinity, i.e. parallel to the y axis, the edge's left or rightness
 *              compared to the vertex can be easily determined
 *
 * @param firstdgePoint First vertex of the edge
 * @param secondEdgePoint Second vertex of the edge
 * @param nonEdgePoint The vertex whose orientation to the edge is to be determined
 * @return True if the edge is a top or left edge
 */
bool isTopOrLeftEdge(GzCoord firstEdgePoint, GzCoord secondEdgePoint, GzCoord nonEdgePoint)
{
    // If only right or left can be determined
    if(firstEdgePoint[0] - secondEdgePoint[0] == 0)
    {
        // The point is a right point
        if(nonEdgePoint[0] > firstEdgePoint[0])
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    // Sort points by x and perform the operation
    if(firstEdgePoint[0] > secondEdgePoint[0])
    {
        size_t gzCoordSize = sizeof(float) * 3;
        GzCoord temp;

        memcpy(temp, firstEdgePoint, gzCoordSize);
        memcpy(firstEdgePoint, secondEdgePoint, gzCoordSize);
        memcpy(secondEdgePoint, temp, gzCoordSize);
    }

    float dy = secondEdgePoint[1] - firstEdgePoint[1];
    float dx = secondEdgePoint[0] - firstEdgePoint[0];
    float m = dy / dx;
    float c = firstEdgePoint[1] - m * firstEdgePoint[0];

    float d = dy * nonEdgePoint[0]  - dx * nonEdgePoint[1] + dx * c;

    // Refer function documentation for details
    if(d > 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

/**
 * Sorts the vertices.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * Sort based on ascending Y coordinate. If the Y coordinates are equal, sort by ascending X coordinate.
 *
 * @param threePointArray The list of vertices to be sorted
 */
void sortVertices(GzCoord threePointArray[3], GzCoord* threeNormalArray, GzCoord* gzCoordUVList)
{
    size_t gzCoordSize = sizeof(float) * 3;
    GzCoord temp;

    if( (threePointArray[0][1] > threePointArray[1][1]) || ((threePointArray[0][1] == threePointArray[1][1]) && (threePointArray[0][0] > threePointArray[1][0])) )
    {
        memcpy(temp, threePointArray[0], gzCoordSize);
        memcpy(threePointArray[0], threePointArray[1], gzCoordSize);
        memcpy(threePointArray[1], temp, gzCoordSize);

        if(threeNormalArray != NULL)
        {
            memcpy(temp, threeNormalArray[0], gzCoordSize);
            memcpy(threeNormalArray[0], threeNormalArray[1], gzCoordSize);
            memcpy(threeNormalArray[1], temp, gzCoordSize);
        }

        if(gzCoordUVList != NULL)
        {
            memcpy(temp, gzCoordUVList[0], gzCoordSize);
            memcpy(gzCoordUVList[0], gzCoordUVList[1], gzCoordSize);
            memcpy(gzCoordUVList[1], temp, gzCoordSize);
        }
    }

    if( (threePointArray[0][1] > threePointArray[2][1]) || ((threePointArray[0][1] == threePointArray[2][1]) && (threePointArray[0][0] > threePointArray[2][0])) )
    {
        memcpy(temp, threePointArray[0], gzCoordSize);
        memcpy(threePointArray[0], threePointArray[2], gzCoordSize);
        memcpy(threePointArray[2], temp, gzCoordSize);

        if(threeNormalArray != NULL)
        {
            memcpy(temp, threeNormalArray[0], gzCoordSize);
            memcpy(threeNormalArray[0], threeNormalArray[2], gzCoordSize);
            memcpy(threeNormalArray[2], temp, gzCoordSize);
        }

        if(gzCoordUVList != NULL)
        {
            memcpy(temp, gzCoordUVList[0], gzCoordSize);
            memcpy(gzCoordUVList[0], gzCoordUVList[2], gzCoordSize);
            memcpy(gzCoordUVList[2], temp, gzCoordSize);
        }
    }

    if( (threePointArray[1][1] > threePointArray[2][1]) || ((threePointArray[1][1] == threePointArray[2][1]) && (threePointArray[1][0] > threePointArray[2][0])) )
    {
        memcpy(temp, threePointArray[2], gzCoordSize);
        memcpy(threePointArray[2], threePointArray[1], gzCoordSize);
        memcpy(threePointArray[1], temp, gzCoordSize);

        if(threeNormalArray != NULL)
        {
            memcpy(temp, threeNormalArray[2], gzCoordSize);
            memcpy(threeNormalArray[2], threeNormalArray[1], gzCoordSize);
            memcpy(threeNormalArray[1], temp, gzCoordSize);
        }

        if(gzCoordUVList != NULL)
        {
            memcpy(temp, gzCoordUVList[2], gzCoordSize);
            memcpy(gzCoordUVList[2], gzCoordUVList[1], gzCoordSize);
            memcpy(gzCoordUVList[1], temp, gzCoordSize);
        }
    }
}

/**
 * Determines the orientation of the point to the edge.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * None.
 *
 * @param point Point whose orientation to the edge is determined
 * @param headLinePoint First point of the edge
 * @param secondLinePoint Second point of the edge
 * @return The value computed using the 
 */
int pointOrientationToLine(GzCoord point, GzCoord headLinePoint, GzCoord tailLinePoint)
{
    float dX = headLinePoint[0] - tailLinePoint[0];
    float dY = headLinePoint[1] - tailLinePoint[1];

    float A = dY;
    float B = -dX;
    float C = dX * tailLinePoint[1] - dY * tailLinePoint[0];

    float orientation = dY * point[0] - dX * point[1] + C;

    int returnValue = 0;

    if(orientation > 0)
    {
        returnValue = 1;
    }
    else if(orientation < 0)
    {
        returnValue = -1;
    }
    else
    {
        returnValue = 0;
    }

    return returnValue;
}

/**
 * Checks if a pixel is to be shaded.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * Sorted vertices of the triangle are passed in the threePointArray.
 *
 * If a point to be shaded lies on a whole number x and y, determine whether it needs to be displayed
 * using the pointOrientationToLine.
 *
 * Consider A to be the point which we need to check if to be displayed. Consider PQR as the vertices
 * of the triangle. Say, A lies on the edge PQ the triangle. We pass the parameters to
 * pointOrientationToLine() as A, P, Q.
 *
 * @param testPoint The point that is tested if shadable
 * @param threePointArray The sorted coordinates of the triangle to be shaded
 * @return True if the pixel is shaded else false.
 */
bool isPixelShadable(GzCoord testPoint, GzCoord threePointArray[3])
{
    bool isShade = false;

    int firstOrientation = pointOrientationToLine(testPoint, threePointArray[0], threePointArray[1]);
    if(firstOrientation == 0)
    {
        if(isTopOrLeftEdge(threePointArray[0], threePointArray[1], threePointArray[2]) == true)
        {
            return false;
        }
    }

    int secondOrientation = pointOrientationToLine(testPoint, threePointArray[1], threePointArray[2]);
    if(secondOrientation == 0)
    {
        if(isTopOrLeftEdge(threePointArray[1], threePointArray[2], threePointArray[0]) == true)
        {
            return false;
        }
    }

    int thirdOrientation = pointOrientationToLine(testPoint, threePointArray[2], threePointArray[0]);
    if(thirdOrientation == 0)
    {
        if(isTopOrLeftEdge(threePointArray[2], threePointArray[0], threePointArray[1]) == true)
        {
            return false;
        }
    }

    // Determine if the orientation is zero and handle that case here
    if(firstOrientation >= 0 && secondOrientation >= 0 && thirdOrientation >= 0)
    {
        isShade = true;
    }
    else if(firstOrientation <= 0 && secondOrientation <= 0 && thirdOrientation <= 0)
    {
        isShade = true;
    }

    return isShade;
}

/**
 * Computes the color at a point by adding the specular, diffuse and ambient components
 * of all light sources with the texture interaction at that point.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * The light and eye are assumed to be very far from the point.
 *
 * @param outputColor Output Parameter. Color computed at a point
 * @param render The render object containing the material properties and lighting information
 * @param normal The normal at the point of interest
 * @param u The texture parameter
 * @param v The texture parameter
 */
void ComputePhongTextureShading(GzCoord outputColor, GzRender* render, GzCoord normal, float u, float v,GzCoord orivert)
{
    GzCoord Kd = {0};

    render->tex_fun(u, v, Kd);

    GzCoord transformedNormal = {0};
    Matrix4x4MultiplyBy4x1(transformedNormal, render->Xnorm[render->matlevel], normal);
    NormalizeVector(transformedNormal);

    GzCoord eyeVector = {0, 0, -1};
    NormalizeVector(eyeVector);

    // Compute specular and diffuse material shading based on the lighting
    for(unsigned int lightIndex=0; lightIndex <= render->numlights; ++lightIndex)
    {
        bool isLightVisible = isAffectedByLight(render,&render->lights[lightIndex], orivert);

        if(isLightVisible == false)
        {
            continue;
        }

        NormalizeVector(render->lights[lightIndex].direction);

        GzCoord lightVector = {0};
        CopyVector(lightVector, render->lights[lightIndex].direction);
        NormalizeVector(lightVector);

        float normalDotLight    = VectorDotProduct(transformedNormal, lightVector);
        float normalDotEye      = VectorDotProduct(transformedNormal, eyeVector);

        GzCoord intermediateNormalVector = {0};
        CopyVector(intermediateNormalVector, transformedNormal);
        MultiplyScalarToVector(normalDotLight * 2, intermediateNormalVector);

        GzCoord reflectionVector = {0};
        CreateVector(reflectionVector, intermediateNormalVector, render->lights[lightIndex].direction);

        float reflectionDotEye = VectorDotProduct(reflectionVector, eyeVector);

        // Compute cases where the eye and the light do not see each other
        if(normalDotLight * normalDotEye < 0)
        {
            continue;
        }
        else if(normalDotLight < 0)
        {
            MultiplyScalarToVector(-1.0f, transformedNormal);

            normalDotLight  = -normalDotLight;
            normalDotEye    = -normalDotEye;
        }

        for(unsigned int colorIndex=0; colorIndex < 3; ++colorIndex)
        {
            // Diffuse light at point
            float totalDiffuseLight = render->lights[lightIndex].color[colorIndex] * normalDotLight;
            totalDiffuseLight = totalDiffuseLight < 0 ? 0 : totalDiffuseLight;

            // Specular light at point
            float totalSpecularLight = render->lights[lightIndex].color[colorIndex] * pow(reflectionDotEye, render->spec);
            totalSpecularLight = totalSpecularLight < 0 ? 0 : totalSpecularLight;

            // Net light at point
            outputColor[colorIndex] +=
                Kd[colorIndex] * totalDiffuseLight +
                render->Ks[colorIndex] * totalSpecularLight;
        }
    }

    // Add ambient material shading based on diffuse material properties
    for(unsigned int colorIndex=0; colorIndex < 3; ++colorIndex)
    {
        outputColor[colorIndex] += Kd[colorIndex] * render->ambientlight.color[colorIndex];

        // Clamp the color to 1.0f if it exceeds 1.0f
        outputColor[colorIndex] > 1.0f ? 1.0f : outputColor[colorIndex];
    }
}

/**
 * Computes the color at a point by adding the specular, diffuse and ambient components
 * of all light sources at that point.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * The light and eye are assumed to be very far from the point.
 *
 * @param outputColor Output Parameter. Color computed at a point
 * @param render The render object containing the material properties and lighting information
 * @param normal The normal at the point of interest
 */
void ComputeShading(GzCoord outputColor, GzRender* render, GzCoord normal,GzCoord currentpixelcoord, bool applyMaterialProperties = true)
{
    GzCoord transformedNormal = {0};
    Matrix4x4MultiplyBy4x1(transformedNormal, render->Xnorm[render->matlevel], normal);
    NormalizeVector(transformedNormal);

    GzCoord eyeVector = {0, 0, -1};
    NormalizeVector(eyeVector);

    // Compute specular and diffuse material shading based on the lighting
    for(unsigned int lightIndex=0; lightIndex <= render->numlights; ++lightIndex)
    {
        bool isLightVisible = isAffectedByLight(render,&render->lights[lightIndex],currentpixelcoord);

        if(isLightVisible == false)
        {
            continue;
        }

        NormalizeVector(render->lights[lightIndex].direction);

        GzCoord lightVector = {0};
        CopyVector(lightVector, render->lights[lightIndex].direction);
        NormalizeVector(lightVector);

        float normalDotLight    = VectorDotProduct(transformedNormal, lightVector);
        float normalDotEye      = VectorDotProduct(transformedNormal, eyeVector);

        GzCoord intermediateNormalVector = {0};
        CopyVector(intermediateNormalVector, transformedNormal);
        MultiplyScalarToVector(normalDotLight * 2, intermediateNormalVector);

        GzCoord reflectionVector = {0};
        CreateVector(reflectionVector, intermediateNormalVector, render->lights[lightIndex].direction);

        float reflectionDotEye = VectorDotProduct(reflectionVector, eyeVector);

        // Compute cases where the eye and the light do not see each other
        if(normalDotLight * normalDotEye < 0)
        {
            continue;
        }
        else if(normalDotLight < 0)
        {
            MultiplyScalarToVector(-1.0f, transformedNormal);

            normalDotLight  = -normalDotLight;
            normalDotEye    = -normalDotEye;
        }

        for(unsigned int colorIndex=0; colorIndex < 3; ++colorIndex)
        {
            // Diffuse light at point
            float totalDiffuseLight = render->lights[lightIndex].color[colorIndex] * normalDotLight;
            totalDiffuseLight = totalDiffuseLight < 0 ? 0 : totalDiffuseLight;

            // Specular light at point
            float totalSpecularLight = render->lights[lightIndex].color[colorIndex] * pow(reflectionDotEye, render->spec);
            totalSpecularLight = totalSpecularLight < 0 ? 0 : totalSpecularLight;

            if(applyMaterialProperties == true)
            {
                // Net light at point
                outputColor[colorIndex] +=
                    render->Kd[colorIndex] * totalDiffuseLight +
                    render->Ks[colorIndex] * totalSpecularLight;
            }
            else
            {
                outputColor[colorIndex] += totalDiffuseLight + totalSpecularLight;
            }
        }
    }

    // Add ambient material shading based on ambient lighting
    for(unsigned int colorIndex=0; colorIndex < 3; ++colorIndex)
    {
        if(applyMaterialProperties == true)
        {
            outputColor[colorIndex] += render->Ka[colorIndex] * render->ambientlight.color[colorIndex];
        }
        else
        {
            outputColor[colorIndex] += render->ambientlight.color[colorIndex];
        }

        // Clamp the color to 1.0f if it exceeds 1.0f
        outputColor[colorIndex] > 1.0f ? 1.0f : outputColor[colorIndex];
    }
}

/**
 * Given three parameters at the vertices of a triangle, the method computes an interpolated paramter
 * inside the triangle at a given (x,y) location.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * The xIndex and yIndex are within the perimeter of the triangle.
 *
 * @param interpolationOutput Output Parameter. Interpolated output of the parameters at the vertices of the triangle
 * @param triangleVertices The vertices of the triangle within which interpolation is required
 * @param interpolationInput The interpolation parameters at each of the vertices of the triangle
 * @param xIndex X index of the point where the interpolation is computed
 * @param yIndex Y index of the point where the interpolation is computed
 */
void interpolateTriangle(float interpolationOutput[3], GzCoord triangleVertices[3], float interpolationInput[3][3], float xIndex, float yIndex)
{
    float dy01 = triangleVertices[1][1] - triangleVertices[0][1];
    float dx01 = triangleVertices[1][0] - triangleVertices[0][0];

    float m01 = 0;
    float c01 = 0;
    bool  is01ParallelToYAxis = false;

    if(dx01 == 0)
    {
        is01ParallelToYAxis = true;
    }
    else
    {
        m01 = dy01 / dx01;
        c01 = triangleVertices[0][1] - (m01 * triangleVertices[0][0]);
    }

    float dy2p = yIndex - triangleVertices[2][1];
    float dx2p = xIndex - triangleVertices[2][0];

    float x01 = 0.0f, y01 = 0.0f;

    float mdp = 0;
    float cdp = 0;
    bool  isdpParallelToYAxis = false;

    if(dx2p == 0)
    {
        isdpParallelToYAxis = true;
    }
    else
    {
        mdp = (dy2p / dx2p);
        cdp = yIndex - (mdp * xIndex);
    }

    // Assuming non parallel lines, i.e. m01 != mdp
    x01 = (c01 - cdp) / (mdp - m01);
    y01 = mdp * x01 + cdp;

    // Compute distances and ratios
    // Distance formula
    float distance01 = sqrt(
        (triangleVertices[1][1] - triangleVertices[0][1]) * (triangleVertices[1][1] - triangleVertices[0][1]) +
        (triangleVertices[1][0] - triangleVertices[0][0]) * (triangleVertices[1][0] - triangleVertices[0][0]) );

    // Distance formula
    float distance0to01 = sqrt(
        (y01 - triangleVertices[0][1]) * (y01 - triangleVertices[0][1]) +
        (x01 - triangleVertices[0][0]) * (x01 - triangleVertices[0][0]) );

    float ratio0to01 = distance0to01 / distance01;

    // Distance formula
    float distance01to2 = sqrt(
        (y01 - triangleVertices[2][1]) * (y01 - triangleVertices[2][1]) +
        (x01 - triangleVertices[2][0]) * (x01 - triangleVertices[2][0]) );

    // Distance formula
    float distance2toP = sqrt(
        (yIndex - triangleVertices[2][1]) * (yIndex - triangleVertices[2][1]) +
        (xIndex - triangleVertices[2][0]) * (xIndex - triangleVertices[2][0]) );

    float pointRatio = distance2toP / distance01to2;

    // Get interpolation of point on 0,1 line. Extend it to p
    float intepolation01[3] = {0};

    for(unsigned int index = 0; index < 3; ++index)
    {
        intepolation01[index]  = interpolationInput[1][index] * ratio0to01 + interpolationInput[0][index] * (1 - ratio0to01);
        interpolationOutput[index]   = intepolation01[index] * pointRatio + interpolationInput[2][index] * (1 - pointRatio);
    }
}

/**
 * Draws a triangle at the given coordinates.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * Compute the bounding rectangle of the triangle. For each of the points at a (X, Y) in the
 * triangle where X, Y belongs to the set of integers, determine if the pixel is to be shaded.
 * Also make the decision based on the Z value at that point.
 *
 * @param render The render to which the triangle is drawn
 * @param firstVertex The first vertex of the the triangle that is drawn
 * @param secondVertex The second vertex of the the triangle that is drawn
 * @param thirdVertex The third vertex of the the triangle that is drawn
 */
void drawTriangle(GzRender* render, GzCoord firstVertex, GzCoord secondVertex, GzCoord thirdVertex, GzCoord threeNormalArray[3], GzTextureIndex uvList[3])
{
    GzCoord originalFirstVertex, originalSecondVertex, originalThirdVertex;
    GzCoord gzCoordUvList[3] = {0};
    for(unsigned int index = 0; index < 3; ++index)
    {
        originalFirstVertex[index]  = firstVertex[index];
        originalSecondVertex[index] = secondVertex[index];
        originalThirdVertex[index]  = thirdVertex[index];

        // Convert the UV list to a GzCoord with z = 0 since interpolation would be easy
        for(unsigned int uvIndex = 0; uvIndex < 2; ++uvIndex)
        {
            gzCoordUvList[index][uvIndex] = uvList[index][uvIndex];
        }
    }

    OriginalVertexAtWorldSpace(render, originalFirstVertex, originalSecondVertex, originalThirdVertex);

    Matrix4x4MultiplyBy4x1(firstVertex,  render->Ximage[render->matlevel], originalFirstVertex);
    Matrix4x4MultiplyBy4x1(secondVertex, render->Ximage[render->matlevel], originalSecondVertex);
    Matrix4x4MultiplyBy4x1(thirdVertex,  render->Ximage[render->matlevel], originalThirdVertex);

    GzCoord Vertex1AtW,Vertex2AtW,Vertex3AtW ;

    VertexBackAtWorldSpace(render, firstVertex, &Vertex1AtW);
    VertexBackAtWorldSpace(render, secondVertex, &Vertex2AtW);
    VertexBackAtWorldSpace(render, thirdVertex, &Vertex3AtW);

    int boundRectangleLeftX     = (int) floor( min( min(firstVertex[0], secondVertex[0]), thirdVertex[0] ) );
    int boundRectangleRightX    = (int) ceil( max( max(firstVertex[0], secondVertex[0]), thirdVertex[0] ) );
    int boundRectangleTopY      = (int) floor( min( min(firstVertex[1], secondVertex[1]), thirdVertex[1] ) );
    int boundRectangleBottomY   = (int) ceil( max( max(firstVertex[1], secondVertex[1]), thirdVertex[1] ) );

    // If the z value of any vertex is less than 0, we need not display the triangle
    if (firstVertex[2] < 0 || secondVertex[2] < 0 || thirdVertex[2] < 0)
    {
        return;
    }

    GzCoord threePointArray[3];
    for(unsigned int coordIndex = 0; coordIndex <= 2; ++coordIndex)
    {
        threePointArray[0][coordIndex] = firstVertex[coordIndex];
        threePointArray[1][coordIndex] = secondVertex[coordIndex];
        threePointArray[2][coordIndex] = thirdVertex[coordIndex];
    }

    float planeD = 0;
    GzCoord normal;
    computeNormal(normal, planeD, threePointArray);
    sortVertices(threePointArray, threeNormalArray, gzCoordUvList);

    // Color shading is disabled for Shadow Mapping
    // This is required only for Color Shading
    GzColor threeColorArray[3] = {0};
    //for(unsigned int vertexIndex=0; vertexIndex < 3; ++vertexIndex)
    //{
    //    if(render->tex_fun == NULL)
    //    {
    //        ComputeShading(threeColorArray[vertexIndex], render, threeNormalArray[vertexIndex]);
    //    }
    //    else
    //    {
    //        ComputeShading(threeColorArray[vertexIndex], render, threeNormalArray[vertexIndex], false);
    //    }
    //}

    // Transform UV of vertices from image to perspective space
    const float zMax = INT_MAX;
    float vzPrime = 0;

    GzCoord screenUVArray[3] = {0};
    for(unsigned int index = 0; index < 3; ++index)
    {
        vzPrime = threePointArray[index][2] / (zMax - threePointArray[index][2]);
        GzConvertImageToScreen(screenUVArray[index], gzCoordUvList[index], vzPrime);
    }

    // Check if each point in the bounding rectangle can be displayed and display it
    for(int yIndex = boundRectangleTopY; yIndex <= boundRectangleBottomY; ++yIndex)
    {
        for(int xIndex = boundRectangleLeftX; xIndex <= boundRectangleRightX; ++xIndex)
        {
            // Compute the Z value using the plane equation
            GzDepth z = getZForVertex(normal, planeD, xIndex, yIndex);

            GzCoord currentPixelIndex;
            currentPixelIndex[0] = (float) xIndex + render->aliasX;
            currentPixelIndex[1] = (float) yIndex + render ->aliasY;
            currentPixelIndex[2] = (float) z;

            float currentPixelX = currentPixelIndex[0];
            float currentPixelY = currentPixelIndex[1];

            bool isShade = isPixelShadable(currentPixelIndex, threePointArray);

            if(isShade == true)
            {
                GzColor colorp = {0};

                if(render->tex_fun == NULL)
                {
                    if(render->interp_mode == GZ_COLOR)
                    {
                        interpolateTriangle(colorp, threePointArray, threeColorArray, currentPixelX, currentPixelY);
                    }
                    else if(render->interp_mode == GZ_NORMAL)
                    {
                        GzCoord interpolatedNormal = {0};
                        interpolateTriangle(interpolatedNormal, threePointArray, threeNormalArray, currentPixelX, currentPixelY);
                        ComputeShading(colorp, render, interpolatedNormal,currentPixelIndex);
                    }
                    else
                    {
                        CopyVector(colorp, render->flatcolor);
                    }
                }
                else
                {
                    // Interpolate the UV of perspective space
                    GzCoord interpolatedPerspectiveUV = {0};
                    interpolateTriangle(interpolatedPerspectiveUV, threePointArray, screenUVArray, currentPixelX, currentPixelY);

                    // Convert UV back to image space
                    vzPrime = z / (zMax - z);
                    GzCoord imageUV = {0};
                    GzConvertScreenToImage(imageUV, interpolatedPerspectiveUV, vzPrime);

                    if(render->interp_mode == GZ_COLOR)
                    {
                        // Interpolate the color
                        interpolateTriangle(colorp, threePointArray, threeColorArray, currentPixelX, currentPixelY);

                        GzCoord Kt = {0};
                        render->tex_fun(imageUV[0], imageUV[1], Kt);

                        for(unsigned int index = 0; index < 3; ++index)
                        {
                            colorp[index] = colorp[index] * Kt[index];
                        }
                    }
                    else if(render->interp_mode == GZ_NORMAL)
                    {
                        GzCoord interpolatedNormal = {0};
                        interpolateTriangle(interpolatedNormal, threePointArray, threeNormalArray, currentPixelX, currentPixelY);
                        ComputePhongTextureShading(colorp, render, interpolatedNormal, imageUV[0], imageUV[1],currentPixelIndex);
                    }
                    else
                    {
                        // Handle texture flat shading here
                        memcpy(colorp, 0x00, sizeof(colorp));
                    }
                }

                GzDepth existingZ = 0;
                GzIntensity r = 0, g = 0, b = 0, a = 0;

                GzGetDisplay(render->display, xIndex, yIndex, &r, &g, &b, &a, &existingZ);

                //If the current Z is lesser than existing Z, it should be displayed.
                if(z <= existingZ)
                {
                    r = (GzIntensity) (colorp[0] * 4095);
                    g = (GzIntensity) (colorp[1] * 4095);
                    b = (GzIntensity) (colorp[2] * 4095);
                    GzPutDisplay(render->display, xIndex, yIndex, r, g, b, 0xFFFF, z);
                }
            }
        }
    }
}

//TODO:
int GzUpdateShadowBuffers(GzRender* render, int numParts, GzToken* nameList, GzPointer* valueList)
{
    if(render == NULL || render-> open != 1)
    {
        return GZ_FAILURE;
    }

    GzCoord         vertexList[3]   = {0};
    GzCoord         normalList[3]   = {0};
    GzTextureIndex  uvList[3]       = {0};

    for(unsigned int index = 0; index < (unsigned int)numParts; ++index)
    {
        switch(nameList[index])
        {
            case GZ_POSITION:
            {
                // The vertexList is of size 3 i.e. GzCoord array of size 3.
                // Each GzCoord has components x, y and z.
                GzCoord *tempVertexList = (GzCoord*) valueList[index];
                memcpy(vertexList, tempVertexList, sizeof(GzCoord[3]));
                break;
            }
        }
    }

    updateLightBuffers(render, vertexList);

    return GZ_SUCCESS;
}

//TODO:
void updateLightBuffers(GzRender *render, GzCoord vertexList[3])
{
    for(unsigned int lightIndex = 0; lightIndex < 1; ++lightIndex)
    {
        GzLight *light = &(render->lights[lightIndex]);

        if(light->isCastShadows == false)
        {
            continue;
        }

        GzMatrix Xlbm, temp, temp2;
        Matrix4x4MultiplyBy4x4(temp, render->Xlbp, light->Xpl);
        Matrix4x4MultiplyBy4x4(temp2, temp, light->Xlw);
        Matrix4x4MultiplyBy4x4(Xlbm, temp2, render->Xwm[render->matlevelsXwm]);

        // Code modified from renderer
        float *firstVertex  = vertexList[0];
        float *secondVertex = vertexList[1];
        float *thirdVertex  = vertexList[2];

        GzCoord originalFirstVertex, originalSecondVertex, originalThirdVertex;
        for(unsigned int index = 0; index < 3; ++index)
        {
            originalFirstVertex[index]  = firstVertex[index];
            originalSecondVertex[index] = secondVertex[index];
            originalThirdVertex[index]  = thirdVertex[index];
        }

        Matrix4x4MultiplyBy4x1(firstVertex,  Xlbm, originalFirstVertex);
        Matrix4x4MultiplyBy4x1(secondVertex, Xlbm, originalSecondVertex);
        Matrix4x4MultiplyBy4x1(thirdVertex,  Xlbm, originalThirdVertex);

        int boundRectangleLeftX     = (int) floor( min( min(firstVertex[0], secondVertex[0]), thirdVertex[0] ) );
        int boundRectangleRightX    = (int) ceil( max( max(firstVertex[0], secondVertex[0]), thirdVertex[0] ) );
        int boundRectangleTopY      = (int) floor( min( min(firstVertex[1], secondVertex[1]), thirdVertex[1] ) );
        int boundRectangleBottomY   = (int) ceil( max( max(firstVertex[1], secondVertex[1]), thirdVertex[1] ) );

        // If the z value of any vertex is less than 0, we need not consider the triangle since one point is behind the light source
        if (firstVertex[2] < 0 || secondVertex[2] < 0 || thirdVertex[2] < 0)
        {
            continue;
        }

        GzCoord threePointArray[3];
        for(unsigned int coordIndex = 0; coordIndex <= 2; ++coordIndex)
        {
            threePointArray[0][coordIndex] = firstVertex[coordIndex];
            threePointArray[1][coordIndex] = secondVertex[coordIndex];
            threePointArray[2][coordIndex] = thirdVertex[coordIndex];
        }

        float planeD = 0;
        GzCoord normal;
        computeNormal(normal, planeD, threePointArray);
        sortVertices(threePointArray, NULL, NULL);

        // Check if each point in the bounding rectangle can be displayed and display it
        for(int yIndex = boundRectangleTopY; yIndex <= boundRectangleBottomY; ++yIndex)
        {
            for(int xIndex = boundRectangleLeftX; xIndex <= boundRectangleRightX; ++xIndex)
            {
                // Compute the Z value using the plane equation
                float z = getFloatZForVertex(normal, planeD, xIndex, yIndex);

                GzCoord currentPixelIndex;
                currentPixelIndex[0] = (float) xIndex;
                currentPixelIndex[1] = (float) yIndex;

                float currentPixelX = currentPixelIndex[0];
                float currentPixelY = currentPixelIndex[1];

                bool isStorable = isPixelShadable(currentPixelIndex, threePointArray);
                if(xIndex<LIGHT_RES_X && xIndex>=0 && yIndex<LIGHT_RES_Y && yIndex>=0)
                {
                    if(isStorable == true)
                    {
                        float existingZ = light->shadowZ[LIGHT_RES_X * yIndex + xIndex];

                       //If the current Z is lesser than existing Z, it should be stored.
                        if(z <= existingZ)
                        {
                            light->shadowZ[LIGHT_RES_X * yIndex + xIndex] = z;
                        }
                    }
                }
            }
        }
    }
}

//TODO:
void OriginalVertexAtWorldSpace(GzRender * render, GzCoord firstVertex,GzCoord secondVertex, GzCoord thirdVertex)
{
    GzCoord new1W,new2W,new3W;

    Matrix4x4MultiplyBy4x1(new1W,render->Xwm[render->matlevelsXwm],firstVertex);
    Matrix4x4MultiplyBy4x1(new2W,render->Xwm[render->matlevelsXwm],secondVertex);
    Matrix4x4MultiplyBy4x1(new3W,render->Xwm[render->matlevelsXwm],thirdVertex);
}

//TODO:
void VertexAtPerpectiveSpace(GzRender * render, GzCoord ScreenSpaceVertex,GzCoord PerpVertex)
{
    float d = 1 / tan( (render->camera.FOV / 2) * PI / 180.0);

    float xres = (render->display->xres+1)/2;
    float yres = (render->display->yres+1)/2;

    PerpVertex[0] = (ScreenSpaceVertex[0]-xres)/xres;
    PerpVertex[1] = (ScreenSpaceVertex[1]-yres)/-yres;
    PerpVertex[2] = ScreenSpaceVertex[2]*d/MAXINT;
}

//TODO:
void inversePerspective(GzCoord initialVertex, GzCoord transformedVertex,float zval)
{
    float zValue, z1;
    zValue = zval;

    z1 = zValue/(MAXINT - zValue);

    for (int i =0 ;i< 3; i++)
    {
        transformedVertex[i] =initialVertex[i]*(z1+1);
    }
}

//TODO:
void VertexBackAtWorldSpace(GzRender * render, GzCoord ScreenSpaceVertex,GzCoord *outputVertex)
{
    GzCoord PerspVertex;
    GzCoord CameraVertex;
    GzCoord WorldVertex;

    //Screen to Perspective Transformation
    VertexAtPerpectiveSpace(render, ScreenSpaceVertex,PerspVertex);

    //Perspective to Image Transformation
    inversePerspective(PerspVertex,CameraVertex,ScreenSpaceVertex[2]);

    //Image to World Transformation
    Matrix4x4MultiplyBy4x1(WorldVertex,render->camera.Xwi,CameraVertex);

    memcpy(outputVertex,WorldVertex,sizeof(GzCoord));
}

//TODO:
/*a point from screen space to world space, then from world space to Light space to LB space.
    Then calculate the z buffer of the point. 
    If it is equal to the z buffer of the material, then true else false */
bool isAffectedByLight(GzRender* render,GzLight* light,GzCoord orivert)
{
    if(light->isCastShadows==false)
    {
        return true;
    }

    GzCoord newLvert,newLBvert,newVWvert, newVertInP, newVertInLb;

    VertexBackAtWorldSpace(render, orivert, &newVWvert);
    Matrix4x4MultiplyBy4x1(newLvert, light->Xlw, newVWvert);
    Matrix4x4MultiplyBy4x1(newVertInP, light->Xpl, newLvert);
    Matrix4x4MultiplyBy4x1(newVertInLb, render->Xlbp, newVertInP);

    unsigned int xIndex = (floor) (newVertInLb[0] + 0.5);
    unsigned int yIndex = (floor) (newVertInLb[1] + 0.5);

    if(xIndex < 0 || xIndex >= LIGHT_RES_X || yIndex < 0 || yIndex >= LIGHT_RES_Y)
    {
        return false;
    }

    float existingZ = light->shadowZ[yIndex * LIGHT_RES_X + xIndex];

    float difference = abs(existingZ - newVertInLb[2]);

    if(difference <= 100000.0f)
    {
        return true;
    }
    else
    {
        return false;
    }
}
