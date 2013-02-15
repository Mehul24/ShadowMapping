#include "StdAfx.h"
#include <math.h>
#include "gzMath.h"

#ifndef PI
#define PI 3.14159265
#endif

/**
 * Copies a vector.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * None.
 *
 * @param dest The destination vector
 * @param source The source vector
 */
void CopyVector(GzCoord dest, GzCoord source)
{
    for(unsigned int index = 0; index < 3; ++index)
    {
        dest[index] = source[index];
    }
}

/**
 * Returns the magnitude of a vector.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * None.
 *
 * @param vector The vector whose magnitude is computed
 * @return Magnitude of the vector
 */
float VectorMagnitude(GzCoord vector)
{
    return sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
}

/**
 * Normalizes a vector such that the magnitude of the vector will be 1 and direction will remain the same.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * None.
 *
 * @param vector Input Output Parameter. The vector that is normalized
 */
void NormalizeVector(GzCoord vector)
{
    float magnitude = VectorMagnitude(vector);

    for(unsigned int index = 0; index < 3; ++index)
    {
        vector[index] = vector[index] / magnitude;
    }
}

/**
 * Computes the cross product of the vector.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * None.
 *
 * @param result The resulting cross product vector
 * @param The first vector in the cross product computation
 * @param The first vector in the cross product computation
 */
void VectorCrossProduct(GzCoord result, GzCoord vectorOne, GzCoord vectorTwo)
{
    result[0] = vectorOne[1] * vectorTwo[2] - vectorOne[2] * vectorTwo[1];
    result[1] = vectorOne[2] * vectorTwo[0] - vectorOne[0] * vectorTwo[2];
    result[2] = vectorOne[0] * vectorTwo[1] - vectorOne[1] * vectorTwo[0];
}

/**
 * Computes the dot product of the vector.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * None.
 *
 * @param vectorOne The first vector in the dot product computation
 * @param vectorTwo The second vector in the dot product computation
 * @return The dot product of the vector
 */
float VectorDotProduct(GzCoord vectorOne, GzCoord vectorTwo)
{
    float dotProduct = 0.0f;

    for(unsigned int index = 0; index < 3; ++index)
    {
        dotProduct += vectorOne[index] * vectorTwo[index];
    }

    return dotProduct;
}

/**
 * Creates a vector from two given points.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * None.
 *
 * @param vector The vector computed
 * @param headPoint The head point of the vector
 * @param tailPoint The tail point of the vector
 */
void CreateVector(GzCoord vector, GzCoord headPoint, GzCoord tailPoint)
{
    for(unsigned int index = 0; index < 3; ++index)
    {
        vector[index] = headPoint[index] - tailPoint[index];
    }
}

/**
 * Multiplies a scalar value to a given vector.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * None.
 *
 * @param scalar The scalar value that is multiplied
 * @param vector The vector to which the scalar value is multiplied
 */
void MultiplyScalarToVector(float scalar, GzCoord vector)
{
    for(unsigned int index = 0; index < 3; ++index)
    {
        vector[index] *= scalar;
    }
}

/**
 * Creates a 4x4 identity matrix.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * None.
 *
 * @param matrix Contains the identity matrix after the computation.
 */
void GzMakeIdentityMatrix(GzMatrix matrix)
{
    for(unsigned int column = 0; column < 4; ++column)
    {
        for(unsigned int row = 0; row < 4; ++row)
        {
            if(column == row)
            {
                matrix[column][row] = 1;
            }
            else
            {
                matrix[column][row] = 0;
            }
        }
    }
}

/**
 * Copies a 4x4 matrix
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * None.
 *
 * @param dest Destination matrix of the copy operation.
 * @param source Source matrix of the copy operation.
 */
void MatrixCopy(GzMatrix dest, GzMatrix source)
{
    for(unsigned int column = 0; column < 4; ++column)
    {
        for(unsigned int row = 0; row < 4; ++row)
        {
            dest[column][row] = source[column][row];
        }
    }
}

/**
 * Multiplies two 4x4 matrices.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * None.
 *
 * @param dest Destination matrix of the multiplication operation.
 * @param sourceOne First matrix of the multiplication.
 * @param sourceTwo Second matrix of the multiplication.
 */
void Matrix4x4MultiplyBy4x4(GzMatrix dest, GzMatrix sourceOne, GzMatrix sourceTwo)
{
    dest[0][0] = sourceOne[0][0] * sourceTwo[0][0] + sourceOne[0][1] * sourceTwo[1][0] + sourceOne[0][2] * sourceTwo[2][0] + sourceOne[0][3] * sourceTwo[3][0];
    dest[0][1] = sourceOne[0][0] * sourceTwo[0][1] + sourceOne[0][1] * sourceTwo[1][1] + sourceOne[0][2] * sourceTwo[2][1] + sourceOne[0][3] * sourceTwo[3][1];
    dest[0][2] = sourceOne[0][0] * sourceTwo[0][2] + sourceOne[0][1] * sourceTwo[1][2] + sourceOne[0][2] * sourceTwo[2][2] + sourceOne[0][3] * sourceTwo[3][2];
    dest[0][3] = sourceOne[0][0] * sourceTwo[0][3] + sourceOne[0][1] * sourceTwo[1][3] + sourceOne[0][2] * sourceTwo[2][3] + sourceOne[0][3] * sourceTwo[3][3];
    dest[1][0] = sourceOne[1][0] * sourceTwo[0][0] + sourceOne[1][1] * sourceTwo[1][0] + sourceOne[1][2] * sourceTwo[2][0] + sourceOne[1][3] * sourceTwo[3][0];
    dest[1][1] = sourceOne[1][0] * sourceTwo[0][1] + sourceOne[1][1] * sourceTwo[1][1] + sourceOne[1][2] * sourceTwo[2][1] + sourceOne[1][3] * sourceTwo[3][1];
    dest[1][2] = sourceOne[1][0] * sourceTwo[0][2] + sourceOne[1][1] * sourceTwo[1][2] + sourceOne[1][2] * sourceTwo[2][2] + sourceOne[1][3] * sourceTwo[3][2];
    dest[1][3] = sourceOne[1][0] * sourceTwo[0][3] + sourceOne[1][1] * sourceTwo[1][3] + sourceOne[1][2] * sourceTwo[2][3] + sourceOne[1][3] * sourceTwo[3][3];
    dest[2][0] = sourceOne[2][0] * sourceTwo[0][0] + sourceOne[2][1] * sourceTwo[1][0] + sourceOne[2][2] * sourceTwo[2][0] + sourceOne[2][3] * sourceTwo[3][0];
    dest[2][1] = sourceOne[2][0] * sourceTwo[0][1] + sourceOne[2][1] * sourceTwo[1][1] + sourceOne[2][2] * sourceTwo[2][1] + sourceOne[2][3] * sourceTwo[3][1];
    dest[2][2] = sourceOne[2][0] * sourceTwo[0][2] + sourceOne[2][1] * sourceTwo[1][2] + sourceOne[2][2] * sourceTwo[2][2] + sourceOne[2][3] * sourceTwo[3][2];
    dest[2][3] = sourceOne[2][0] * sourceTwo[0][3] + sourceOne[2][1] * sourceTwo[1][3] + sourceOne[2][2] * sourceTwo[2][3] + sourceOne[2][3] * sourceTwo[3][3];
    dest[3][0] = sourceOne[3][0] * sourceTwo[0][0] + sourceOne[3][1] * sourceTwo[1][0] + sourceOne[3][2] * sourceTwo[2][0] + sourceOne[3][3] * sourceTwo[3][0];
    dest[3][1] = sourceOne[3][0] * sourceTwo[0][1] + sourceOne[3][1] * sourceTwo[1][1] + sourceOne[3][2] * sourceTwo[2][1] + sourceOne[3][3] * sourceTwo[3][1];
    dest[3][2] = sourceOne[3][0] * sourceTwo[0][2] + sourceOne[3][1] * sourceTwo[1][2] + sourceOne[3][2] * sourceTwo[2][2] + sourceOne[3][3] * sourceTwo[3][2];
    dest[3][3] = sourceOne[3][0] * sourceTwo[0][3] + sourceOne[3][1] * sourceTwo[1][3] + sourceOne[3][2] * sourceTwo[2][3] + sourceOne[3][3] * sourceTwo[3][3];
};

/**
 * Multiplies a 4x4 matrix with a 3x1 matrix.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 *   The 3x1 matrix is padded with a 4th element 1.0.
 *   The output point must have 1.0 at the 4th element. This is achieved by dividing each component of the multiplied component with the 4th element computed.
 *
 * @param dest Destination point of the multiplication operation.
 * @param sourceOne The matrix for the multiplication.
 * @param sourceTwo The point for the multiplication.
 */
void Matrix4x4MultiplyBy4x1(GzCoord dest, GzMatrix sourceOne, GzCoord sourceTwo, bool isVector)
{
    float dimension = 1.0f;

    if(isVector == false)
    {
        // GzCoord consists information of 3 vertices. Hence we do not need append the fourth point
        dimension = sourceOne[3][0] * sourceTwo[0] + sourceOne[3][1] * sourceTwo[1] + sourceOne[3][2] * sourceTwo[2] + sourceOne[3][3] * 1;
    }

    dest[0] = (sourceOne[0][0] * sourceTwo[0] + sourceOne[0][1] * sourceTwo[1] + sourceOne[0][2] * sourceTwo[2] + sourceOne[0][3] * 1) / dimension;
    dest[1] = (sourceOne[1][0] * sourceTwo[0] + sourceOne[1][1] * sourceTwo[1] + sourceOne[1][2] * sourceTwo[2] + sourceOne[1][3] * 1) / dimension;
    dest[2] = (sourceOne[2][0] * sourceTwo[0] + sourceOne[2][1] * sourceTwo[1] + sourceOne[2][2] * sourceTwo[2] + sourceOne[2][3] * 1) / dimension;
}

//TODO: Document
void MatrixMultiplyByScalar(GzMatrix matrix, float scalar)
{
    for(unsigned int column = 0; column < 4; ++column)
    {
        for(unsigned int row = 0; row < 4; ++row)
        {
            matrix[column][row] *= scalar;
        }
    }
}

//TODO: Document
float DegreeToRadian(float degree)
{
    return degree * 0.0174532925;
}

/**
 * Computes cosine and sine of a given angle
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * None.
 *
 * @param degree The angle for the computation in degrees.
 * @param sinDegree Output parameter containing the sine of the input degree angle.
 * @param sourceTwo Output parameter containing the cosine of the input degree angle.
 */
void ComputeCosAndSin(float degree, float& sinDegree, float& cosDegree)
{
    sinDegree = sin(degree * PI / 180);
    cosDegree = cos(degree * PI / 180);
}

/**
 * Computes the normal of a plane through 3 given points and also the D factor.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * None.
 *
 * @param normal Output parameter that contains the normal of the plane
 * @param planeD Output parameter that contains the D constant of the plane
 * @param threePointArray Three vertices whose plane normal and D constant is determined
 */
void computeNormal(GzCoord& normal, float& planeD, const GzCoord threePointArray[3])
{
    GzCoord firstEdge;
    firstEdge[0] = threePointArray[1][0] - threePointArray[0][0];
    firstEdge[1] = threePointArray[1][1] - threePointArray[0][1];
    firstEdge[2] = threePointArray[1][2] - threePointArray[0][2];

    GzCoord secondEdge;
    secondEdge[0] = threePointArray[2][0] - threePointArray[1][0];
    secondEdge[1] = threePointArray[2][1] - threePointArray[1][1];
    secondEdge[2] = threePointArray[2][2] - threePointArray[1][2];

    normal[0] = firstEdge[1] * secondEdge[2] - firstEdge[2] * secondEdge[1];
    normal[1] = firstEdge[2] * secondEdge[0] - firstEdge[0] * secondEdge[2];
    normal[2] = firstEdge[0] * secondEdge[1] - firstEdge[1] * secondEdge[0];

    planeD = - (normal[0] * threePointArray[0][0] + normal[1] * threePointArray[0][1] + normal[2] * threePointArray[0][2]);
}

/**
 * Returns the interpolated Z value for the given point on the triangle plane.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * Use the equation of the plane to compute Z: Ax + By + Cz + D = 0
 *
 * @param normal A normal to the triangle plane
 * @param planeD The plane constant as above
 * @param x The x coordinate of the point whose z is determined
 * @param y The y coordinate of the point whose z is determined
 * @return Z coordinate of the point
 */
int getZForVertex(GzCoord normal, float planeD, int x, int y)
{
    int z = (int) ( -(normal[0] * x + normal[1] * y + planeD) / normal[2] );
    return z;
}

/**
 * Returns the interpolated float Z value for the given point on the triangle plane.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * Use the equation of the plane to compute Z: Ax + By + Cz + D = 0
 *
 * @param normal A normal to the triangle plane
 * @param planeD The plane constant as above
 * @param x The x coordinate of the point whose z is determined
 * @param y The y coordinate of the point whose z is determined
 * @return Z coordinate of the point
 */
float getFloatZForVertex(GzCoord normal, float planeD, int x, int y)
{
    float z = -(normal[0] * x + normal[1] * y + planeD) / normal[2];
    return z;
}

/**
 * Computes the rotation of a matrix along the X axis.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * None.
 *
 * @param degree The angle for the computation of the rotation matrix in degrees.
 * @param mat Output Parameter. The computed modelling matrix.
 * @return GZ_SUCCESS value
 */
int GzRotXMat(float degree, GzMatrix mat)
{
    float cosDegree, sinDegree;
    ComputeCosAndSin(degree, sinDegree, cosDegree);

    GzMakeIdentityMatrix(mat);

    mat[1][1] = cosDegree;
    mat[2][2] = cosDegree;
    mat[1][2] = -(sinDegree);
    mat[2][1] = sinDegree;

    return GZ_SUCCESS;
}

/**
 * Computes the rotation of a matrix along the Y axis.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * None.
 *
 * @param degree The angle for the computation of the rotation matrix in degrees.
 * @param mat Output Parameter. The computed modelling matrix.
 * @return GZ_SUCCESS value
 */
int GzRotYMat(float degree, GzMatrix mat)
{
    float cosDegree, sinDegree;
    ComputeCosAndSin(degree, sinDegree, cosDegree);

    GzMakeIdentityMatrix(mat);

    mat[0][0] = cosDegree;
    mat[2][2] = cosDegree;
    mat[0][2] = sinDegree;
    mat[2][0] = -(sinDegree);

    return GZ_SUCCESS;
}

/**
 * Computes the rotation of a matrix along the Z axis.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * None.
 *
 * @param degree The angle for the computation of the rotation matrix in degrees.
 * @param mat Output Parameter. The computed modelling matrix.
 * @return GZ_SUCCESS value
 */
int GzRotZMat(float degree, GzMatrix mat)
{
    float cosDegree, sinDegree;
    ComputeCosAndSin(degree, sinDegree, cosDegree);

    GzMakeIdentityMatrix(mat);

    mat[0][0] = cosDegree;
    mat[1][1] = cosDegree;
    mat[1][0] = sinDegree;
    mat[0][1] = -(sinDegree);

    return GZ_SUCCESS;
}

/**
 * Computes the translation of a matrix for a given set of coordinates.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * None.
 *
 * @param translate The translation coordinates.
 * @param mat Output Parameter. The computed modelling matrix.
 * @return GZ_SUCCESS value
 */
int GzTrxMat(GzCoord translate, GzMatrix mat)
{
    GzMakeIdentityMatrix(mat);

    for(unsigned int index = 0; index < 3; ++index)
    {
        mat[index][3] = translate[index];
    }

    return GZ_SUCCESS;
}

/**
 * Computes the scaling matrix for a given set of coordinates.
 * <P>
 * NOTES/ASSUMPTIONS:<BR>
 * None.
 *
 * @param scaling The scaling coordinates.
 * @param mat Output Parameter. The computed modelling matrix.
 * @return GZ_SUCCESS value
 */
int GzScaleMat(GzCoord scale, GzMatrix mat)
{
    GzMakeIdentityMatrix(mat);

    for(unsigned int index = 0; index < 3; ++index)
    {
        mat[index][index] = scale[index];
    }

    return GZ_SUCCESS;
}
