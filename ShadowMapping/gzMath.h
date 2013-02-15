#ifndef __GZ_MATH_H
#define __GZ_MATH_H

#include "Gz.h"

int     GzRotXMat(float degree, GzMatrix mat);
int     GzRotYMat(float degree, GzMatrix mat);
int     GzRotZMat(float degree, GzMatrix mat);
int     GzTrxMat(GzCoord translate, GzMatrix mat);
int     GzScaleMat(GzCoord scale, GzMatrix mat);

void    CopyVector(GzCoord dest, GzCoord source);
float   VectorMagnitude(GzCoord vector);
void    NormalizeVector(GzCoord vector);
void    VectorCrossProduct(GzCoord result, GzCoord vectorOne, GzCoord vectorTwo);
float   VectorDotProduct(GzCoord vectorOne, GzCoord vectorTwo);
void    CreateVector(GzCoord vector, GzCoord headPoint, GzCoord tailPoint);
void    MultiplyScalarToVector(float scalar, GzCoord vector);

void    GzMakeIdentityMatrix(GzMatrix matrix);
void    MatrixCopy(GzMatrix dest, GzMatrix source);
void    Matrix4x4MultiplyBy4x4(GzMatrix dest, GzMatrix sourceOne, GzMatrix sourceTwo);
void    Matrix4x4MultiplyBy4x1(GzCoord dest, GzMatrix sourceOne, GzCoord sourceTwo, bool isVector = false);
void    MatrixMultiplyByScalar(GzMatrix matrix, float scalar);

float   DegreeToRadian(float degree);
void    ComputeCosAndSin(float degree, float& sinDegree, float& cosDegree);

void    computeNormal(GzCoord& normal, float& planeD, const GzCoord threePointArray[3]);
int     getZForVertex(GzCoord normal, float planeD, int x, int y);
float   getFloatZForVertex(GzCoord normal, float planeD, int x, int y);

#endif
