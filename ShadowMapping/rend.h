#include "disp.h" /* include your own disp.h file (e.g. hw1)*/

#define PI 3.14159265

/* Camera defaults */
#define	DEFAULT_FOV		45.0
#define	DEFAULT_IM_Z	(-10.0)  /* world coords for image plane origin */
#define	DEFAULT_IM_Y	(5.0)    /* default look-at point = 0,0,0 */
#define	DEFAULT_IM_X	(-10.0)

#define	DEFAULT_AMBIENT	{0.1, 0.1, 0.1}
#define	DEFAULT_DIFFUSE	{0.7, 0.6, 0.5}
#define	DEFAULT_SPECULAR	{0.2, 0.3, 0.4}
#define	DEFAULT_SPEC		32

#define	MATLEVELS	100		/* how many matrix pushes allowed */
#define	MAX_LIGHTS	10		/* how many lights allowed */

#ifndef GZRENDER
#define GZRENDER

/* define a renderer */
typedef struct
{
  GzRenderClass renderClass;
  GzDisplay     *display;
  short         open;
  GzCamera      camera;
  short         matlevel;           /* top of stack - current xform */
  GzMatrix      Ximage[MATLEVELS];  /* stack of xforms (Xsm) */
  GzMatrix      Xnorm[MATLEVELS];   /* xforms for norms (Xim) */
  GzMatrix      Xsp;                /* NDC to screen (pers-to-screen) */
  GzColor       flatcolor;          /* color state for flat shaded triangles */
  int           interp_mode;
  int           numlights;
  GzLight       lights[MAX_LIGHTS];
  GzLight       ambientlight;
  GzColor       Ka, Kd, Ks;
  float         spec;       /* specular power */
  GzTexture     tex_fun;    /* tex_fun(float u, float v, GzColor color) */

  float         aliasX;
  float         aliasY;

  GzMatrix      Xwm[MATLEVELS];
  short         matlevelsXwm;

  GzMatrix      Xlbp;
}  GzRender;
#endif

// Function declaration
int GzNewRender(GzRender** render, GzRenderClass renderClass, GzDisplay *display);
int GzFreeRender(GzRender* render);
int GzBeginRender(GzRender* render);
int GzPutAttribute(GzRender* render, int numAttributes, GzToken* nameList, GzPointer* valueList);
int GzPutTriangle(GzRender* render, int numParts, GzToken* nameList, GzPointer* valueList);

int     GzUpdateShadowBuffers(GzRender* render, int numParts, GzToken* nameList, GzPointer* valueList);
void    updateLightBuffers(GzRender* render, GzCoord triangleVertices[3]);

bool    isTopOrLeftEdge(GzCoord firstEdgePoint, GzCoord secondEdgePoint, GzCoord nonEdgePoint);
void    sortVertices(GzCoord threePointArray[3]);
int     pointOrientationToLine(GzCoord point, GzCoord headLinePoint, GzCoord tailLinePoint);
bool    isPixelShadable(GzCoord testPoint, GzCoord threePointArray[3]);
void    computeNormal(GzCoord& normal, float& planeD, const GzCoord threePointArray[3]);
int     getZForVertex(GzCoord normal, float planeD, int x, int y);
float   computeDistance(GzCoord* threePointArray, int firstVertexIndex, int secondVertexIndex);
void    getLargestEdgeVertices(GzCoord* threePointArray, int& firstVertexIndex, int& secondVertexIndex);
void    drawTriangle(GzRender* render, GzCoord firstVertex, GzCoord secondVertex, GzCoord thirdVertex, GzCoord threeNormalList[3], GzTextureIndex uvList[3]);

int GzPutCamera(GzRender *render, GzCamera *camera);
int GzPushMatrix(GzRender *render, GzMatrix matrix);
int GzPushMatrix(GzRender *render, GzMatrix matrix, bool isModelMatrix);
int GzPopMatrix(GzRender *render);

// Object Translation
int GzRotXMat(float degree, GzMatrix mat);
int GzRotYMat(float degree, GzMatrix mat);
int GzRotZMat(float degree, GzMatrix mat);
int GzTrxMat(GzCoord translate, GzMatrix mat);
int GzScaleMat(GzCoord scale, GzMatrix mat);

// Matrix functions
void GzMakeIdentityMatrix(GzMatrix matrix);

void OriginalVertexAtWorldSpace(GzRender * render, GzCoord firstVertex,GzCoord secondVertex, GzCoord thirdVertex);
void VertexAtPerpectiveSpace(GzRender * render, GzCoord ScreenSpaceVertex,GzCoord PerpVertex);
void inversePerspective(GzCoord initialVertex, GzCoord transformedVertex,float zval );
void VertexBackAtWorldSpace(GzRender * render, GzCoord ScreenSpaceVertex,GzCoord *outputVertex);

bool isAffectedByLight(GzRender* render,GzLight* light,GzCoord orivert);