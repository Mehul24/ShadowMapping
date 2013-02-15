/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
//#include<math.h>

#define zmaxbyd 677099008.00

typedef struct vertice
{
	float x,y,z;
	float nx,ny,nz;
	float u, v;
};

typedef struct dda{
	vertice start,end,current;
	float dy,dx,dz,mxy,mzy,mzx;
		
	GzColor sColor,eColor,cColor;
	GzCoord sNormal,eNormal,cNormal;
	GzCoord mcy,mcx,mny,mnx;

	float muy,mvy,mux,mvx;
	float dyp;
	bool drawline;
};


vertice StoreVertex(float x, float y, float z)
{
	vertice temp;
	temp.x=x;
	temp.y=y;
	temp.z=z;
	return temp;
}

void StoreNormal(vertice * temp, float x, float y, float z)
{
	temp->nx=x;
	temp->ny=y;
	temp->nz=z;
	//return temp;
}

void StoreUVCoord(vertice * temp, float u, float v)
{
	temp->u = u;
	temp->v = v;
}

int compare(vertice *a, vertice* b)
{
	if (a->y > b->y)
	{
		return -1;
	}
	else if (a->y < b->y)
	{
		return 1;
	}
	else 
	{
		if (a->x > b->x)
		{
			return -1;
		}
		else
		{
			return 1;
		}
	}
}

void SortVertice(vertice *a1,vertice *a2,vertice *a3)
{
	vertice a[3];
	vertice temp;

	a[0] = *a1;
	a[1] = *a2;
	a[2] = *a3;

	for(int i=0;i<3;i++)
	{
		for(int j=i+1; j<3;j++)
		{
			if(compare(&a[i], &a[j]) < 0)
			{
				temp = a[j];
				a[j] = a[i];
				a[i] = temp;
			}
		}
	}

		*a1 = a[0];
		*a2 = a[1];
		*a3 = a[2];
}


void filldda(dda *line,vertice a1,vertice a2)
{
	line->start = a1;
	line->current = line->start;
	line->end = a2;
	
	line->dy = a1.y - a2.y;
	line->dx = a1.x - a2.x;
	line->dz = a1.z - a2.z;

	if(line->dy ==  0)
	{
		line->mxy = line->mzy=0;
	}
	else
	{
	line->mxy = line->dx / line->dy;
	line->mzy = line->dz / line->dy ;
	}
	line->mzx = 0;
}

void fillddaColor(dda *line, GzColor c1, GzColor c2)
{
	memcpy(line->sColor,c1,sizeof(GzColor));
	memcpy(line->eColor,c2,sizeof(GzColor));
	memcpy(line->cColor,c1,sizeof(GzColor));

	if(line->dy == 0)
	{
		line->mcy[0]= line->mcy[1]=line->mcy[2]=0;
	}
	else
	{
		line->mcy[0] = (line->sColor[0]-line->eColor[0])/line->dy;
		line->mcy[1] = (line->sColor[1]-line->eColor[1])/line->dy;
		line->mcy[2] = (line->sColor[2]-line->eColor[2])/line->dy;
	}
	line->mcx[0] = line->mcx[1] = line->mcx[2] = 0;
}

void fillddaNormal(dda *line, GzCoord N1, GzCoord N2)
{
	memcpy(line->sNormal,N1,sizeof(GzCoord));
	memcpy(line->cNormal,N1,sizeof(GzCoord));
	memcpy(line->eNormal,N2,sizeof(GzCoord));

	if(line->dy == 0)
	{
		line->mny[0]= line->mny[1]=line->mny[2]=0;
	}
	else
	{
		line->mny[0] = (line->sNormal[0] - line->eNormal[0])/line->dy;
		line->mny[1] = (line->sNormal[1] - line->eNormal[1])/line->dy;
		line->mny[2] = (line->sNormal[2] - line->eNormal[2])/line->dy;
	}
	line->mnx[0] = line->mnx[1] = line->mnx[2] = 0;
}

void addUVCoord(dda*l , vertice a1, vertice a2)
{
	l->start.u = l->current.u = a1.u;
	l->end.u = a2.u;

	l->start.v = l->current.v = a1.v;
	l->end.v = a2.v;

	if(l->dy == 0)
		l->muy = l->mvy = 0;
	else
	{	
		
		l->muy = (l->start.u - l->end.u)/l->dy;
		l->mvy = (l->start.v - l->end.v)/l->dy;

	}
	l->mux = l->mvx =0;
}

/* NOT part of API - just for general assistance */

short	ctoi(float color)		/* convert float color to GzIntensity short */
{
  return(short)((int)(color * ((1 << 12) - 1)));
}

float degreeToRad(float degree)
{
	return degree*0.0174532925;
}

int GzRotXMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along x axis
// Pass back the matrix using mat value
	float Sin = sin(degree*0.0174532);
	float Cos = cos(degree*0.0174532);

	GzMatrix Rotx =
	{
		1, 0, 0, 0,
		0, Cos, -Sin, 0,
		0, Sin, Cos, 0,
		0, 0, 0, 1
	};

	memcpy(mat, Rotx, sizeof(GzMatrix));

	return GZ_SUCCESS;
}


int GzRotYMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along y axis
// Pass back the matrix using mat value

	float Sin = sin(degree*0.0174532);
	float Cos = cos(degree*0.0174532);

	GzMatrix Roty =
	{
		Cos, 0, Sin, 0,
		0, 1, 0, 0,
		-Sin, 0, Cos, 0,
		0, 0, 0, 1
	};

	memcpy(mat, Roty, sizeof(GzMatrix));

	return GZ_SUCCESS;
}


int GzRotZMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along z axis
// Pass back the matrix using mat value

	float Sin = sin(degree*0.0174532);
	float Cos = cos(degree*0.0174532);

	GzMatrix Rotz =
	{
		Cos, -Sin, 0, 0,
		Sin, Cos, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1
	};

	memcpy(mat, Rotz, sizeof(GzMatrix));

	return GZ_SUCCESS;
}


int GzTrxMat(GzCoord translate, GzMatrix mat)
{
// Create translation matrix
// Pass back the matrix using mat value
	GzMatrix trans =
	{
		1, 0, 0, translate[0],
		0, 1, 0, translate[1],
		0, 0, 1, translate[2],
		0, 0, 0, 1
	};

	memcpy(mat, trans, sizeof(GzMatrix));
	
	return GZ_SUCCESS;
}


int GzScaleMat(GzCoord scale, GzMatrix mat)
{
// Create scaling matrix
// Pass back the matrix using mat value

	GzMatrix sc =
	{
		scale[0], 0, 0, 0,
		0, scale[1], 0, 0,
		0, 0, scale[2], 0,
		0, 0, 0, 1
	};

	memcpy(mat, sc, sizeof(GzMatrix));

	return GZ_SUCCESS;
}


//----------------------------------------------------------
// Begin main functions

int GzNewRender(GzRender **render, GzRenderClass renderClass, GzDisplay	*display)
{
/*  
- malloc a renderer struct 
- keep closed until all inits are done 
- setup Xsp and anything only done once 
- span interpolator needs pointer to display 
- check for legal class GZ_Z_BUFFER_RENDER 
- init default camera 
*/ 
	(*render) = new GzRender();
	(*render)->open = false;

	(*render)->renderClass = renderClass;
	(*render)->display = display;

	(*render)->camera.FOV = DEFAULT_FOV;
	(*render)->camera.position[0] = DEFAULT_IM_X;
	(*render)->camera.position[1] = DEFAULT_IM_Y;
	(*render)->camera.position[2] = DEFAULT_IM_Z;
	(*render)->camera.lookat[0] = 0;
	(*render)->camera.lookat[1] = 0;
	(*render)->camera.lookat[2] = 0;
	(*render)->camera.worldup[0] = 0;
	(*render)->camera.worldup[1] = 1;
	(*render)->camera.worldup[2] = 0;
	(*render)->matlevel = -1;

	float d = 1/tan(degreeToRad((*render)->camera.FOV/2));
	float zmax= MAXINT;
		
	GzMatrix Xsp1 = 
	{
		display->xres/2, 0, 0, display->xres/2,
		0, -display->yres/2,0, display->yres/2,
		0, 0,zmax/d , 0,
		0, 0, 0, 1
	};
	memcpy((*render)->Xsp, Xsp1,sizeof(Xsp1));


	if(!(*render)->renderClass == GZ_Z_BUFFER_RENDER)
		return GZ_FAILURE;
	return GZ_SUCCESS;

}


int GzFreeRender(GzRender *render)
{
/* 
-free all renderer resources
*/

	free(render);
	return GZ_SUCCESS;
}


int GzBeginRender(GzRender *render)
{
/*  
- set up for start of each frame - clear frame buffer 
- compute Xiw and projection xform Xpi from camera definition 
- init Ximage - put Xsp at base of stack, push on Xpi and Xiw 
- now stack contains Xsw and app can push model Xforms if it want to. 
*/ 

	float d = 1/tan((render->camera.FOV*0.01745329)/2);
	
	GzMatrix Xpi1 = 
	{
		1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		0, 0, 1/d, 1
	};
	memcpy((render)->camera.Xpi , Xpi1, sizeof(Xpi1));

	GzCoord Cl;
	GzCoord CamZ;
	Cl[0] = (render)->camera.lookat[0] - (render)->camera.position[0];
	Cl[1] = (render)->camera.lookat[1] - (render)->camera.position[1];
	Cl[2] = (render)->camera.lookat[2] - (render)->camera.position[2];

	float modCl = sqrt(Cl[0]*Cl[0] + Cl[1]*Cl[1] + Cl[2]*Cl[2]);
		
	CamZ[0] = Cl[0]/modCl;
	CamZ[1] = Cl[1]/modCl;
	CamZ[2] = Cl[2]/modCl;

	GzCoord CamY,Up1;
	float modUp1,UpZ;
	UpZ = (render)->camera.worldup[0] *CamZ[0] +  (render)->camera.worldup[1] *CamZ[1] +  (render)->camera.worldup[2] *CamZ[2];
	Up1[0] = (render)->camera.worldup[0] - UpZ*CamZ[0];
	Up1[1] = (render)->camera.worldup[1] - UpZ*CamZ[1];
	Up1[2] = (render)->camera.worldup[2] - UpZ*CamZ[2];

	modUp1 = sqrt(Up1[0]*Up1[0] + Up1[1]*Up1[1] + Up1[2]*Up1[2]);
	CamY[0] = Up1[0]/modUp1;
	CamY[1] = Up1[1]/modUp1;
	CamY[2] = Up1[2]/modUp1;

	GzCoord CamX;
	CamX[0] = CamY[1]*CamZ[2] - CamY[2]*CamZ[1];
	CamX[1] = CamY[2]*CamZ[0] - CamY[0]*CamZ[2];
	CamX[2] = CamY[0]*CamZ[1] - CamY[1]*CamZ[0];

	float xc,yc,zc;
	xc = (render)->camera.position[0]*CamX[0] +(render)->camera.position[1]*CamX[1]+ (render)->camera.position[2]*CamX[2];
	yc = (render)->camera.position[0]*CamY[0] +(render)->camera.position[1]*CamY[1]+ (render)->camera.position[2]*CamY[2];
	zc = (render)->camera.position[0]*CamZ[0] +(render)->camera.position[1]*CamZ[1]+ (render)->camera.position[2]*CamZ[2];

	GzMatrix Xiw1 =
	{
		CamX[0], CamX[1], CamX[2],-xc,
		CamY[0], CamY[1], CamY[2],-yc ,
		CamZ[0], CamZ[1], CamZ[2],-zc ,
		0, 0, 0, 1
	};

	memcpy(render->camera.Xiw,Xiw1,sizeof(Xiw1));

	GzPushMatrix(render, render->Xsp);
	GzPushMatrix(render,render->camera.Xpi);
	GzPushMatrix(render,render->camera.Xiw);
	
	render->open = true;
	return GZ_SUCCESS;
}

int GzPutCamera(GzRender *render, GzCamera *camera)
{
/*
- overwrite renderer camera structure with new camera definition
*/	render->camera = *camera;

	float d = 1/tan(degreeToRad((render)->camera.FOV/2));
	float zmax= MAXINT;
		
	GzMatrix Xsp1 = 
	{
		render->display->xres/2, 0, 0, render->display->xres/2,
		0, -render->display->yres/2,0, render->display->yres/2,
		0, 0,zmax/d , 0,
		0, 0, 0, 1
	};
	memcpy((render)->Xsp, Xsp1,sizeof(Xsp1));
	

	return GZ_SUCCESS;	
}



void MultiplyMatrix(GzMatrix a, GzMatrix b, GzMatrix output)
{
	int i,j,k;
	for (i=0; i<4; i++) 
	{
		for (k=0; k<4; k++) 
		{
			output[i][k]=0;
			for (j=0; j<4; j++) 
			{
				output[i][k]+=a[i][j]*b[j][k]; 
			}
		}
	}
}

void MultiplyMatrixVertex(GzMatrix a, vertice *outputVertex)
{
	int i ,j ;
	float output[4];
	float b[4];
	b[0] = outputVertex->x;
	b[1] = outputVertex->y;
	b[2] = outputVertex->z;
	b[3] = 1;

	for (i=0; i<4; i++) 
	{
			output[i]=0;
			for (j=0; j<4; j++) 
			{
				output[i]+=a[i][j]*b[j]; 
			}
	}

	outputVertex->x = (output[0] /output[3]);
	outputVertex->y = (output[1] /output[3]);
	outputVertex->z = (output[2] /output[3]);
}

void MultiplyMatrixNormal(GzMatrix a, vertice *outputVertex)
{
	int i ,j ;
	float output[4];
	float b[4];
	b[0] = outputVertex->nx;
	b[1] = outputVertex->ny;
	b[2] = outputVertex->nz;
	b[3] = 1;
	for (i=0; i<4; i++) 
	{
			output[i]=0;
			for (j=0; j<4; j++) 
			{
				output[i]+=a[i][j]*b[j]; 
			}
	}

	outputVertex->nx = (output[0] /output[3]);
	outputVertex->ny = (output[1] /output[3]);
	outputVertex->nz = (output[2] /output[3]);
}


void MakeIdentityNormal(GzMatrix a)
{
	GzMatrix I={
		1,0,0,0,
		0,1,0,0,
		0,0,1,0,
		0,0,0,1
	};
	memcpy(a,I,sizeof(GzMatrix));

}

int GzPushMatrix(GzRender *render, GzMatrix	matrix)
{
/*
- push a matrix onto the Ximage stack
- check for stack overflow
*/
	GzMatrix out,out1;
	GzMatrix normMat;

	if(render->matlevel == -1)
	{
		render->matlevel++;
		memcpy(render->Ximage[render->matlevel], matrix,sizeof(GzMatrix));
		MakeIdentityNormal(render->Xnorm[render->matlevel]);

	}
	else
	{
		render->matlevel++;
		
		MultiplyMatrix(render->Ximage[render->matlevel-1], matrix, out);
		memcpy(render->Ximage[render->matlevel], out,sizeof(GzMatrix));
		
		//SKIP FOR XPI MATRIX.	Explicitlly pushing XPI as Identity
		if(render->matlevel == 1)
			MakeIdentityNormal(render->Xnorm[render->matlevel]);
		else
		{
			memcpy(normMat,matrix,sizeof(GzMatrix));
			normMat[0][3] = normMat[1][3] = normMat[2][3] = 0 ;
			normMat[3][0] = normMat[3][1] = normMat[3][2] = 0 ;

			float a = normMat[0][0];
			float b = normMat[0][1];
			float c	 = normMat[0][2];
			float k = sqrt (a*a + b*b + c*c);
			for (int i=0;i<4;i++)
			{
				for (int j=0;j<4;j++)
					normMat[i][j] /= k;
			}

			MultiplyMatrix(render->Xnorm[render->matlevel-1], normMat, out1);
			memcpy(render->Xnorm[render->matlevel], out1,sizeof(GzMatrix));
		}
	}

	return GZ_SUCCESS;
}

int GzPopMatrix(GzRender *render)
{
/*
- pop a matrix off the Ximage stack
- check for stack underflow
*/
	if(render->matlevel == -1 )
		return GZ_FAILURE;
	render->matlevel--;

	return GZ_SUCCESS;
}


int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList, 
	GzPointer	*valueList) /* void** valuelist */
{
/*
- set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
- later set shaders, interpolaters, texture maps, and lights
*/
	int index;
	GzLight ambLig;
	for(index=0;index<numAttributes;index++)
	{
		switch(nameList[index])
		{
			case GZ_RGB_COLOR:
				render->flatcolor[0] = ((float*)valueList[index])[0];
				render->flatcolor[1] = ((float*)valueList[index])[1];
				render->flatcolor[2] = ((float*)valueList[index])[2];
				break;
		
			case GZ_AMBIENT_LIGHT:
				//ambiLig = 
				render->ambientlight = *((GzLight*)valueList[index]);
				break;

			case GZ_DIRECTIONAL_LIGHT:
				render->lights[render->numlights] = *((GzLight*)valueList[index]);
				render->numlights++;
				break;
			
			case GZ_DIFFUSE_COEFFICIENT:
				render->Kd[0] = ((float*)valueList[index])[0];
				render->Kd[1] = ((float*)valueList[index])[1];
				render->Kd[2] = ((float*)valueList[index])[2];
				break;
			
			case GZ_AMBIENT_COEFFICIENT:
				render->Ka[0] = ((float*)valueList[index])[0];
				render->Ka[1] = ((float*)valueList[index])[1];
				render->Ka[2] = ((float*)valueList[index])[2];
				break;
			
			case GZ_SPECULAR_COEFFICIENT:
				render->Ks[0] = ((float*)valueList[index])[0];
				render->Ks[1] = ((float*)valueList[index])[1];
				render->Ks[2] = ((float*)valueList[index])[2];
				break;

			case  GZ_DISTRIBUTION_COEFFICIENT:
				render->spec = ((float*)valueList[index])[0];
				break;

			case GZ_INTERPOLATE:
				render->interp_mode = *((int*)valueList[index]);		//(int*) typecasts a void pointer to int *. an * before it assign the value to an int
				break;

			case GZ_TEXTURE_MAP:
				render->tex_fun = (GzTexture)valueList[index];
				break;

			case GZ_AASHIFTX:
				render->aashiftX = *(float*)valueList[index];
				break;
			
			case GZ_AASHIFTY:
				render->aashiftY = *(float*)valueList[index];
				break;
		
		}
	}
	
	return GZ_SUCCESS;
}

void interpolateColorAlongY(dda *l,float dely)
{
	l->cColor[0] = l->cColor[0] + l->mcy[0]*dely;
	l->cColor[1] = l->cColor[1] + l->mcy[1]*dely;
	l->cColor[2] = l->cColor[2] + l->mcy[2]*dely;

}

void interpolateColorAlongX(dda * l,float delx)
{
	l->cColor[0] = l->cColor[0] + l->mcx[0]*delx;
	l->cColor[1] = l->cColor[1] + l->mcx[1]*delx;
	l->cColor[2] = l->cColor[2] + l->mcx[2]*delx;

}

void interpolateNormalAlongY(dda *l,float dely)
{
	l->cNormal[0] = l->cNormal[0] + l->mny[0]*dely;
	l->cNormal[1] = l->cNormal[1] + l->mny[1]*dely;
	l->cNormal[2] = l->cNormal[2] + l->mny[2]*dely;

}

void interpolateNormalAlongX(dda * l,float delx)
{
	l->cNormal[0] = l->cNormal[0] + l->mnx[0]*delx;
	l->cNormal[1] = l->cNormal[1] + l->mnx[1]*delx;
	l->cNormal[2] = l->cNormal[2] + l->mnx[2]*delx;
}


void interpolateUVAlongY(dda *l,float dely)
{
	l->current.u = l->current.u + l->muy*dely;
	l->current.v = l->current.v + l->mvy*dely;
}

void interpolateUVAlongX(dda *l,float delx)
{
	l->current.u = l->current.u + l->mux*delx;
	l->current.v = l->current.v + l->mvx*delx;
}

void CreateNormalVector(vertice vert, float *normal)
{
	float length = sqrt(vert.nx * vert.nx + vert.ny * vert.ny + vert.nz * vert.nz);
	normal[0] = vert.nx / length;
	normal[1] = vert.ny / length;
	normal[2] = vert.nz / length;
}

float DotProduct(GzCoord a, GzCoord b)
{
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

float MagnitudeOfVector(GzCoord a)
{
	return sqrt(DotProduct(a, a));
}

void NormalizeVector(float *coord)
{
	float mag = sqrt(coord[0] * coord[0] + coord[1] * coord[1] + coord[2] * coord[2]);
	coord[0] /= mag;
	coord[1] /= mag;
	coord[2] /= mag;
}

void CopyGzCoord(float *newCoord, float *oldCoord)
{
	newCoord[0] = oldCoord[0];
	newCoord[1] = oldCoord[1];
	newCoord[2] = oldCoord[2];
}

void CopyGzColor(float *toColor, float *fromColor)
{
	toColor[0] = fromColor[0];
	toColor[1] = fromColor[1];
	toColor[2] = fromColor[2];
}

void Normalize(float *vector)
{
	float magni = sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]);
	vector[0] /=magni;
	vector[1] /=magni;
	vector[2] /=magni;

}

void ComputeR(GzCoord R,GzCoord Normal,GzLight Light)
{
	float nl2 = 2*DotProduct(Normal,Light.direction);
	GzCoord n;
	memcpy(n,Normal,sizeof(GzCoord));
	for(int i =0;i<3;i++)
		n[i] = n[i] *nl2;
	
	GzCoord R1;
	for(int i =0;i<3;i++)
		R1[i] = n[i] - Light.direction[i];

	memcpy(R,R1,sizeof(GzCoord));
}


void ComputeColor(GzRender *render, GzColor *color, GzCoord *normal)
{
//Ks * sumn (le(R.E)^s) 
//	+ Kd * sumn (le(N.L))
//	+ ka*la
	GzCoord eye = {0,0,-1};
	GzColor leres = {0, 0, 0};
	GzColor dl = {0, 0, 0};
	GzCoord nNormal = {0, 0, 0};
	CopyGzCoord(nNormal, (float*)normal);
	NormalizeVector(nNormal);
	for(int i=0;i<render->numlights;i++)
	{
		NormalizeVector((float*)render->lights[i].direction);
		float nl = DotProduct(nNormal, render->lights[i].direction);
		float ne = DotProduct(nNormal, eye);
		if(nl < 0 && ne < 0)
		{
			float *normalFloat = (float*)nNormal;
			for(int j=0;j<3;j++)
				normalFloat[j] *= -1;
		}
		if((nl < 0 && ne < 0) || (nl > 0 && ne > 0))
		{
			GzCoord rCoord;
			ComputeR(rCoord,nNormal, render->lights[i]);
			NormalizeVector((float*)rCoord);

			float rePowerS = pow(DotProduct(rCoord, eye), render->spec);
			for(int j=0;j<3;j++)
				leres[j] += render->lights[i].color[j] * rePowerS;
			
			float nl = DotProduct(nNormal, render->lights[i].direction);
			for(int j=0;j<3;j++)
				dl[j] += render->lights[i].color[j] * nl;
			
		}
	}

	for(int i=0;i<3;i++)
		leres[i] *= render->Ks[i];
	
	for(int i=0;i<3;i++)
		dl[i] *= render->Kd[i];
	
	GzColor al;
	for(int i=0;i<3;i++)
		al[i] = render->Ka[i] * render->ambientlight.color[i];
	
	float *finalC = (float*)color;
	for(int i=0;i<3;i++)
		finalC[i] = leres[i] + dl[i] + al[i];
	
	finalC[0] = min(1.0, leres[0] + dl[0] + al[0]);
	finalC[1] = min(1.0, leres[1] + dl[1] + al[1]);
	finalC[2] = min(1.0, leres[2] + dl[2] + al[2]);

}



/*
void ComputeColor(GzRender *render, GzColor *color, GzCoord *normal)
{
//Ks * sumn (le(R.E)^s) 
//	+ Kd * sumn (le(N.L))
//	+ ka*la

	GzColor dL={0,0,0},sL={0,0,0},aL={0};
	GzCoord eye = {0,0,-1};
	GzCoord repows= {0,0,0};
	GzCoord nNormal = {0};
	
	memcpy(nNormal,*normal,sizeof(GzCoord));
	Normalize(nNormal);

	for(int noL=0;noL<render->numlights;noL++)
	{
		Normalize(render->lights[noL].direction);
		float nl,ne;
		nl = DotProduct(nNormal,render->lights[noL].direction);
		ne = DotProduct(nNormal,eye);

		if(nl<0 && ne <0)
		{
			nNormal[0] *= -1;
			nNormal[1] *= -1;
			nNormal[2] *= -1;
		}
		
		if((nl>0 && ne >0) || (nl<0 && ne <0))
		{
			GzCoord R;
			ComputeR(R,nNormal,render->lights[noL]);
			Normalize(R);
			
			float  repows = pow(DotProduct(R,eye),render->spec);
			for(int i =0;i<3;i++)
				sL[i] += render->lights[noL].color[i] * repows;
			
			for(int i =0;i<3;i++)
				dL[i] += render->lights[noL].color[i] * nl;
		}

	}
	for(int i =0;i<3;i++)
		sL[i] *=render->Ks[i];
		
	for(int i =0;i<3;i++)
		dL[i] *=render->Kd[i];
		
	for(int i =0;i<3;i++)
		aL[i] = render->Ka[i] * render->ambientlight.color[i];

	
	float *finalColor = (float*)color;
	for(int i =0;i<3;i++)
		finalColor[i] = sL[i] + dL[i] + aL[i];
	
	finalColor[0] = min(1.0, finalColor[0]);
	finalColor[1] = min(1.0, finalColor[0]);
	finalColor[2] = min(1.0, finalColor[0]);

}
*/


int pos(float a)
{
	if (a<0)
		return -a;
	return a;
}


void PerspectiveInterpolation(GzRender * render,vertice *a)
{
	float zmax = MAXINT;
	float zp;

	zp = a->z/(zmax - a->z);
	a->u = a->u/(zp+1);
	a->v = a->v/(zp+1);
	
}

void ReversePerspectiveInterpolation(GzRender * render, dda * l,float *u,float *v)
{
	float zmax = MAXINT;

	float zp;	
	zp = l->current.z/(zmax - l->current.z);
	*u = l->current.u *(zp+1);
	*v = l->current.v *(zp+1);
}



void GouraudRasterize(GzRender	*render,vertice a1, vertice a2, vertice a3, dda l12, dda l13,dda l23)
{
	int dir=0, j,flag=0;
	float dely , delx;
	GzIntensity r,g,b,a;
	GzDepth z;
	bool horzLine = false;
	dda l;
	GzColor currentColor;

	if(a1.y == a2.y || a2.y ==  a3.y )
			horzLine = true;

	if(l12.mxy < l13.mxy)
	{
		l12.drawline = true;
		l23.drawline = true;
		l13.drawline = false;
	}
	else
	{
		l12.drawline = false;
		l23.drawline = false;
		l13.drawline = true;
	}

		//i = ceil(l13.current.y);
	dely = ceil(l13.current.y) - l13.current.y;
	l13.current.x = l13.current.x + l13.mxy*dely;
	l13.current.y = l13.current.y +dely;
	l13.current.z = l13.current.z +l13.mzy*dely;

	l12.current.x = l12.current.x + l12.mxy*dely;
	l12.current.y = l12.current.y +dely;
	l12.current.z = l12.current.z +l12.mzy*dely;

	//Gouraud
	if(render->interp_mode == GZ_COLOR)
	{
		interpolateColorAlongY(&l13,dely);
		interpolateColorAlongY(&l12,dely);
	}

	//TEXTURE INTERPOLATION
	interpolateUVAlongY(&l13,dely);
	interpolateUVAlongY(&l12,dely);

	while(l13.current.y < l12.end.y)
	//for(i; i<(l12.end.y) ; i++)
	{
		if(l12.drawline == true)
		{
			if(render->interp_mode == GZ_COLOR)
			{
				memcpy(l.sColor,l12.cColor,sizeof(GzColor));
				memcpy(l.cColor,l12.cColor,sizeof(GzColor));
				memcpy(l.eColor,l13.cColor,sizeof(GzColor));
			}

			l.start.x = l.current.x = l12.current.x ;
			l.end.x = l13.current.x ;
			l.start.z = l.current.z = l12.current.z ;
			l.end.z = l13.current.z ;

			//Texture
			l.start.u = l.current.u = l12.current.u;
			l.end.u = l13.current.u;
			l.start.v = l.current.v = l12.current.v;
			l.end.v = l13.current.v;

		}
		else
		{
			if(render->interp_mode == GZ_COLOR)
			{
				memcpy(l.sColor,l13.cColor,sizeof(GzColor));
				memcpy(l.cColor,l13.cColor,sizeof(GzColor));
				memcpy(l.eColor,l12.cColor,sizeof(GzColor));
			}

			l.start.x = l.current.x = l13.current.x ;
			l.end.x = l12.current.x ;
			l.start.z = l.current.z = l13.current.z ;
			l.end.z = l12.current.z ;

			//Texture
			l.start.u = l.current.u = l13.current.u;
			l.end.u = l12.current.u;
			l.start.v = l.current.v = l13.current.v;
			l.end.v = l12.current.v;
		}

		if((l.start.x - l.end.x)!= 0)
		{
			
			if(render->interp_mode == GZ_COLOR)
			{
				l.mcx[0] = (l.eColor[0] - l.sColor[0])/(l.end.x - l.start.x ) ;
				l.mcx[1] = (l.eColor[1] - l.sColor[1])/(l.end.x - l.start.x ) ;
				l.mcx[2] = (l.eColor[2] - l.sColor[2])/(l.end.x - l.start.x ) ;
			}

			//TEXTURE
			l.mux = (l.start.u - l.end.u) /(l.start.x - l.end.x);
			l.mvx = (l.start.v - l.end.v) /(l.start.x - l.end.x);

			l.mzx = (l.start.z - l.end.z)/ (l.start.x - l.end.x) ;
				
		}
		else
		{
			l.mzx = 0;
			l.mcx[0] = l.mcx[1] = l.mcx[2] = 0;
			l.mux = l.mvx = 0;
		}


		delx = ceil(l.start.x) - l.start.x ;
		l.current.z = l.current.z + l.mzx * delx;
		l.current.x = l.current.x + delx;

		//Gouraud
		if(render->interp_mode == GZ_COLOR)
			interpolateColorAlongX(&l,delx);

		//TEXTURE
		interpolateUVAlongX(&l,delx);
	
		for(j=ceil(l.start.x); j<(l.end.x) ;j++)
		{
				
			GzColor C,textColor;
			float u,v;

			memcpy(render->flatcolor,l.cColor,sizeof(GzColor));
			
			ReversePerspectiveInterpolation(render,&l,&u,&v);
			render->tex_fun(u,v,textColor);
			for(int i=0;i<3;i++)
				render->flatcolor[i] *=textColor[i];
				
			GzGetDisplay(render->display,j,(int)l13.current.y,&r,&g,&b,&a,&z);
			if((l.start.x == l.current.x) && (l12.drawline == true))
			{
				j++;
				continue;
			}
			else if(l.current.z <=z)
			{
				z = ceil(l.current.z);
				GzPutDisplay(render->display,j,(int)l13.current.y,ctoi(render->flatcolor[0]),ctoi(render->flatcolor[1]),ctoi(render->flatcolor[2]),1,z);
			}
			delx=1;
			l.current.z = l.current.z + l.mzx * delx;
			l.current.x = l.current.x + delx;
				
			//Gouraud
			if(render->interp_mode == GZ_COLOR)
				interpolateColorAlongX(&l,delx);

			//TEXTURE
			interpolateUVAlongX(&l,delx);
		}

		dely = 1 ;
		l13.current.x = l13.current.x + l13.mxy*dely;
		l13.current.y = l13.current.y +dely;
		l13.current.z = l13.current.z +l13.mzy*dely;

		l12.current.x = l12.current.x + l12.mxy*dely;
		l12.current.y = l12.current.y +dely;
		l12.current.z = l12.current.z +l12.mzy*dely;

		if(render->interp_mode == GZ_COLOR)
		{
			//Gouraud
			interpolateColorAlongY(&l12,dely);
			interpolateColorAlongY(&l13,dely);
		}

		//TEXTURE
		interpolateUVAlongY(&l13,dely);
		interpolateUVAlongY(&l12,dely);

	}	//END OF WHILE
		
	//i= ceil(l23.start.y);
	dely = l13.current.y - l23.current.y;
	l23.current.x = l23.current.x + l23.mxy*dely;
	l23.current.y = l23.current.y +dely;
	l23.current.z = l23.current.z +l23.mzy*dely;

	//Gouraud
	if(render->interp_mode == GZ_COLOR)
		interpolateColorAlongY(&l23,dely);

	//TEXTURE
	interpolateUVAlongY(&l23,dely);

	//for(i; i<ceil(l13.end.y) ; i++)
	while(l13.current.y < l23.end.y)
	{
	
		if(l23.drawline == true)
		{
			if(render->interp_mode == GZ_COLOR)
			{
				memcpy(l.sColor,l23.cColor,sizeof(GzColor));
				memcpy(l.cColor,l23.cColor,sizeof(GzColor));
				memcpy(l.eColor,l13.cColor,sizeof(GzColor));
			}

			l.start.x = l.current.x = l23.current.x ;
			l.end.x = l13.current.x ;
			l.start.z = l.current.z = l23.current.z ;
			l.end.z = l13.current.z ;

			//Texture
			l.start.u = l.current.u = l23.current.u;
			l.end.u = l13.current.u;
			l.start.v = l.current.v = l23.current.v;
			l.end.v = l13.current.v;	
			
		}
		else
		{
			if(render->interp_mode == GZ_COLOR)
			{
				memcpy(l.sColor,l13.cColor,sizeof(GzColor));
				memcpy(l.cColor,l13.cColor,sizeof(GzColor));
				memcpy(l.eColor,l23.cColor,sizeof(GzColor));
			}

			l.start.x = l.current.x = l13.current.x ;
			l.end.x = l23.current.x ;
			l.start.z = l.current.z = l13.current.z ;
			l.end.z = l23.current.z ;

			//Texture
			l.start.u = l.current.u = l13.current.u;
			l.end.u = l23.current.u;
			l.start.v = l.current.v = l13.current.v;
			l.end.v = l23.current.v;
			
		}

		if((l.start.x - l.end.x)!= 0)
		{
			
			if(render->interp_mode == GZ_COLOR)
			{
				l.mcx[0] = (l.eColor[0] - l.sColor[0])/(l.end.x - l.start.x ) ; ;
				l.mcx[1] = (l.eColor[1] - l.sColor[1])/(l.end.x - l.start.x ) ; ;
				l.mcx[2] = (l.eColor[2] - l.sColor[2])/(l.end.x - l.start.x ) ; ;
			}

			//TEXTURE
			l.mux = (l.start.u - l.end.u) /(l.start.x - l.end.x);
			l.mvx = (l.start.v - l.end.v) /(l.start.x - l.end.x);
			
			l.mzx = (l.start.z - l.end.z)/ (l.start.x - l.end.x) ;
		}
		else
		{
			l.mzx = 0;
			l.mcx[0] = l.mcx[1] = l.mcx[2] = 0;
		}

		delx = ceil(l.start.x) - l.start.x ;
		l.current.x = l.current.x + delx;
		l.current.z = l.current.z + l.mzx * delx;

		if(render->interp_mode == GZ_COLOR)
			interpolateColorAlongX(&l,delx);

		//Texture
		interpolateUVAlongX(&l,delx);

		for(j=ceil(l.start.x); j<(l.end.x) ;j++)
		{
			GzColor C,textColor;
			float u,v;

			memcpy(render->flatcolor,l.cColor,sizeof(GzColor));
			
			ReversePerspectiveInterpolation(render,&l,&u,&v);
			render->tex_fun(u,v,textColor);
			for(int i=0;i<3;i++)
				render->flatcolor[i] *=textColor[i];
			
			GzGetDisplay(render->display,j,(int)l13.current.y,&r,&g,&b,&a,&z);
			if(l.start.x == l.current.x && (l12.drawline == true))
			{
				j++;
				continue;
			}
			else if(l.current.z <=z)
			{
				z = (GzDepth)ceil(l.current.z);
				GzPutDisplay(render->display,j,(int)l13.current.y,ctoi(render->flatcolor[0]),ctoi(render->flatcolor[1]),ctoi(render->flatcolor[2]),1,z);
			}
			delx=1;
			l.current.x = l.current.x + delx;
			l.current.z = l.current.z + l.mzx * delx;
			
			if(render->interp_mode == GZ_COLOR)
				interpolateColorAlongX(&l,delx);

			//Texture
			interpolateUVAlongX(&l,delx);
						
		}

		dely = 1 ;
		l13.current.x = l13.current.x + l13.mxy*dely;
		l13.current.y = l13.current.y + dely ;
		l13.current.z = l13.current.z +l13.mzy*dely;

		l23.current.x = l23.current.x + l23.mxy*dely;
		l23.current.y = l23.current.y +dely;
		l23.current.z = l23.current.z +l23.mzy*dely;

		if(render->interp_mode == GZ_COLOR)
		{
			interpolateColorAlongY(&l13,dely);
			interpolateColorAlongY(&l23,dely);
		}

		//TEXTURE
		interpolateUVAlongY(&l13,dely);
		interpolateUVAlongY(&l23,dely);

	}
}


void PhongRasterize(GzRender*render,vertice a1, vertice a2, vertice a3, dda l12, dda l13,dda l23)
{
	int dir=0, j,flag=0;
	float dely , delx;
	GzIntensity r,g,b,a;
	GzDepth z;
	bool horzLine = false;
	dda l;
	GzColor currentColor;

	if(a1.y == a2.y || a2.y ==  a3.y )
			horzLine = true;

	if(l12.mxy < l13.mxy)
	{
		l12.drawline = true;
		l23.drawline = true;
		l13.drawline = false;
	}
	else
	{
		l12.drawline = false;
		l23.drawline = false;
		l13.drawline = true;
	}

		//i = ceil(l13.current.y);
	dely = ceil(l13.current.y) - l13.current.y;
	l13.current.x = l13.current.x + l13.mxy*dely;
	l13.current.y = l13.current.y +dely;
	l13.current.z = l13.current.z +l13.mzy*dely;

	l12.current.x = l12.current.x + l12.mxy*dely;
	l12.current.y = l12.current.y +dely;
	l12.current.z = l12.current.z +l12.mzy*dely;

	//PHONG
	if(render->interp_mode == GZ_NORMALS)
	{
		interpolateNormalAlongY(&l13,dely);
		interpolateNormalAlongY(&l12,dely);
	}

	//TEXTURE INTERPOLATION
	
	interpolateUVAlongY(&l13,dely);
	interpolateUVAlongY(&l12,dely);

	while(l13.current.y < l12.end.y)
	//for(i; i<(l12.end.y) ; i++)
	{
		if(l12.drawline == true)
		{
			if(render->interp_mode == GZ_NORMALS)
			{
				memcpy(l.sNormal,l12.cNormal,sizeof(GzCoord));
				memcpy(l.cNormal,l12.cNormal,sizeof(GzCoord));
				memcpy(l.eNormal,l13.cNormal,sizeof(GzCoord));
			}

			l.start.x = l.current.x = l12.current.x ;
			l.end.x = l13.current.x ;
			l.start.z = l.current.z = l12.current.z ;
			l.end.z = l13.current.z ;

			//Texture
			l.start.u = l.current.u = l12.current.u;
			l.end.u = l13.current.u;
			l.start.v = l.current.v = l12.current.v;
			l.end.v = l13.current.v;

		}
		else
		{
			if(render->interp_mode == GZ_NORMALS)
			{
				memcpy(l.sNormal,l13.cNormal,sizeof(GzCoord));
				memcpy(l.cNormal,l13.cNormal,sizeof(GzCoord));
				memcpy(l.eNormal,l12.cNormal,sizeof(GzCoord));
			}

			l.start.x = l.current.x = l13.current.x ;
			l.end.x = l12.current.x ;
			l.start.z = l.current.z = l13.current.z ;
			l.end.z = l12.current.z ;


			//Texture
			l.start.u = l.current.u = l13.current.u;
			l.end.u = l12.current.u;
			l.start.v = l.current.v = l13.current.v;
			l.end.v = l12.current.v;

		}
				
		if((l.start.x - l.end.x)!= 0)
		{
			
			if(render->interp_mode == GZ_NORMALS)
			{
				l.mnx[0] = (l.eNormal[0] - l.sNormal[0])/(l.end.x - l.start.x );
				l.mnx[1] = (l.eNormal[1] - l.sNormal[1])/(l.end.x - l.start.x ) ;
				l.mnx[2] = (l.eNormal[2] - l.sNormal[2])/(l.end.x - l.start.x ) ;
			}

			////TEXTURE
			l.mux = (l.start.u - l.end.u) /(l.start.x - l.end.x);
			l.mvx = (l.start.v - l.end.v) /(l.start.x - l.end.x);

			l.mzx = (l.start.z - l.end.z)/ (l.start.x - l.end.x) ;
				
		}
		else
		{
			l.mzx = 0;
			l.mnx[0] = l.mnx[1] = l.mnx[2] = 0;
			l.mux = l.mvx = 0;
		}


		delx = ceil(l.start.x) - l.start.x ;
		l.current.z = l.current.z + l.mzx * delx;
		l.current.x = l.current.x + delx;

		//PHONG
		if(render->interp_mode == GZ_NORMALS)
			interpolateNormalAlongX(&l,delx);
	
		//TEXTURE
		interpolateUVAlongX(&l,delx);

		for(j=ceil(l.start.x); j<(l.end.x) ;j++)
		{
			GzColor C,textColor;
			float u,v;
			
			ReversePerspectiveInterpolation(render,&l,&u,&v);

			if(render->tex_fun!=NULL)
			{
				render->tex_fun(u,v,textColor);
				memcpy(render->Ka,textColor,sizeof(GzColor));
				memcpy(render->Kd,textColor,sizeof(GzColor));
			}

			ComputeColor(render, &C, &l.cNormal);
			memcpy(render->flatcolor,C,sizeof(GzColor));
				
			GzGetDisplay(render->display,j,(int)l13.current.y,&r,&g,&b,&a,&z);
			if((l.start.x == l.current.x) && (l12.drawline == true))
			{
				j++;
				continue;
			}
			else if(l.current.z <=z)
			{
				z = ceil(l.current.z);
				GzPutDisplay(render->display,j,(int)l13.current.y,ctoi(render->flatcolor[0]),ctoi(render->flatcolor[1]),ctoi(render->flatcolor[2]),1,z);
			}
			delx=1;
			l.current.z = l.current.z + l.mzx * delx;
			l.current.x = l.current.x + delx;
				
			//PHONG
			if(render->interp_mode == GZ_NORMALS)
				interpolateNormalAlongX(&l,delx);

			//TEXTURE
			interpolateUVAlongX(&l,delx);

		}

		dely = 1 ;
		l13.current.x = l13.current.x + l13.mxy*dely;
		l13.current.y = l13.current.y +dely;
		l13.current.z = l13.current.z +l13.mzy*dely;

		l12.current.x = l12.current.x + l12.mxy*dely;
		l12.current.y = l12.current.y +dely;
		l12.current.z = l12.current.z +l12.mzy*dely;

		if(render->interp_mode == GZ_NORMALS)
		{
			//Phong
			interpolateNormalAlongY(&l12,dely);
			interpolateNormalAlongY(&l13,dely);
		}

		//TEXTURE
		interpolateUVAlongY(&l13,dely);
		interpolateUVAlongY(&l12,dely);
	}	//END OF WHILE
		
	//i= ceil(l23.start.y);
	dely = l13.current.y - l23.current.y;
	l23.current.x = l23.current.x + l23.mxy*dely;
	l23.current.y = l23.current.y +dely;
	l23.current.z = l23.current.z +l23.mzy*dely;

	//Phong
	if(render->interp_mode == GZ_NORMALS)
		interpolateNormalAlongY(&l23,dely);

	//TEXTURE
	interpolateUVAlongY(&l23,dely);
	

	//for(i; i<ceil(l13.end.y) ; i++)
	while(l13.current.y < l23.end.y)
	{
	
		if(l23.drawline == true)
		{
			if(render->interp_mode == GZ_NORMALS)
			{
				memcpy(l.sNormal,l23.cNormal,sizeof(GzCoord));
				memcpy(l.cNormal,l23.cNormal,sizeof(GzCoord));
				memcpy(l.eNormal,l13.cNormal,sizeof(GzCoord));
			}


			//Texture
			l.start.u = l.current.u = l23.current.u;
			l.end.u = l13.current.u;
			l.start.v = l.current.v = l23.current.v;
			l.end.v = l13.current.v;	


			l.start.x = l.current.x = l23.current.x ;
			l.end.x = l13.current.x ;
			l.start.z = l.current.z = l23.current.z ;
			l.end.z = l13.current.z ;
			
		}
		else
		{
			if(render->interp_mode == GZ_NORMALS)
			{
				memcpy(l.sNormal,l13.cNormal,sizeof(GzCoord));
				memcpy(l.cNormal,l13.cNormal,sizeof(GzCoord));
				memcpy(l.eNormal,l23.cNormal,sizeof(GzCoord));
			}


			//Texture
			l.start.u = l.current.u = l13.current.u;
			l.end.u = l23.current.u;
			l.start.v = l.current.v = l13.current.v;
			l.end.v = l23.current.v;

			l.start.x = l.current.x = l13.current.x ;
			l.end.x = l23.current.x ;
			l.start.z = l.current.z = l13.current.z ;
			l.end.z = l23.current.z ;
		}

		
		if((l.start.x - l.end.x)!= 0)
		{
			if(render->interp_mode == GZ_NORMALS)
			{
				l.mnx[0] = (l.eNormal[0] - l.sNormal[0])/(l.end.x - l.start.x ) ;
				l.mnx[1] = (l.eNormal[1] - l.sNormal[1])/(l.end.x - l.start.x ) ;
				l.mnx[2] = (l.eNormal[2] - l.sNormal[2])/(l.end.x - l.start.x ) ;
			}
			
			
			//TEXTURE
			l.mux = (l.start.u - l.end.u) /(l.start.x - l.end.x);
			l.mvx = (l.start.v - l.end.v) /(l.start.x - l.end.x);
			
			l.mzx = (l.start.z - l.end.z)/ (l.start.x - l.end.x) ;
		}
		else
		{
			l.mzx = 0;
			l.mnx[0] = l.mnx[1] = l.mnx[2] = 0;
			
		}

		delx = ceil(l.start.x) - l.start.x ;
		l.current.x = l.current.x + delx;
		l.current.z = l.current.z + l.mzx * delx;

		if(render->interp_mode == GZ_NORMALS)
			interpolateNormalAlongX(&l,delx);

		//Texture
		interpolateUVAlongX(&l,delx);

		for(j=ceil(l.start.x); j<(l.end.x) ;j++)
		{
			GzColor C,textColor;
			float u,v;
			
			ReversePerspectiveInterpolation(render,&l,&u,&v);

			if(render->tex_fun!=NULL)
			{
				render->tex_fun(u,v,textColor);
				memcpy(render->Ka,textColor,sizeof(GzColor));
				memcpy(render->Kd,textColor,sizeof(GzColor));
			}
			ComputeColor(render, &C, &l.cNormal);
			memcpy(render->flatcolor,C,sizeof(GzColor));
			
			GzGetDisplay(render->display,j,(int)l13.current.y,&r,&g,&b,&a,&z);
			if(l.start.x == l.current.x && (l12.drawline == true))
			{
				j++;
				continue;
			}
			else if(l.current.z <=z)
			{
				z = (GzDepth)ceil(l.current.z);
				GzPutDisplay(render->display,j,(int)l13.current.y,ctoi(render->flatcolor[0]),ctoi(render->flatcolor[1]),ctoi(render->flatcolor[2]),1,z);
			}
			delx=1;
			l.current.x = l.current.x + delx;
			l.current.z = l.current.z + l.mzx * delx;
			
			if(render->interp_mode == GZ_NORMALS)
				interpolateNormalAlongX(&l,delx);

			//Texture
			interpolateUVAlongX(&l,delx);
		}

		dely = 1 ;
		l13.current.x = l13.current.x + l13.mxy*dely;
		l13.current.y = l13.current.y + dely ;
		l13.current.z = l13.current.z +l13.mzy*dely;

		l23.current.x = l23.current.x + l23.mxy*dely;
		l23.current.y = l23.current.y +dely;
		l23.current.z = l23.current.z +l23.mzy*dely;

		if(render->interp_mode == GZ_NORMALS)
		{
			interpolateNormalAlongY(&l13,dely);
			interpolateNormalAlongY(&l23,dely);
		}

		//TEXTURE
		interpolateUVAlongY(&l13,dely);
		interpolateUVAlongY(&l23,dely);

	}
}


void Rasterize(GzRender*render,vertice a1, vertice a2, vertice a3, dda l12, dda l13,dda l23)
{
	int dir=0, j,flag=0;
	float dely , delx;
	GzIntensity r,g,b,a;
	GzDepth z;
	bool horzLine = false;
	dda l;
	GzColor currentColor;

	if(a1.y == a2.y || a2.y ==  a3.y )
			horzLine = true;

	if(l12.mxy < l13.mxy)
	{
		l12.drawline = true;
		l23.drawline = true;
		l13.drawline = false;
	}
	else
	{
		l12.drawline = false;
		l23.drawline = false;
		l13.drawline = true;
	}

		//i = ceil(l13.current.y);
	dely = ceil(l13.current.y) - l13.current.y;
	l13.current.x = l13.current.x + l13.mxy*dely;
	l13.current.y = l13.current.y +dely;
	l13.current.z = l13.current.z +l13.mzy*dely;

	l12.current.x = l12.current.x + l12.mxy*dely;
	l12.current.y = l12.current.y +dely;
	l12.current.z = l12.current.z +l12.mzy*dely;

	//PHONG
	if(render->interp_mode == GZ_NORMALS)
	{
		interpolateNormalAlongY(&l13,dely);
		interpolateNormalAlongY(&l12,dely);
	}
	
	//Gouraud
	if(render->interp_mode == GZ_COLOR)
	{
		interpolateColorAlongY(&l13,dely);
		interpolateColorAlongY(&l12,dely);
	}

	//TEXTURE INTERPOLATION
	
	interpolateUVAlongY(&l13,dely);
	interpolateUVAlongY(&l12,dely);

	while(l13.current.y < l12.end.y)
	//for(i; i<(l12.end.y) ; i++)
	{
		if(l12.drawline == true)
		{
			if(render->interp_mode == GZ_NORMALS)
			{
				memcpy(l.sNormal,l12.cNormal,sizeof(GzCoord));
				memcpy(l.cNormal,l12.cNormal,sizeof(GzCoord));
				memcpy(l.eNormal,l13.cNormal,sizeof(GzCoord));
			}

			if(render->interp_mode == GZ_COLOR)
			{
				memcpy(l.sColor,l12.cColor,sizeof(GzColor));
				memcpy(l.cColor,l12.cColor,sizeof(GzColor));
				memcpy(l.eColor,l13.cColor,sizeof(GzColor));
			}

			l.start.x = l.current.x = l12.current.x ;
			l.end.x = l13.current.x ;
			l.start.z = l.current.z = l12.current.z ;
			l.end.z = l13.current.z ;

			//Texture
			l.start.u = l.current.u = l12.current.u;
			l.end.u = l13.current.u;
			l.start.v = l.current.v = l12.current.v;
			l.end.v = l13.current.v;

		}
		else
		{
			if(render->interp_mode == GZ_NORMALS)
			{
				memcpy(l.sNormal,l13.cNormal,sizeof(GzCoord));
				memcpy(l.cNormal,l13.cNormal,sizeof(GzCoord));
				memcpy(l.eNormal,l12.cNormal,sizeof(GzCoord));
			}
			
			if(render->interp_mode == GZ_COLOR)
			{
				memcpy(l.sColor,l13.cColor,sizeof(GzColor));
				memcpy(l.cColor,l13.cColor,sizeof(GzColor));
				memcpy(l.eColor,l12.cColor,sizeof(GzColor));
			}
			
			l.start.x = l.current.x = l13.current.x ;
			l.end.x = l12.current.x ;
			l.start.z = l.current.z = l13.current.z ;
			l.end.z = l12.current.z ;


			//Texture
			l.start.u = l.current.u = l13.current.u;
			l.end.u = l12.current.u;
			l.start.v = l.current.v = l13.current.v;
			l.end.v = l12.current.v;

		}

		if((l.start.x - l.end.x)!= 0)
		{
			
			if(render->interp_mode == GZ_NORMALS)
			{
				l.mnx[0] = (l.eNormal[0] - l.sNormal[0])/(l.end.x - l.start.x );
				l.mnx[1] = (l.eNormal[1] - l.sNormal[1])/(l.end.x - l.start.x ) ;
				l.mnx[2] = (l.eNormal[2] - l.sNormal[2])/(l.end.x - l.start.x ) ;
			}

			if(render->interp_mode == GZ_COLOR)
			{
				l.mcx[0] = (l.eColor[0] - l.sColor[0])/(l.end.x - l.start.x ) ;
				l.mcx[1] = (l.eColor[1] - l.sColor[1])/(l.end.x - l.start.x ) ;
				l.mcx[2] = (l.eColor[2] - l.sColor[2])/(l.end.x - l.start.x ) ;
			}

			
			//TEXTURE
			l.mux = (l.start.u - l.end.u) /(l.start.x - l.end.x);
			l.mvx = (l.start.v - l.end.v) /(l.start.x - l.end.x);

			l.mzx = (l.start.z - l.end.z)/ (l.start.x - l.end.x) ;
				
		}
		else
		{
			l.mzx = 0;
			l.mnx[0] = l.mnx[1] = l.mnx[2] = 0;
			l.mux = l.mvx = 0;
		}


		delx = ceil(l.start.x) - l.start.x ;
		l.current.z = l.current.z + l.mzx * delx;
		l.current.x = l.current.x + delx;

		//PHONG
		if(render->interp_mode == GZ_NORMALS)
			interpolateNormalAlongX(&l,delx);

		//Gouraud
		if(render->interp_mode == GZ_COLOR)
			interpolateColorAlongX(&l,delx);	
			
		//TEXTURE
		interpolateUVAlongX(&l,delx);

		for(j=ceil(l.start.x); j<(l.end.x) ;j++)
		{
			GzColor C,textColor;
			
			render->tex_fun(l.current.u,l.current.v,textColor);
			
			if(render->interp_mode == GZ_COLOR)
			{
				memcpy(render->Ka,textColor,sizeof(GzColor));
				memcpy(render->Kd,textColor,sizeof(GzColor));
				memcpy(render->Ks,textColor,sizeof(GzColor));
			}
			if(render->interp_mode == GZ_NORMALS)
			{
				memcpy(render->Ka,textColor,sizeof(GzColor));
				memcpy(render->Kd,textColor,sizeof(GzColor));
			}

			ComputeColor(render, &C, &l.cNormal);
			memcpy(render->flatcolor,C,sizeof(GzColor));
				
			GzGetDisplay(render->display,j,(int)l13.current.y,&r,&g,&b,&a,&z);
			if((l.start.x == l.current.x) && (l12.drawline == true))
			{
				j++;
				continue;
			}
			else if(l.current.z <=z)
			{
				z = ceil(l.current.z);
				GzPutDisplay(render->display,j,(int)l13.current.y,ctoi(render->flatcolor[0]),ctoi(render->flatcolor[1]),ctoi(render->flatcolor[2]),1,z);
			}
			delx=1;
			l.current.z = l.current.z + l.mzx * delx;
			l.current.x = l.current.x + delx;
				
			//PHONG
			if(render->interp_mode == GZ_NORMALS)
				interpolateNormalAlongX(&l,delx);

			//Gouraud
			if(render->interp_mode == GZ_COLOR)
				interpolateColorAlongX(&l,delx);	
				
			//TEXTURE
			interpolateUVAlongX(&l,delx);

		}

		dely = 1 ;
		l13.current.x = l13.current.x + l13.mxy*dely;
		l13.current.y = l13.current.y +dely;
		l13.current.z = l13.current.z +l13.mzy*dely;

		l12.current.x = l12.current.x + l12.mxy*dely;
		l12.current.y = l12.current.y +dely;
		l12.current.z = l12.current.z +l12.mzy*dely;

		if(render->interp_mode == GZ_NORMALS)
		{
			//Phong
			interpolateNormalAlongY(&l12,dely);
			interpolateNormalAlongY(&l13,dely);
		}

		if(render->interp_mode == GZ_COLOR)
		{
			//Gouraud
			interpolateColorAlongY(&l12,dely);
			interpolateColorAlongY(&l13,dely);
		}
		
		//TEXTURE
		interpolateUVAlongY(&l13,dely);
		interpolateUVAlongY(&l12,dely);
	}	//END OF WHILE
		
	//i= ceil(l23.start.y);
	dely = l13.current.y - l23.current.y;
	l23.current.x = l23.current.x + l23.mxy*dely;
	l23.current.y = l23.current.y +dely;
	l23.current.z = l23.current.z +l23.mzy*dely;

	//Phong
	if(render->interp_mode == GZ_NORMALS)
		interpolateNormalAlongY(&l23,dely);

	//Gouraud
	if(render->interp_mode == GZ_COLOR)
		interpolateColorAlongY(&l23,dely);	
		
	//TEXTURE
	interpolateUVAlongY(&l23,dely);
	

	//for(i; i<ceil(l13.end.y) ; i++)
	while(l13.current.y < l23.end.y)
	{
	
		if(l23.drawline == true)
		{
			if(render->interp_mode == GZ_NORMALS)
			{
				memcpy(l.sNormal,l23.cNormal,sizeof(GzCoord));
				memcpy(l.cNormal,l23.cNormal,sizeof(GzCoord));
				memcpy(l.eNormal,l13.cNormal,sizeof(GzCoord));
			}

			
			if(render->interp_mode == GZ_COLOR)
			{
				memcpy(l.sColor,l23.cColor,sizeof(GzColor));
				memcpy(l.cColor,l23.cColor,sizeof(GzColor));
				memcpy(l.eColor,l13.cColor,sizeof(GzColor));
			}

			//Texture
			l.start.u = l.current.u = l23.current.u;
			l.end.u = l13.current.u;
			l.start.v = l.current.v = l23.current.v;
			l.end.v = l13.current.v;	


			l.start.x = l.current.x = l23.current.x ;
			l.end.x = l13.current.x ;
			l.start.z = l.current.z = l23.current.z ;
			l.end.z = l13.current.z ;
			
		}
		else
		{
			if(render->interp_mode == GZ_NORMALS)
			{
				memcpy(l.sNormal,l13.cNormal,sizeof(GzCoord));
				memcpy(l.cNormal,l13.cNormal,sizeof(GzCoord));
				memcpy(l.eNormal,l23.cNormal,sizeof(GzCoord));
			}

			if(render->interp_mode == GZ_COLOR)
			{
				memcpy(l.sColor,l13.cColor,sizeof(GzColor));
				memcpy(l.cColor,l13.cColor,sizeof(GzColor));
				memcpy(l.eColor,l23.cColor,sizeof(GzColor));
			}
			
			//Texture
			l.start.u = l.current.u = l13.current.u;
			l.end.u = l23.current.u;
			l.start.v = l.current.v = l13.current.v;
			l.end.v = l23.current.v;

			l.start.x = l.current.x = l13.current.x ;
			l.end.x = l23.current.x ;
			l.start.z = l.current.z = l13.current.z ;
			l.end.z = l23.current.z ;
		}

		if((l.start.x - l.end.x)!= 0)
		{
			if(render->interp_mode == GZ_NORMALS)
			{
				l.mnx[0] = (l.eNormal[0] - l.sNormal[0])/(l.end.x - l.start.x ) ;
				l.mnx[1] = (l.eNormal[1] - l.sNormal[1])/(l.end.x - l.start.x ) ;
				l.mnx[2] = (l.eNormal[2] - l.sNormal[2])/(l.end.x - l.start.x ) ;
			}
			
			if(render->interp_mode == GZ_COLOR)
			{
				l.mcx[0] = (l.eColor[0] - l.sColor[0])/(l.end.x - l.start.x ) ; ;
				l.mcx[1] = (l.eColor[1] - l.sColor[1])/(l.end.x - l.start.x ) ; ;
				l.mcx[2] = (l.eColor[2] - l.sColor[2])/(l.end.x - l.start.x ) ; ;
			}
			
			
			//TEXTURE
			l.mux = (l.start.u - l.end.u) /(l.start.x - l.end.x);
			l.mvx = (l.start.v - l.end.v) /(l.start.x - l.end.x);
			
			l.mzx = (l.start.z - l.end.z)/ (l.start.x - l.end.x) ;
		}
		else
		{
			l.mzx = 0;
			l.mnx[0] = l.mnx[1] = l.mnx[2] = 0;
			l.mux = l.mvx = 0;
		}

		delx = ceil(l.start.x) - l.start.x ;
		l.current.x = l.current.x + delx;
		l.current.z = l.current.z + l.mzx * delx;

		if(render->interp_mode == GZ_NORMALS)
			interpolateNormalAlongX(&l,delx);

		if(render->interp_mode == GZ_COLOR)
			interpolateColorAlongX(&l,delx);
		
		//Texture
		interpolateUVAlongX(&l,delx);

		for(j=ceil(l.start.x); j<(l.end.x) ;j++)
		{
			GzColor C,textColor;
			
			render->tex_fun(l.current.u,l.current.v,textColor);
			if(render->interp_mode == GZ_COLOR)
			{
				memcpy(render->Ka,textColor,sizeof(GzColor));
				memcpy(render->Kd,textColor,sizeof(GzColor));
				memcpy(render->Ks,textColor,sizeof(GzColor));
			}
			if(render->interp_mode == GZ_NORMALS)
			{
				memcpy(render->Ka,textColor,sizeof(GzColor));
				memcpy(render->Kd,textColor,sizeof(GzColor));
			}

			ComputeColor(render, &C, &l.cNormal);
			memcpy(render->flatcolor,C,sizeof(GzColor));
			
			GzGetDisplay(render->display,j,(int)l13.current.y,&r,&g,&b,&a,&z);
			if(l.start.x == l.current.x && (l12.drawline == true))
			{
				j++;
				continue;
			}
			else if(l.current.z <=z)
			{
				z = (GzDepth)ceil(l.current.z);
				GzPutDisplay(render->display,j,(int)l13.current.y,ctoi(render->flatcolor[0]),ctoi(render->flatcolor[1]),ctoi(render->flatcolor[2]),1,z);
			}
			delx=1;
			l.current.x = l.current.x + delx;
			l.current.z = l.current.z + l.mzx * delx;
			
			if(render->interp_mode == GZ_NORMALS)
				interpolateNormalAlongX(&l,delx);

			if(render->interp_mode == GZ_COLOR)
				interpolateColorAlongX(&l,delx);
			
			//Texture
			interpolateUVAlongX(&l,delx);
		}

		dely = 1 ;
		l13.current.x = l13.current.x + l13.mxy*dely;
		l13.current.y = l13.current.y + dely ;
		l13.current.z = l13.current.z +l13.mzy*dely;

		l23.current.x = l23.current.x + l23.mxy*dely;
		l23.current.y = l23.current.y +dely;
		l23.current.z = l23.current.z +l23.mzy*dely;

		if(render->interp_mode == GZ_NORMALS)
		{
			interpolateNormalAlongY(&l13,dely);
			interpolateNormalAlongY(&l23,dely);
		}
		
		if(render->interp_mode == GZ_COLOR)
		{
			interpolateColorAlongY(&l13,dely);
			interpolateColorAlongY(&l23,dely);
		}

		//TEXTURE
		interpolateUVAlongY(&l13,dely);
		interpolateUVAlongY(&l23,dely);

	}
}


void AAShift(GzRender * render, vertice * a)
{
	a->x = a->x + render->aashiftX;
	a->y = a->y + render->aashiftY;

}


int GzPutTriangle(GzRender	*render, int numParts, GzToken *nameList, 
				  GzPointer	*valueList)
/* numParts : how many names and values */
{
/*  
- pass in a triangle description with tokens and values corresponding to 
      GZ_POSITION:3 vert positions in model space 
- Xform positions of verts  
- Clip - just discard any triangle with verts behind view plane 
       - test for triangles with all three verts off-screen 
- invoke triangle rasterizer  
*/ 
	int dir=0, i, j,index,flag=0;
	float dely , delx;
	dda l12,l23,l13,l;
	GzIntensity r,g,b,a;
	GzDepth z;
	//3 vertices
	vertice a1,a2,a3;
	//3 dummy vertices used for multiplication
	vertice v1,v2,v3;
	
	for(index=0;index<numParts;index++)
	{
		switch(nameList[index])
		{		
		case GZ_POSITION:
			a1 = StoreVertex(((GzCoord*)valueList[0])[0][0],((GzCoord*)valueList[0])[0][1],((GzCoord*)valueList[0])[0][2]);
			a2 = StoreVertex(((GzCoord*)valueList[0])[1][0],((GzCoord*)valueList[0])[1][1],((GzCoord*)valueList[0])[1][2]);
			a3 = StoreVertex(((GzCoord*)valueList[0])[2][0],((GzCoord*)valueList[0])[2][1],((GzCoord*)valueList[0])[2][2]);
			break;
		
		case GZ_NORMAL:
			StoreNormal(&a1, ((GzCoord*)valueList[1])[0][0],((GzCoord*)valueList[1])[0][1],((GzCoord*)valueList[1])[0][2]);
			StoreNormal(&a2, ((GzCoord*)valueList[1])[1][0],((GzCoord*)valueList[1])[1][1],((GzCoord*)valueList[1])[1][2]);
			StoreNormal(&a3, ((GzCoord*)valueList[1])[2][0],((GzCoord*)valueList[1])[2][1],((GzCoord*)valueList[1])[2][2]);
			break;	

		case GZ_TEXTURE_INDEX:
			StoreUVCoord(&a1, ((GzTextureIndex*)valueList[2])[0][0],((GzTextureIndex*)valueList[2])[0][1]);
			StoreUVCoord(&a2, ((GzTextureIndex*)valueList[2])[1][0],((GzTextureIndex*)valueList[2])[1][1]);
			StoreUVCoord(&a3, ((GzTextureIndex*)valueList[2])[2][0],((GzTextureIndex*)valueList[2])[2][1]);
			break;	
		}

	}

	MultiplyMatrixVertex(render->Ximage[render->matlevel], &a1);
	MultiplyMatrixVertex(render->Ximage[render->matlevel], &a2);
	MultiplyMatrixVertex(render->Ximage[render->matlevel], &a3);

	MultiplyMatrixNormal(render->Xnorm[render->matlevel], &a1);
	MultiplyMatrixNormal(render->Xnorm[render->matlevel], &a2);
	MultiplyMatrixNormal(render->Xnorm[render->matlevel], &a3);
	
	//Skip the triangle if z value of any vertice is less than 0
	if(a1.z <0 || a2.z <0 || a3.z<0)
		return GZ_SUCCESS;

	SortVertice(&a1,&a2,&a3);
	
	AAShift(render,&a1);
	AAShift(render,&a2);
	AAShift(render,&a3);

	filldda(&l12,a1,a2);
	filldda(&l13,a1,a3);
	filldda(&l23,a2,a3);

	PerspectiveInterpolation(render,&a1);
	PerspectiveInterpolation(render,&a2);
	PerspectiveInterpolation(render,&a3);

	addUVCoord(&l12,a1,a2);
	addUVCoord(&l23,a2,a3);
	addUVCoord(&l13,a1,a3);

	//CASE FLAT SHADING
	GzCoord normal;
	if(render->interp_mode == GZ_NONE)
	{
		GzColor C1;
		CreateNormalVector(a1, (float*)&normal);
		ComputeColor(render, &C1, &normal);
		render->flatcolor[0] = C1[0];
		render->flatcolor[1] = C1[1];
		render->flatcolor[2] = C1[2];
	}
	//CASE GOURAUD SHADING
	else if(render->interp_mode == GZ_COLOR)
	{
		GzColor identity = {1,1,1};

			//calculate color at vertices
		memcpy(render->Ka,identity,sizeof(GzColor));
		memcpy(render->Ks,identity,sizeof(GzColor));
		memcpy(render->Kd,identity,sizeof(GzColor));

		GzColor c1,c2,c3;
		CreateNormalVector(a1, (float*)&normal);
		ComputeColor(render, &c1, &normal);
		CreateNormalVector(a2, (float*)&normal);
		ComputeColor(render, &c2, &normal);
		CreateNormalVector(a3, (float*)&normal);
		ComputeColor(render, &c3, &normal);

		fillddaColor(&l12, c1,c2);
		fillddaColor(&l13, c1,c3);
		fillddaColor(&l23, c2,c3);

		GouraudRasterize(render,a1,a2,a3,l12,l13,l23);
		//Rasterize(render,a1,a2,a3,l12,l13,l23);
		
	}
	//PHONG SHADING
	else if(render->interp_mode == GZ_NORMAL)
	{
		//calculate color at vertices
		GzCoord Normal1, Normal2, Normal3;
		CreateNormalVector(a1, (float*)&Normal1);
		CreateNormalVector(a2, (float*)&Normal2);
		CreateNormalVector(a3, (float*)&Normal3);
		
		fillddaNormal(&l12, Normal1,Normal2);
		fillddaNormal(&l13, Normal1,Normal3);
		fillddaNormal(&l23, Normal2,Normal3);
		
		PhongRasterize(render,a1,a2,a3,l12,l13,l23);
		//Rasterize(render,a1,a2,a3,l12,l13,l23);

	}

	return GZ_SUCCESS;
}


