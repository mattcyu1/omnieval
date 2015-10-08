#include "map_utils.h"

// ==============================================================================
//                               GLOBAL VAR DEFS
// ==============================================================================
point G_CENT_POINTS[] = {
    { 0.5f, 0.5f },
};
float* G_INP_LUT = NULL;
float* G_OUT_LUT = NULL;
int    G_NUMPTS_LUT = 100000;
float  G_iK[3][3];
float  G_R[3][3];
bool   G_ACSFLAG = false;
float* G_LAT_INP_LUT = NULL;
float* G_LAT_OUT_LUT = NULL;
int    G_LAT_NUMPTS_LUT = 400;
float  G_MAP_OFF_Y = 0.;
float  G_MAP_OFF_X = 0.;
float  G_MAP_SF_Y  = 0.;
float  G_MAP_SF_X  = 0.;
const float G_OV  = 0.02;
const float G_OVX = G_OV/3.;

// ==============================================================================
//                               MAPPING LOOKUP TABLES
// ==============================================================================
void setupLatLUT(const char* fname){
    FILE *fp = fopen(fname,"r");
    G_LAT_INP_LUT = (float*) malloc (sizeof(float)*G_LAT_NUMPTS_LUT);
    G_LAT_OUT_LUT = (float*) malloc (sizeof(float)*G_LAT_NUMPTS_LUT);
    for(int z=0; z<G_LAT_NUMPTS_LUT; z++)
	fscanf(fp, "%f %f", &G_LAT_INP_LUT[z], &G_LAT_OUT_LUT[z]);
    fclose(fp);
}

float latLUTLookup(float x){
    x = x*PI/180.0f;
    float step = PI/(G_LAT_NUMPTS_LUT-1);
    float idx = (x+PIBY2)/step;
    int lowidx = (int)floor(idx);
    int higidx = (int)ceil(idx);
    if(lowidx==higidx)
	return G_LAT_OUT_LUT[lowidx];

    float d = idx-lowidx;
    return lerp(G_LAT_OUT_LUT[lowidx],G_LAT_OUT_LUT[higidx],d);
}


void setupCos2(){
    G_INP_LUT = (float*) malloc (sizeof(float)*G_NUMPTS_LUT);
    G_OUT_LUT = (float*) malloc (sizeof(float)*G_NUMPTS_LUT);

    float lut_step = PI/G_NUMPTS_LUT;
    for (int i=0; i<G_NUMPTS_LUT; i++){
	float x = -PIBY2+lut_step*i;
	G_INP_LUT[i] = x;
	G_OUT_LUT[i] = 2/PI*(x+sin(x)*cos(x));
    }
}

/*
 * bSearch_LUT
 * --------------------------------------------------
 * Descrp: Given that a lookup table is monotically increasing,
 *         we use binary search to find the idx corresponding
 *         to the value closest to p
 * Inputs: (float) target val
 *         (int)   current lower bound
 *         (int)   current upper bound
 * Return: (int)   idx
 */
int bSearch_LUT(float p, int minVal, int maxVal){
    int idx = (maxVal+minVal)/2;
    
    if(G_OUT_LUT[idx]==p)
	return idx;
    else if (minVal==maxVal)
	return idx;
    else if (minVal==idx)
	return idx;
    else if(G_OUT_LUT[idx] < p)
	return bSearch_LUT(p,idx,maxVal);
    else
	return bSearch_LUT(p,minVal,idx);
}

float LUTinv(float p){
    // find closest idx
    int idx0 = bSearch_LUT(p,0,G_NUMPTS_LUT);

    // then bilinear interpolate
    int idx1;
    float x0 = G_OUT_LUT[idx0];
    if (x0==p)
	return G_INP_LUT[idx0];
    else if (x0 > p)
	idx1 = idx0-1;
    else
	idx1 = idx0+1;
    float x1 = G_OUT_LUT[idx1];

    float diff = fabs(x1-x0);
    float diff0 = fabs(x0-p);
    float alph = 1.0-diff0/diff;
    return alph*G_INP_LUT[idx0]+(1-alph)*G_INP_LUT[idx1];
}

// ==============================================================================
//                               MAPPING TILING HELPER
// ==============================================================================
void setMapOffset(float2 t){
    G_MAP_OFF_X = t.x;
    G_MAP_OFF_Y = t.y;
}

void setMapScalefactor(float2 t){
    G_MAP_SF_X = t.x;
    G_MAP_SF_Y = t.y;
}

// ==============================================================================
//                               MAPPING VIEWPORT HELPER
// ==============================================================================
void clearTrData(){
    G_R[0][0] =  1.0f;    G_R[0][1] =  0.0f;    G_R[0][2] =  0.0f;
    G_R[1][0] =  0.0f;    G_R[1][1] =  1.0f;    G_R[1][2] =  0.0f;
    G_R[2][0] =  0.0f;    G_R[2][1] =  0.0f;    G_R[2][2] =  1.0f;
}

void setRotationMat(float p, float t){
    float phi = p*M_PI/180.0;
    float tht = -t*M_PI/180.0;
    
    G_R[0][0] =  cosf(tht);    G_R[0][1] = sinf(tht)*sinf(phi);    G_R[0][2] = sinf(tht)*cosf(phi);
    G_R[1][0] =       0.0f;    G_R[1][1] =           cosf(phi);    G_R[1][2] =          -sinf(phi);
    G_R[2][0] = -sinf(tht);    G_R[2][1] = cosf(tht)*sinf(phi);    G_R[2][2] = cosf(tht)*cosf(phi);
}

void setRotationMat(float *rt, float *up, float *lk){
    normalize(rt);
    normalize(up);
    normalize(lk);

    G_R[0][0] = rt[0];    G_R[0][1] = rt[1];    G_R[0][2] = rt[2];
    G_R[1][0] = up[0];    G_R[1][1] = up[1];    G_R[1][2] = up[2];
    G_R[2][0] = lk[0];    G_R[2][1] = lk[1];    G_R[2][2] = lk[2];
}

/*
 * setIntrMat
 * ------------------------------------
 * Descrp: Camera K matrix
 * Inputs: Viewport width w, height h
 * Outputs: Computes K and inv(K)
 */
void setIntrMat(int dw, int dh, int sw, int sh, float fvx, float fvy){
    float fovx = TWOPI * fvx/360.0;
    float fovy = TWOPI * fvy/360.0;

    float fx = dw/2 * 1/tanf(fovx/2);
    float fy = dh/2 * 1/tanf(fovy/2);

    float K[3][3]={{fx,0,dw/2.0f},{0,-fy,dh/2.0f},{0,0,1}};
    invert3x3(K);
}

/*
 * invert3x3
 * ------------------------------------
 * Descrp: computes the inverse of a 3x3 matrix
 * Inputs: K
 * Outputs: set global variable G_iK=inv(K)
 */
int invert3x3(float K[3][3]){
    float det = K[0][0] * (K[1][1] * K[2][2] - K[2][1] * K[1][2]) -
	        K[0][1] * (K[1][0] * K[2][2] - K[1][2] * K[2][0]) +
	        K[0][2] * (K[1][0] * K[2][1] - K[1][1] * K[2][0]);

    if ( fabs(det) < 1e-2 ){
        memset( G_iK, 0, 9*sizeof(float) );
        return 1;
    }

    float idet = 1 / det;
    G_iK[0][0] = (K[1][1] * K[2][2] - K[2][1] * K[1][2]) * idet;
    G_iK[0][1] = (K[0][2] * K[2][1] - K[0][1] * K[2][2]) * idet;
    G_iK[0][2] = (K[0][1] * K[1][2] - K[0][2] * K[1][1]) * idet;
    G_iK[1][0] = (K[1][2] * K[2][0] - K[1][0] * K[2][2]) * idet;
    G_iK[1][1] = (K[0][0] * K[2][2] - K[0][2] * K[2][0]) * idet;
    G_iK[1][2] = (K[1][0] * K[0][2] - K[0][0] * K[1][2]) * idet;
    G_iK[2][0] = (K[1][0] * K[2][1] - K[2][0] * K[1][1]) * idet;
    G_iK[2][1] = (K[2][0] * K[0][1] - K[0][0] * K[2][1]) * idet;
    G_iK[2][2] = (K[0][0] * K[1][1] - K[1][0] * K[0][1]) * idet;

    return 0;
}

int readNextTrData(FILE* ft, int nf, bool& bTrdata){
    int tr = 0;     // frame num
    float ts = 0;   // time stamp

    if ( feof(ft) || !fscanf(ft, "%d,", &tr) ){
        bTrdata = false;
        printf("Frame: %d\n", nf);
        printf ("Track file shorter than yuv..using last track\n");
        return 1;
    }

    // read frame num
    if(tr != nf){
        printf("Frame: %d\n", nf);
        printf ("Track frame num does not match yuv frame num..using last track\n");
        return 1;
    }

    // read time stamp
    fscanf(ft, "%f,", &ts);

    // read vectors
    float rt[3], up[3], lk[3];
    fscanf(ft, "%f,", &rt[0]); fscanf(ft, "%f,", &rt[1]); fscanf(ft, "%f,", &rt[2]);
    fscanf(ft, "%f,", &up[0]); fscanf(ft, "%f,", &up[1]); fscanf(ft, "%f,", &up[2]);
    fscanf(ft, "%f,", &lk[0]); fscanf(ft, "%f,", &lk[1]); fscanf(ft, "%f,", &lk[2]);

    // rotation matrix
    setRotationMat(rt, up, lk);

    return 0;
}

/*
 * sphToCart
 * -------------------------------------
 * Descrp: Converts a spherical point to cartesian coordinates assuming radius of sphere is 1
 * Inputs: (float2) sph point
 * Return: (float3) cart point
 */
float3 sphToCart(float2 sph){
    float lat = sph.x*PI/180.0f;
    float lon = sph.y*PI/180.0f;
    
    float3 out;
    out.x =  sinf(lon) * cosf(lat);
    out.y =              sinf(lat);
    out.z = -cosf(lon) * cosf(lat);
    return out;
}

// ==============================================================================
//                               MAPPINGS
// ==============================================================================
// sph2[mapping]
int sph2rect(int* f, float* i, float* j, const image* img, const float* v, int b){
    int w = img->w;
    int h = img->h;

    float x, y, z;
    x = v[0]; y = v[1]; z = v[2];
    
    float phi = acosf(y);
    float theta = atan2f(x, -z);

    *f = 0;
    *i = h * (phi / PI);
    *j = w * (0.5f + theta / TWOPI);

    return 1;
}

int sph2eqar(int *f, float *i, float *j, const image* img, const float *v, int b){
    int w = img->w;
    int h = img->h;

    float x, y, z;
    x = v[0]; y = v[1]; z = v[2];

    float phi = acosf(y) - PI/2;    // range: -pi/2 to pi/2
    float theta = atan2f(x, -z);    // range: -pi to pi

    *f = 0;
    *i = h * ((1 + sinf(phi))/2);
    *j = w * (0.5f + theta / TWOPI);

    return 1;
}

int sph2merc(int *f, float *i, float *j, const image* img, const float *v, int b){
    int w = img->w;
    int h = img->h;

    float x, y, z;
    x = v[0]; y = v[1]; z = v[2];

    float phi = acosf(y) - PI/2;    // range: -pi/2 to pi/2
    float theta = atan2f(x, -z);    // range: -pi to pi

    *f = 0;
    *i = h * (3.0+fmin(log(tan(phi)+1.0/cos(phi)),3.0))/6;
    *j = w * (0.5f + theta / TWOPI);

    return 1;
}

int sph2cos2(int *f, float *i, float *j, const image* img, const float *v, int b){
    int w = img->w;
    int h = img->h;

    float x, y, z;
    x = v[0]; y = v[1]; z = v[2];

    float phi = acosf(y) - PI/2;    // range: -pi/2 to pi/2
    float theta = atan2f(x, -z);    // range: -pi to pi

    *f = 0;
    *i = h * ((1 + 2/PI*(phi+sin(phi)*cos(phi)))/2);
    *j = w * (0.5f + theta / TWOPI);

    return 1;
}

// Note that dyad gets broken because of the lack of ability to change w
// should probably change this to the old verison which can handle modifying w/h
int sph2dyad(int *f, float *i, float *j, const image* img, const float *v, int b){
    int w = img->w;
    int h = img->h;

    float x, y, z;
    x = v[0]; y = v[1]; z = v[2];
    float phi = acosf(y)-PIBY2;

    if(fabs(phi) <= PI/3.0)
	return sph2rect(f,i,j,img,v,b);
    return sph2rect(f,i,j,img,v,b);
}

int sph2mult(int *f, float *i, float *j, const image* img, const float *v, int b){
    int w,h;
    const float X = fabsf(v[0]);
    const float Y = fabsf(v[1]);
    const float Z = fabsf(v[2]);
    float x,y;

    if     (v[2] < 0 && Z >= X && Z >= Y) { *f = 0; x =  v[0] / Z; y = -v[1] / Z; w = img[0].w; h = img[0].h;} // pFront
    else if(v[0] > 0 && X >= Y && X >= Z) { *f = 1; x =  v[2] / X; y = -v[1] / X; w = img[1].w; h = img[1].h;} // pRight
    else if(v[2] > 0 && Z >= X && Z >= Y) { *f = 2; x = -v[0] / Z; y = -v[1] / Z; w = img[2].w; h = img[2].h;} // pBack
    else if(v[0] < 0 && X >= Y && X >= Z) { *f = 3; x = -v[2] / X; y = -v[1] / X; w = img[3].w; h = img[3].h;} // pLeft
    else if(v[1] > 0 && Y >= X && Y >= Z) { *f = 4; x =  v[0] / Y; y = -v[2] / Y; w = img[4].w; h = img[4].h;} // pTop
    else if(v[1] < 0 && Y >= X && Y >= Z) { *f = 5; x =  v[0] / Y; y =  v[2] / Y; w = img[5].w; h = img[5].h;} // pBottom
    else return 0;

    *i = h*(y+1.0f) / 2.0f;
    *j = w*(x+1.0f) / 2.0f;

    return 1;
}

int sph2bmul(int *f, float *i, float *j, const image* img, const float *v, int b){
    int w,h;
    const float X = fabsf(v[0]);
    const float Y = fabsf(v[1]);
    const float Z = fabsf(v[2]);
    float x,y;

    if(b==1){
	if     (v[2] < 0 && Z >= X && Z >= Y) { *f = 0; x =  v[0] / Z; y = -v[1] / Z; w = img[0].w; h = img[0].h;} // pFront
	else if(v[0] > 0 && X >= Y && X >= Z) { *f = 1; x =  v[2] / X; y = -v[1] / X; w = img[1].w; h = img[1].h;} // pRight
	else if(v[2] > 0 && Z >= X && Z >= Y) { *f = 2; x = -v[0] / Z; y = -v[1] / Z; w = img[2].w; h = img[2].h;} // pBack
	else if(v[0] < 0 && X >= Y && X >= Z) { *f = 3; x = -v[2] / X; y = -v[1] / X; w = img[3].w; h = img[3].h;} // pLeft
	else if(v[1] > 0 && Y >= X && Y >= Z) { *f = 4; x =  v[0] / Y; y = -v[2] / Y; w = img[4].w; h = img[4].h;} // pTop
	else if(v[1] < 0 && Y >= X && Y >= Z) { *f = 5; x =  v[0] / Y; y =  v[2] / Y; w = img[5].w; h = img[5].h;} // pBottom
    }
    else if(b==2){
	if(Y<Z && Y<X){
	    if     (v[2] <0 && v[0]>0)        { *f = 0; x =  v[0] / Z; y = -v[1] / Z; w = img[0].w; h = img[0].h;} // pFront
	    else if(v[2] >0 && v[0]>0)        { *f = 1; x =  v[2] / X; y = -v[1] / X; w = img[1].w; h = img[1].h;} // pRight
	    else if(v[2] >0 && v[0]<0)        { *f = 2; x = -v[0] / Z; y = -v[1] / Z; w = img[2].w; h = img[2].h;} // pBack
	    else if(v[2] <0 && v[0]<0)        { *f = 3; x = -v[2] / X; y = -v[1] / X; w = img[3].w; h = img[3].h;} // pLeft
	}
	else if(v[1]>0)                       { *f = 4; x =  v[0] / Y; y = -v[2] / Y; w = img[4].w; h = img[4].h;} // pTop
	else if(v[1]<0)                       { *f = 5; x =  v[0] / Y; y =  v[2] / Y; w = img[5].w; h = img[5].h;} // pBottom
    }
    else if (b==3){
	if(Y<Z && Y<X){
	    if     (v[2] <0 && v[0]>0)        { *f = 1; x =  v[2] / X; y = -v[1] / X; w = img[1].w; h = img[1].h;} // pRight
	    else if(v[2] >0 && v[0]>0)        { *f = 2; x = -v[0] / Z; y = -v[1] / Z; w = img[2].w; h = img[2].h;} // pBack
	    else if(v[2] >0 && v[0]<0)        { *f = 3; x = -v[2] / X; y = -v[1] / X; w = img[3].w; h = img[3].h;} // pLeft
	    else if(v[2] <0 && v[0]<0)        { *f = 0; x =  v[0] / Z; y = -v[1] / Z; w = img[0].w; h = img[0].h;} // pFront

	}
	else if(v[2]<0 && X<Z)                { *f = 0; x =  v[0] / Z; y = -v[1] / Z; w = img[0].w; h = img[0].h;} // pFront
	else if(v[0]>0 && X>Z)                { *f = 1; x =  v[2] / X; y = -v[1] / X; w = img[1].w; h = img[1].h;} // pRight
	else if(v[2]>0 && X<Z)                { *f = 2; x = -v[0] / Z; y = -v[1] / Z; w = img[2].w; h = img[2].h;} // pBack
	else if(v[0]<0 && X>Z)                { *f = 3; x = -v[2] / X; y = -v[1] / X; w = img[3].w; h = img[3].h;} // pLeft  
    }
    else return 0;

    *i = h*(y+(1.+G_OV*2.)) / (2.+G_OV*4.);
    *j = w*(x+(1.+G_OV*2.)) / (2.+G_OV*4.);

    return 1;
}

int sph2trec(int *f, float *i, float *j, const image* img, const float *v, int b){
    int w,h;
    float phi   = acosf(v[1])/PI;
    float theta = atan2f(v[0],-v[2]);
    float pmin  = 0.; 

    if     (phi <= 1./6) { *f = 0; w = img[0].w; h = img[0].h; pmin=0./6;} // pTop
    else if(phi <= 2./6) { *f = 1; w = img[1].w; h = img[1].h; pmin=1./6;} // p1
    else if(phi <= 3./6) { *f = 2; w = img[2].w; h = img[2].h; pmin=2./6;} // p2
    else if(phi <= 4./6) { *f = 3; w = img[3].w; h = img[3].h; pmin=3./6;} // p3
    else if(phi <= 5./6) { *f = 4; w = img[4].w; h = img[4].h; pmin=4./6;} // p4
    else if(phi <= 6./6) { *f = 5; w = img[5].w; h = img[5].h; pmin=5./6;} // pBot
    else return 0;

    *i = h * (phi-pmin)/(1./6);
    *j = w * (0.5f + theta / TWOPI);

    return 1;
}

float aBlendBmul(const float *v){
    const float X = fabsf(v[0]);
    const float Y = fabsf(v[1]);
    const float Z = fabsf(v[2]);
    float x,y;

    if(Y<Z && Y<X){
	if     (v[2] <0 && v[0]>0 && X/Z>(1-2*G_OV) && X/Z<(1+2*G_OV)){return 1-(X/Z-(1-2*G_OV))/(4*G_OV);} //1,2
	else if(v[2] >0 && v[0]>0 && X/Z>(1-2*G_OV) && X/Z<(1+2*G_OV)){return   (X/Z-(1-2*G_OV))/(4*G_OV);} //2,3
	else if(v[2] >0 && v[0]<0 && X/Z>(1-2*G_OV) && X/Z<(1+2*G_OV)){return 1-(X/Z-(1-2*G_OV))/(4*G_OV);} //3,4
	else if(v[2] <0 && v[0]<0 && X/Z>(1-2*G_OV) && X/Z<(1+2*G_OV)){return   (X/Z-(1-2*G_OV))/(4*G_OV);} //4,1
    }
    else if(v[1]>0 && Y/Z>(1-2*G_OV) && Y/Z<(1+2*G_OV) && Z/X>1){return (Y/Z-(1-2*G_OV))/(4*G_OV);} //5,0,2
    else if(v[1]>0 && Y/X>(1-2*G_OV) && Y/X<(1+2*G_OV) && Z/X<1){return (Y/X-(1-2*G_OV))/(4*G_OV);} //5,1,3
    else if(v[1]<0 && Y/X>(1-2*G_OV) && Y/X<(1+2*G_OV) && Z/X<1){return (Y/X-(1-2*G_OV))/(4*G_OV);} //6,0,2
    else if(v[1]<0 && Y/Z>(1-2*G_OV) && Y/Z<(1+2*G_OV) && Z/X>1){return (Y/Z-(1-2*G_OV))/(4*G_OV);} //6,1,3
    return 1.;
}

float aBlendBrec(const float *v){
    float phi   = acosf(v[1])/PI;
    if      (phi > (1-G_OV)/6 && phi < (1+G_OV)/6) {return 1.-(phi-(1-G_OV)/6)/(2*G_OV)*6;}
    else if (phi > (2-G_OV)/6 && phi < (2+G_OV)/6) {return 1.-(phi-(2-G_OV)/6)/(2*G_OV)*6;}
    else if (phi > (3-G_OV)/6 && phi < (3+G_OV)/6) {return 1.-(phi-(3-G_OV)/6)/(2*G_OV)*6;}
    else if (phi > (4-G_OV)/6 && phi < (4+G_OV)/6) {return 1.-(phi-(4-G_OV)/6)/(2*G_OV)*6;}
    else if (phi > (5-G_OV)/6 && phi < (5+G_OV)/6) {return 1.-(phi-(5-G_OV)/6)/(2*G_OV)*6;}
    return 1.;
}

float aBlendGrid(const float *v){
    float phi   = acosf(v[1])/PI;
    float tht   = atan2f(v[0],-v[2])/TWOPI+0.5f;
    if      (phi > (1-G_OV)/6 && phi < (1+G_OV)/6) {return 1.-(phi-(1-G_OV)/6)/(2*G_OV)*6;}
    else if (phi > (2-G_OV)/6 && phi < (2+G_OV)/6) {return 1.-(phi-(2-G_OV)/6)/(2*G_OV)*6;}
    else if (phi > (3-G_OV)/6 && phi < (3+G_OV)/6) {return 1.-(phi-(3-G_OV)/6)/(2*G_OV)*6;}
    else if (phi > (4-G_OV)/6 && phi < (4+G_OV)/6) {return 1.-(phi-(4-G_OV)/6)/(2*G_OV)*6;}
    else if (phi > (5-G_OV)/6 && phi < (5+G_OV)/6) {return 1.-(phi-(5-G_OV)/6)/(2*G_OV)*6;}

    else if (tht > (1-G_OVX)/4 && tht < (1+G_OVX)/4) {return 1.-(tht-(1-G_OVX)/4)/(2*G_OVX)*4;}
    else if (tht > (2-G_OVX)/4 && tht < (2+G_OVX)/4) {return 1.-(tht-(2-G_OVX)/4)/(2*G_OVX)*4;}
    else if (tht > (3-G_OVX)/4 && tht < (3+G_OVX)/4) {return 1.-(tht-(3-G_OVX)/4)/(2*G_OVX)*4;}
    return 1.;
}

float aBlendBeqr(const float *v){
    float phi   = acosf(v[1]) - PIBY2;
    float y     = (0.5+sinf(phi)/2.);
    if      (y > (1-G_OV)/6 && y < (1+G_OV)/6) {return 1.-(y-(1-G_OV)/6)/(2*G_OV)*6;}
    else if (y > (2-G_OV)/6 && y < (2+G_OV)/6) {return 1.-(y-(2-G_OV)/6)/(2*G_OV)*6;}
    else if (y > (3-G_OV)/6 && y < (3+G_OV)/6) {return 1.-(y-(3-G_OV)/6)/(2*G_OV)*6;}
    else if (y > (4-G_OV)/6 && y < (4+G_OV)/6) {return 1.-(y-(4-G_OV)/6)/(2*G_OV)*6;}
    else if (y > (5-G_OV)/6 && y < (5+G_OV)/6) {return 1.-(y-(5-G_OV)/6)/(2*G_OV)*6;}
    return 1.;
}

int sph2brec(int *f, float *i, float *j, const image* img, const float *v, int b){
    int w,h;
    float phi   = acosf(v[1])/PI;
    float theta = atan2f(v[0],-v[2]);
    float pmin  = 0.; 
    float pspn  = (1+2*G_OV)/6;

    if(b==1){
	if     (phi <= 1./6) { *f = 0; w = img[0].w; h = img[0].h; pmin=0.  ;       pspn=(1+G_OV)/6;} // pTop
	else if(phi <= 2./6) { *f = 1; w = img[1].w; h = img[1].h; pmin=(1-G_OV)/6;} // p1
	else if(phi <= 3./6) { *f = 2; w = img[2].w; h = img[2].h; pmin=(2-G_OV)/6;} // p2
	else if(phi <= 4./6) { *f = 3; w = img[3].w; h = img[3].h; pmin=(3-G_OV)/6;} // p3
	else if(phi <= 5./6) { *f = 4; w = img[4].w; h = img[4].h; pmin=(4-G_OV)/6;} // p4
	else if(phi <= 6./6) { *f = 5; w = img[5].w; h = img[5].h; pmin=(5-G_OV)/6; pspn=(1+G_OV)/6;} // pBot
    }
    else if(b==2){
	if     (phi <= (1+G_OV)/6) { *f = 0; w = img[0].w; h = img[0].h; pmin=0.;         pspn=(1+G_OV)/6;} // pTop
	else if(phi <= (2+G_OV)/6) { *f = 1; w = img[1].w; h = img[1].h; pmin=(1-G_OV)/6;} // p1
	else if(phi <= (3+G_OV)/6) { *f = 2; w = img[2].w; h = img[2].h; pmin=(2-G_OV)/6;} // p2
	else if(phi <= (4+G_OV)/6) { *f = 3; w = img[3].w; h = img[3].h; pmin=(3-G_OV)/6;} // p3
	else if(phi <= (5+G_OV)/6) { *f = 4; w = img[4].w; h = img[4].h; pmin=(4-G_OV)/6;} // p4
	else if(phi <= 1.)         { *f = 5; w = img[5].w; h = img[5].h; pmin=(5-G_OV)/6; pspn=(1+G_OV)/6;} // pBot
    }
    else if(b==3){
	if     (phi <= (1-G_OV)/6) { *f = 0; w = img[0].w; h = img[0].h; pmin=0.;         pspn=(1+G_OV)/6;} // pTop
	else if(phi <= (2-G_OV)/6) { *f = 1; w = img[1].w; h = img[1].h; pmin=(1-G_OV)/6;} // p1
	else if(phi <= (3-G_OV)/6) { *f = 2; w = img[2].w; h = img[2].h; pmin=(2-G_OV)/6;} // p2
	else if(phi <= (4-G_OV)/6) { *f = 3; w = img[3].w; h = img[3].h; pmin=(3-G_OV)/6;} // p3
	else if(phi <= (5-G_OV)/6) { *f = 4; w = img[4].w; h = img[4].h; pmin=(4-G_OV)/6;} // p4
	else if(phi <= 1.)         { *f = 5; w = img[5].w; h = img[5].h; pmin=(5-G_OV)/6; pspn=(1+G_OV)/6;} // pBot
    }
    else return 0;

    *i = h * (phi-pmin)/pspn;
    *j = w * (0.5f + theta / TWOPI);

    return 1;
}

int sph2grid(int *f, float *i, float *j, const image* img, const float *v, int b){
    int w,h;
    float phi   = acosf(v[1])/PI;
    float tht   = atan2f(v[0],-v[2])/TWOPI+0.5f;
    float pmin  = 0.; 
    float pspn  = (1+2*G_OV)/6;
    float tmin  = 0.;
    float tspn  = (1+2*G_OVX)/4;

    if(b==1){
	if     (phi <= 1./6){ 
	    pmin = 0.;
	    if      (tht<=1./4){ *f = 0; tmin=0.  ;}
	    else if (tht<=2./4){ *f = 1; tmin=(1.-G_OVX)/4;}
	    else if (tht<=3./4){ *f = 2; tmin=(2.-G_OVX)/4;}
	    else if (tht<=1.  ){ *f = 3; tmin=(3.-G_OVX)/4;}
	}
	else if(phi <= 2./6){ 
	    pmin = (1.-G_OV)/6.;
	    if      (tht<=1./4){ *f = 4; tmin=0.  ;}
	    else if (tht<=2./4){ *f = 5; tmin=(1.-G_OVX)/4;}
	    else if (tht<=3./4){ *f = 6; tmin=(2.-G_OVX)/4;}
	    else if (tht<=1.  ){ *f = 7; tmin=(3.-G_OVX)/4;}
	}
	else if(phi <= 3./6){ 
	    pmin = (2.-G_OV)/6.;
	    if      (tht<=1./4){ *f = 8; tmin=0.  ;}
	    else if (tht<=2./4){ *f = 9; tmin=(1.-G_OVX)/4;}
	    else if (tht<=3./4){ *f = 10; tmin=(2.-G_OVX)/4;}
	    else if (tht<=1.  ){ *f = 11; tmin=(3.-G_OVX)/4;}
	}
	else if(phi <= 4./6){ 
	    pmin = (3.-G_OV)/6.;
	    if      (tht<=1./4){ *f = 12; tmin=0.  ;}
	    else if (tht<=2./4){ *f = 13; tmin=(1.-G_OVX)/4;}
	    else if (tht<=3./4){ *f = 14; tmin=(2.-G_OVX)/4;}
	    else if (tht<=1.  ){ *f = 15; tmin=(3.-G_OVX)/4;}
	}
	else if(phi <= 5./6){ 
	    pmin = (4.-G_OV)/6.;
	    if      (tht<=1./4){ *f = 16; tmin=0.  ;}
	    else if (tht<=2./4){ *f = 17; tmin=(1.-G_OVX)/4;}
	    else if (tht<=3./4){ *f = 18; tmin=(2.-G_OVX)/4;}
	    else if (tht<=1.  ){ *f = 19; tmin=(3.-G_OVX)/4;}
	}
	else if(phi <= 1.){ 
	    pmin = (5.-G_OV)/6.;
	    if      (tht<=1./4){ *f = 20; tmin=0.  ;}
	    else if (tht<=2./4){ *f = 21; tmin=(1.-G_OVX)/4;}
	    else if (tht<=3./4){ *f = 22; tmin=(2.-G_OVX)/4;}
	    else if (tht<=1.  ){ *f = 23; tmin=(3.-G_OVX)/4;}
	}

    }
    else if(b==2){
	if     (phi >=(1.-G_OV)/6 && phi <= (1.+G_OV)/6){ 
	    pmin = 0.;
	    if      (tht<=1./4){ *f = 0; tmin=0.  ;}
	    else if (tht<=2./4){ *f = 1; tmin=(1.-G_OVX)/4;}
	    else if (tht<=3./4){ *f = 2; tmin=(2.-G_OVX)/4;}
	    else if (tht<=1.  ){ *f = 3; tmin=(3.-G_OVX)/4;}
	}
	else if(phi >=(1.-G_OV)/6 && phi <= (2.+G_OV)/6){ 
	    pmin = (1.-G_OV)/6.;
	    if      (tht<=1./4){ *f = 4; tmin=0.  ;}
	    else if (tht<=2./4){ *f = 5; tmin=(1.-G_OVX)/4;}
	    else if (tht<=3./4){ *f = 6; tmin=(2.-G_OVX)/4;}
	    else if (tht<=1.  ){ *f = 7; tmin=(3.-G_OVX)/4;}
	}
	else if(phi >=(3.-G_OV)/6 && phi <= (3.+G_OV)/6){ 
	    pmin = (2.-G_OV)/6.;
	    if      (tht<=1./4){ *f = 8; tmin=0.  ;}
	    else if (tht<=2./4){ *f = 9; tmin=(1.-G_OVX)/4;}
	    else if (tht<=3./4){ *f = 10; tmin=(2.-G_OVX)/4;}
	    else if (tht<=1.  ){ *f = 11; tmin=(3.-G_OVX)/4;}
	}
	else if(phi >=(4.-G_OV)/6 && phi <= (4.+G_OV)/6){ 
	    pmin = (3.-G_OV)/6.;
	    if      (tht<=1./4){ *f = 12; tmin=0.  ;}
	    else if (tht<=2./4){ *f = 13; tmin=(1.-G_OVX)/4;}
	    else if (tht<=3./4){ *f = 14; tmin=(2.-G_OVX)/4;}
	    else if (tht<=1.  ){ *f = 15; tmin=(3.-G_OVX)/4;}
	}
	else if(phi >=(5.-G_OV)/6 && phi <= (5.+G_OV)/6){ 
	    pmin = (4.-G_OV)/6.;
	    if      (tht<=1./4){ *f = 16; tmin=0.  ;}
	    else if (tht<=2./4){ *f = 17; tmin=(1.-G_OVX)/4;}
	    else if (tht<=3./4){ *f = 18; tmin=(2.-G_OVX)/4;}
	    else if (tht<=1.  ){ *f = 19; tmin=(3.-G_OVX)/4;}
	}
	else if(tht >= (1.-G_OVX)/4 && tht <= (1.+G_OVX)/4){
	    tmin = 0.;
	    if      (phi<=1./6){ *f = 0; pmin = 0.;} 
	    else if (phi<=2./6){ *f = 4; pmin = (1.-G_OV)/6;}
	    else if (phi<=3./6){ *f = 8; pmin = (2.-G_OV)/6;}
	    else if (phi<=4./6){ *f = 12;pmin = (3.-G_OV)/6;}
	    else if (phi<=5./6){ *f = 16;pmin = (4.-G_OV)/6;}
	    else if (phi<=1.  ){ *f = 20;pmin = (5.-G_OV)/6;}
	}
	else if(tht >= (2.-G_OVX)/4 && tht <= (2.+G_OVX)/4){
	    tmin = (1.-G_OVX)/4;
	    if      (phi<=1./6){ *f = 1; pmin = 0.;} 
	    else if (phi<=2./6){ *f = 5; pmin = (1.-G_OV)/6;}
	    else if (phi<=3./6){ *f = 9; pmin = (2.-G_OV)/6;}
	    else if (phi<=4./6){ *f = 13;pmin = (3.-G_OV)/6;}
	    else if (phi<=5./6){ *f = 17;pmin = (4.-G_OV)/6;}
	    else if (phi<=1.  ){ *f = 21;pmin = (5.-G_OV)/6;}
	}
	else if(tht >= (3.-G_OVX)/4 && tht <= (3.+G_OVX)/4){
	    tmin = (2.-G_OVX)/4;
	    if      (phi<=1./6){ *f = 2; pmin = 0.;} 
	    else if (phi<=2./6){ *f = 6; pmin = (1.-G_OV)/6;}
	    else if (phi<=3./6){ *f = 10;pmin = (2.-G_OV)/6;}
	    else if (phi<=4./6){ *f = 14;pmin = (3.-G_OV)/6;}
	    else if (phi<=5./6){ *f = 18;pmin = (4.-G_OV)/6;}
	    else if (phi<=1.  ){ *f = 22;pmin = (5.-G_OV)/6;}
	}
    }
    else if(b==3){
	if     ((phi >= (1.-G_OV)/6) && (phi <= (1.+G_OV)/6)){ 
	    pmin = (1.-G_OV)/6;
	    if      (tht<=1./4){ *f = 4; tmin=0.;}
	    else if (tht<=2./4){ *f = 5; tmin=(1.-G_OVX)/4;}
	    else if (tht<=3./4){ *f = 6; tmin=(2.-G_OVX)/4;}
	    else if (tht<=1.  ){ *f = 7; tmin=(3.-G_OVX)/4;}
	}
	else if(phi >=(2.-G_OV)/6 && phi <= (2.+G_OV)/6){ 
	    pmin = (2.-G_OV)/6.;
	    if      (tht<=1./4){ *f = 8; tmin=0.  ;}
	    else if (tht<=2./4){ *f = 9; tmin=(1.-G_OVX)/4;}
	    else if (tht<=3./4){ *f = 10; tmin=(2.-G_OVX)/4;}
	    else if (tht<=1.  ){ *f = 11; tmin=(3.-G_OVX)/4;}
	}
	else if(phi >=(3-G_OV)/6 && phi <= (3.+G_OV)/6){ 
	    pmin = (3.-G_OV)/6.;
	    if      (tht<=1./4){ *f = 12; tmin=0.  ;}
	    else if (tht<=2./4){ *f = 13; tmin=(1.-G_OVX)/4;}
	    else if (tht<=3./4){ *f = 14; tmin=(2.-G_OVX)/4;}
	    else if (tht<=1.  ){ *f = 15; tmin=(3.-G_OVX)/4;}
	}
	else if(phi >=(4-G_OV)/6 && phi <= (4.+G_OV)/6){ 
	    pmin = (4.-G_OV)/6.;
	    if      (tht<=1./4){ *f = 16; tmin=0.  ;}
	    else if (tht<=2./4){ *f = 17; tmin=(1.-G_OVX)/4;}
	    else if (tht<=3./4){ *f = 18; tmin=(2.-G_OVX)/4;}
	    else if (tht<=1.  ){ *f = 19; tmin=(3.-G_OVX)/4;}
	}
	else if(phi >=(5-G_OV)/6 && phi <= (5.+G_OV)/6){ 
	    pmin = (5.-G_OV)/6.;
	    if      (tht<=1./4){ *f = 20; tmin=0.  ;}
	    else if (tht<=2./4){ *f = 21; tmin=(1.-G_OVX)/4;}
	    else if (tht<=3./4){ *f = 22; tmin=(2.-G_OVX)/4;}
	    else if (tht<=1.  ){ *f = 23; tmin=(3.-G_OVX)/4;}
	}
	else if(tht >= (1.-G_OVX)/4 && tht <= (1.+G_OVX)/4){
	    tmin = (1.-G_OVX)/4;
	    if      (phi<=1./6){ *f = 1; pmin = 0.;} 
	    else if (phi<=2./6){ *f = 5; pmin = (1.-G_OV)/6;}
	    else if (phi<=3./6){ *f = 9; pmin = (2.-G_OV)/6;}
	    else if (phi<=4./6){ *f = 13;pmin = (3.-G_OV)/6;}
	    else if (phi<=5./6){ *f = 17;pmin = (4.-G_OV)/6;}
	    else if (phi<=1.  ){ *f = 21;pmin = (5.-G_OV)/6;}
	}
	else if(tht >= (2.-G_OVX)/4 && tht <= (2.+G_OVX)/4){
	    tmin = (2.-G_OVX)/4;
	    if      (phi<=1./6){ *f = 2; pmin = 0.;} 
	    else if (phi<=2./6){ *f = 6; pmin = (1.-G_OV)/6;}
	    else if (phi<=3./6){ *f = 10;pmin = (2.-G_OV)/6;}
	    else if (phi<=4./6){ *f = 14;pmin = (3.-G_OV)/6;}
	    else if (phi<=5./6){ *f = 18;pmin = (4.-G_OV)/6;}
	    else if (phi<=1.  ){ *f = 22;pmin = (5.-G_OV)/6;}
	}
	else if(tht >= (3.-G_OVX)/4 && tht <= (3.+G_OVX)/4){
	    tmin = (3.-G_OVX)/4;
	    if      (phi<=1./6){ *f = 3; pmin = 0.;} 
	    else if (phi<=2./6){ *f = 7; pmin = (1.-G_OV)/6;}
	    else if (phi<=3./6){ *f = 11;pmin = (2.-G_OV)/6;}
	    else if (phi<=4./6){ *f = 15;pmin = (3.-G_OV)/6;}
	    else if (phi<=5./6){ *f = 19;pmin = (4.-G_OV)/6;}
	    else if (phi<=1.  ){ *f = 23;pmin = (5.-G_OV)/6;}
	}
    }
    else return 0;

    if      (*f ==0 || *f==4 || *f==8  || *f == 12 || *f == 16 || *f == 20)
	tspn = (1.+G_OVX)/4;
    else if (*f== 3 || *f ==7|| *f==11 || *f == 15 || *f == 19 || *f == 23)
	tspn = (1.+G_OVX)/4;
    if      (*f==0 <=3 || *f >=20)
	pspn = (1.+G_OV)/6;

    w = img[*f].w;
    h = img[*f].h;

    *i = h * (phi-pmin)/pspn;
    *j = w * (tht-tmin)/tspn;

    return 1;
}

int sph2beqr(int *f, float *i, float *j, const image* img, const float *v, int b){
    int w,h;
    float phi   = acosf(v[1]) - PIBY2;
    float y     = (0.5+sinf(phi)/2.);
    float theta = atan2f(v[0],-v[2]);
    float pmin  = 0.; 
    float pspn  = (1+2*G_OV)/6;

    if(b==1){
	if     (y <= 1./6) { *f = 0; w = img[0].w; h = img[0].h; pmin=0.  ;       pspn=(1+G_OV)/6;} // pTop
	else if(y <= 2./6) { *f = 1; w = img[1].w; h = img[1].h; pmin=(1-G_OV)/6;} // p1
	else if(y <= 3./6) { *f = 2; w = img[2].w; h = img[2].h; pmin=(2-G_OV)/6;} // p2
	else if(y <= 4./6) { *f = 3; w = img[3].w; h = img[3].h; pmin=(3-G_OV)/6;} // p3
	else if(y <= 5./6) { *f = 4; w = img[4].w; h = img[4].h; pmin=(4-G_OV)/6;} // p4
	else if(y <= 6./6) { *f = 5; w = img[5].w; h = img[5].h; pmin=(5-G_OV)/6; pspn=(1+G_OV)/6;} // pBot
    }
    else if(b==2){
	if     (y <= (1+G_OV)/6) { *f = 0; w = img[0].w; h = img[0].h; pmin=0.;         pspn=(1+G_OV)/6;} // pTop
	else if(y <= (2+G_OV)/6) { *f = 1; w = img[1].w; h = img[1].h; pmin=(1-G_OV)/6;} // p1
	else if(y <= (3+G_OV)/6) { *f = 2; w = img[2].w; h = img[2].h; pmin=(2-G_OV)/6;} // p2
	else if(y <= (4+G_OV)/6) { *f = 3; w = img[3].w; h = img[3].h; pmin=(3-G_OV)/6;} // p3
	else if(y <= (5+G_OV)/6) { *f = 4; w = img[4].w; h = img[4].h; pmin=(4-G_OV)/6;} // p4
	else if(y <= 1.)         { *f = 5; w = img[5].w; h = img[5].h; pmin=(5-G_OV)/6; pspn=(1+G_OV)/6;} // pBot
    }
    else if(b==3){
	if     (y <= (1-G_OV)/6) { *f = 0; w = img[0].w; h = img[0].h; pmin=0.;         pspn=(1+G_OV)/6;} // pTop
	else if(y <= (2-G_OV)/6) { *f = 1; w = img[1].w; h = img[1].h; pmin=(1-G_OV)/6;} // p1
	else if(y <= (3-G_OV)/6) { *f = 2; w = img[2].w; h = img[2].h; pmin=(2-G_OV)/6;} // p2
	else if(y <= (4-G_OV)/6) { *f = 3; w = img[3].w; h = img[3].h; pmin=(3-G_OV)/6;} // p3
	else if(y <= (5-G_OV)/6) { *f = 4; w = img[4].w; h = img[4].h; pmin=(4-G_OV)/6;} // p4
	else if(y <= 1.)         { *f = 5; w = img[5].w; h = img[5].h; pmin=(5-G_OV)/6; pspn=(1+G_OV)/6;} // pBot
    }
    else return 0;

    *i = h * (y-pmin)/pspn;
    *j = w * (0.5f + theta / TWOPI);

    return 1;
}


int sph2teqr(int *f, float *i, float *j, const image* img, const float *v, int b){
    int w,h;
    float phi   = acosf(v[1]) - PIBY2;
    float y     = (0.5+sinf(phi)/2.);
    float theta = atan2f(v[0],-v[2]);
    float pmin  = 0.; 

    if     (y <= 1./6) { *f = 0; w = img[0].w; h = img[0].h; pmin=0./6;} // pTop
    else if(y <= 2./6) { *f = 1; w = img[1].w; h = img[1].h; pmin=1./6;} // p1
    else if(y <= 3./6) { *f = 2; w = img[2].w; h = img[2].h; pmin=2./6;} // p2
    else if(y <= 4./6) { *f = 3; w = img[3].w; h = img[3].h; pmin=3./6;} // p3
    else if(y <= 5./6) { *f = 4; w = img[4].w; h = img[4].h; pmin=4./6;} // p4
    else if(y <= 6./6) { *f = 5; w = img[5].w; h = img[5].h; pmin=5./6;} // pBot
    else return 0;

    *i = h * (y-pmin)/(1./6);
    *j = w * (0.5f + theta / TWOPI);

    return 1;
}


int sph2cube(int *f, float *i, float *j, const image* img, const float *v, int b){
    int w = img->w;
    int h = img->h;
    
    const float X = fabsf(v[0]);
    const float Y = fabsf(v[1]);
    const float Z = fabsf(v[2]);
    float x;
    float y;
    
    if      (v[0] > 0 && X >= Y && X >= Z) { *f = 0; x = -v[2] / X; y = -v[1] / X; }
    else if (v[0] < 0 && X >= Y && X >= Z) { *f = 1; x =  v[2] / X; y = -v[1] / X; }
    else if (v[1] > 0 && Y >= X && Y >= Z) { *f = 2; x =  v[0] / Y; y =  v[2] / Y; }
    else if (v[1] < 0 && Y >= X && Y >= Z) { *f = 3; x =  v[0] / Y; y = -v[2] / Y; }
    else if (v[2] > 0 && Z >= X && Z >= Y) { *f = 4; x =  v[0] / Z; y = -v[1] / Z; }
    else if (v[2] < 0 && Z >= X && Z >= Y) { *f = 5; x = -v[0] / Z; y = -v[1] / Z; }
    else return 0;

    *i = 1.0f + (h - 2) * (y + 1.0f) / 2.0f;
    *j = 1.0f + (w - 2) * (x + 1.0f) / 2.0f;

    return 1;
}


// [mapping]2sph
int eqar2sph(int f, float i, float j, int h, int w, float *v){
    const float y = (h/2.0f-i)/(h/2.0f);
    const float phi = asinf(y);
    const float theta = TWOPI * j / w - PI;

    v[0] =  sinf(theta) * cosf(phi);
    v[1] =                sinf(phi);
    v[2] = -cosf(theta) * cosf(phi);

    return 1;
}

int teqr2sph(int f, float i, float j, int h, int w, float *v){
    const float phi   = asinf(1-2.*(G_MAP_SF_Y*(i/h)+G_MAP_OFF_Y));
    const float theta = TWOPI * (G_MAP_SF_X*(j/w)+G_MAP_OFF_X) - PI;

    v[0] =  sinf(theta) * cosf(phi);
    v[1] =                sinf(phi);
    v[2] = -cosf(theta) * cosf(phi);

    return 1;
}

int rect2sph(int f, float i, float j, int h, int w, float *v){
    const float lat = PIBY2 - PI * i / h;
    const float lon = TWOPI * j / w - PI;

    v[0] =  sinf(lon) * cosf(lat);
    v[1] =              sinf(lat);
    v[2] = -cosf(lon) * cosf(lat);

    return 1;
}

int trec2sph(int f, float i, float j, int h, int w, float *v){
    const float lat = PIBY2 - PI * (G_MAP_SF_Y*(i / h)+G_MAP_OFF_Y);
    const float lon = TWOPI * (G_MAP_SF_X*(j / w)+G_MAP_OFF_X) - PI;

    v[0] =  sinf(lon) * cosf(lat);
    v[1] =              sinf(lat);
    v[2] = -cosf(lon) * cosf(lat);

    return 1;
}

int brec2sph(int f, float i, float j, int h, int w, float *v){
    return trec2sph(f,i,j,h,w,v);
}

int grid2sph(int f, float i, float j, int h, int w, float *v){
    return trec2sph(f,i,j,h,w,v);
}

int beqr2sph(int f, float i, float j, int h, int w, float *v){
    return teqr2sph(f,i,j,h,w,v);
}

int merc2sph(int f, float i, float j, int h, int w, float *v){
    const float phi = atan(sinh(3-6.0*i/h));
    const float theta = TWOPI * j / w - PI;

    v[0] =  sinf(theta) * cosf(phi);
    v[1] =                sinf(phi);
    v[2] = -cosf(theta) * cosf(phi);

    return 1;
}

int dyad2sph(int f, float i, float j, int h, int w, float *v){
    float lat = PIBY2 - PI * i / h;
    float lon = TWOPI * j / w - PI;

    if (fabsf(lat) <=PI/3.0)
	return rect2sph(f,i,j,h,w,v);
    else if (lon > 0){
	v[0] = 0;
	v[1] = -1;
	v[2] = 0;
	return 1;
    }
    return rect2sph(f,i,j,h,w/2,v);
}

int cube2sph(int f, float i, float j, int h, int w, float *c){
    const int p[6][3][3] = {
        {{  0,  0, -1 }, {  0, -1,  0 }, {  1,  0,  0 }},
        {{  0,  0,  1 }, {  0, -1,  0 }, { -1,  0,  0 }},
        {{  1,  0,  0 }, {  0,  0,  1 }, {  0,  1,  0 }},
        {{  1,  0,  0 }, {  0,  0, -1 }, {  0, -1,  0 }},
        {{  1,  0,  0 }, {  0, -1,  0 }, {  0,  0,  1 }},
        {{ -1,  0,  0 }, {  0, -1,  0 }, {  0,  0, -1 }},
    };

    const float v = 2.0f * i / h - 1.0f;
    const float u = 2.0f * j / w - 1.0f;

    int wxx = p[f][0][0], wxy = p[f][1][0], cx = p[f][2][0];
    int wyx = p[f][0][1], wyy = p[f][1][1], cy = p[f][2][1];
    int wzx = p[f][0][2], wzy = p[f][1][2], cz = p[f][2][2];
    
    float x = wxx * u + wxy * v + cx;
    float y = wyx * u + wyy * v + cy;
    float z = wzx * u + wzy * v + cz;

    c[0] = x; c[1] = y; c[2] = z;
    normalize(c);
    
    return 1;
}

int cos22sph(int f, float i, float j, int h, int w, float *v){
    const float y = (h/2.0f-i)/(h/2.0f);
    const float phi = LUTinv(y);
    const float theta = TWOPI * j / w - PI;

    v[0] =  sinf(theta) * cosf(phi);
    v[1] =                sinf(phi);
    v[2] = -cosf(theta) * cosf(phi);

    return 1;
}

int view2sph(int f, float v, float u, int h, int w, float *p){
    // normalized device coordinates: G_iK * [u;v;1]
    float x2 = G_iK[0][0]*u + G_iK[0][1]*v + G_iK[0][2];
    float y2 = G_iK[1][0]*u + G_iK[1][1]*v + G_iK[1][2];

    // undo perspective division
    float z1 = 1/sqrt(x2*x2+y2*y2+1);
    float x1 = z1*x2;
    float y1 = z1*y2;

    // flip the z-axis
    float p1[3] = {x1, y1, -z1};

    // rotate: p = R * p1
    p[0] = G_R[0][0]*p1[0] + G_R[0][1]*p1[1] + G_R[0][2]*p1[2];
    p[1] = G_R[1][0]*p1[0] + G_R[1][1]*p1[1] + G_R[1][2]*p1[2];
    p[2] = G_R[2][0]*p1[0] + G_R[2][1]*p1[1] + G_R[2][2]*p1[2];

    return 1;
}

// ==============================================================================
//                               INTERPOLATION
// ==============================================================================
void filter_nearest(const image *img, const image *acs, float i, float j, float *p)
{
    const float ii = clamp(i - 0.5f, 0.0f, img->h - 1.0f);
    const float jj = clamp(j - 0.5f, 0.0f, img->w - 1.0f);

    const long  i0 = lrintf(ii);
    const long  j0 = lrintf(jj);

    int k;

    for (k = 0; k < img->c; k++){
        p[k] += img->p[(img->w * i0 + j0) * img->c + k];
    }

    // set accessed pixels
    if (G_ACSFLAG)
        for (k = 0; k < img->c; k++)
            acs->p[(img->w * i0 + j0) * img->c + k] = 1;
}

void filter_linear(const image *img, const image *acs, float i, float j, float *p){
    const float ii = clamp(i - 0.5f, 0.0f, img->h - 1.0f);
    const float jj = clamp(j - 0.5f, 0.0f, img->w - 1.0f);

    const long  i0 = lrintf(floorf(ii)), i1 = lrintf(ceilf(ii));
    const long  j0 = lrintf(floorf(jj)), j1 = lrintf(ceilf(jj));

    const float di = ii - i0;
    const float dj = jj - j0;

    int k;
    for (k = 0; k < img->c; k++)
    {
        p[k] += lerp(lerp(img->p[(img->w * i0 + j0) * img->c + k],
                          img->p[(img->w * i0 + j1) * img->c + k], dj),
                     lerp(img->p[(img->w * i1 + j0) * img->c + k],
                          img->p[(img->w * i1 + j1) * img->c + k], dj), di);

        // update access locations
        if (G_ACSFLAG){
            acs->p[(img->w * i0 + j0) * img->c + k] += (1-di)*(1-dj);
            acs->p[(img->w * i0 + j1) * img->c + k] += (1-di)*dj;
            acs->p[(img->w * i1 + j0) * img->c + k] += di*(1-dj);
            acs->p[(img->w * i1 + j1) * img->c + k] += di*dj;
        }
    }
    return;
}

void filter_bicubic(const image *img, const image *acs, float i, float j, float *p){
    // assume subpixel centroid of (0.5, 0.5) in src image
    i -= 0.5; j -= 0.5;

    // take care of borders by doing linear interpolation
    const int row_LB = 1;
    const int row_HB = img->w-2;
    const int col_LB = 1;
    const int col_HB = img->h-2;

    if( (i <= col_LB) || (i >= col_HB) || (j <= row_LB) || (j >= row_HB) ){
        filter_linear(img, acs, i+0.5, j+0.5, p);
        return;
    }

    // FIND ROWS / COLS
    const long iZ = lrintf(floorf(i)-1);
    const long i0 = iZ+1;
    const long i1 = i0+1;
    const long i2 = i1+1;
    const long jZ = lrintf(floorf(j))-1;
    const long j0 = jZ+1;
    const long j1 = j0+1;
    const long j2 = j1+1;

    // bicubic interpolate in the x direction, then the y direction
    int k;
    for (k=0; k < img->c; k++){
        p[k] += cubicInterp(cubicInterp(img->p[(img->w*iZ+jZ)*img->c+k],
					img->p[(img->w*iZ+j0)*img->c+k],
					img->p[(img->w*iZ+j1)*img->c+k],
					img->p[(img->w*iZ+j2)*img->c+k],j-j0),
			    cubicInterp(img->p[(img->w*i0+jZ)*img->c+k],
					img->p[(img->w*i0+j0)*img->c+k],
					img->p[(img->w*i0+j1)*img->c+k],
					img->p[(img->w*i0+j2)*img->c+k],j-j0),
			    cubicInterp(img->p[(img->w*i1+jZ)*img->c+k],
					img->p[(img->w*i1+j0)*img->c+k],
					img->p[(img->w*i1+j1)*img->c+k],
					img->p[(img->w*i1+j2)*img->c+k],j-j0),
			    cubicInterp(img->p[(img->w*i2+jZ)*img->c+k],
					img->p[(img->w*i2+j0)*img->c+k],
					img->p[(img->w*i2+j1)*img->c+k],
					img->p[(img->w*i2+j2)*img->c+k],j-j0),i-i0);

        // update access locations
        if (G_ACSFLAG){
            float di = i-i0;
            float dj = j-j0;

            acs->p[(img->w*i0+j0)*img->c+k] += (1-di)*(1-dj);
            acs->p[(img->w*i0+j1)*img->c+k] += (1-di)*dj;
            acs->p[(img->w*i1+j0)*img->c+k] += di*(1-dj);
            acs->p[(img->w*i1+j1)*img->c+k] += di*dj;
        }
    }

    return;
}


