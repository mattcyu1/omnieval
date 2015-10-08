#ifndef _MAP_UTILS_H_
#define _MAP_UTILS_H_

#define PIBY2 1.5707963f
#define PI    3.1415927f
#define TWOPI 6.2831853f
#define PIBY3 1.0471976f
#define SIN60 0.8660254f

#include <math.h>
#include "../yuv/yuv_helper.h"

// ==============================================================================
//                               FUNCTION POINTERS
// ==============================================================================
typedef void  (*ifilter)   (const image *, const image*, float, float, float *);
typedef int   (*sph2map)   (int *, float *, float *, const image*, const float *, int);
typedef int   (*map2sph)   (int,   float,   float,   int, int,       float *);
typedef float (*aBlendMap) (const float *v);

// ==============================================================================
//                               STRUCTURES
// ==============================================================================
struct point{
    float i;
    float j;
};
typedef struct point point;

struct pattern{
    int    n;
    point *p;
};
typedef struct pattern pattern;

struct float2{
    float x;
    float y;
};
typedef struct float2 float2;

struct float3{
    float x;
    float y;
    float z;
};
typedef struct float3 float3;

struct int2{
    int x;
    int y;
};
typedef struct int2 int2;

// ==============================================================================
//                               SAMPLING PATTERNS
// ==============================================================================
extern point G_CENT_POINTS[];
const pattern G_CENT_PATTERN = {  1, G_CENT_POINTS};

// ==============================================================================
//                               INLINE HELPER FUNCTIONS
// ==============================================================================
inline float sign(float x){
    float s = (x > 0) ? 1 : ((x < 0) ? -1 : 0);
    return s;
}
inline float lerp(float a, float b, float k){
    return a * (1.f - k) + b * k;
}
inline float clamp(float f, float a, float z){
    if      (f < a) return a;
    else if (f > z) return z;
    else            return f;
}
inline void normalize(float *v) {
    const float k = 1.f / sqrtf(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    v[0] *= k;
    v[1] *= k;
    v[2] *= k;
}
inline float cubicInterp(float aZ, float a0, float a1, float a2, float t){
    float cZ =       2*a0;
    float c0 =  -aZ       +  a1;
    float c1 = 2*aZ -5*a0 +4*a1 -a2;
    float c2 =  -aZ+ 3*a0 -3*a1 +a2;
    
    float v = 0.5f*(cZ+c0*t+c1*t*t+c2*t*t*t);
    return clamp(v, 0, 1.0f);
}

// ==============================================================================
//                               MAPPING HELPERS
// ==============================================================================
extern float* G_INP_LUT;
extern float* G_OUT_LUT;
extern int    G_NUMPTS_LUT;
extern float  G_iK[3][3];
extern float  G_R[3][3];
extern bool   G_ACSFLAG;
extern float* G_LAT_INP_LUT;
extern float* G_LAT_OUT_LUT;
extern int    G_LAT_NUMPTS_LUT;
extern float  G_MAP_OFF_Y;
extern float  G_MAP_OFF_X;
extern float  G_MAP_SF_Y;
extern float  G_MAP_SF_X;

void setupLatLUT(const char* fname);
float latLUTLookup(float x);
void setupCos2();
void setMapOffset(float2 t);
void setMapScalefactor(float2 t);
void clearTrData();
void setRotationMat(float p, float t);
void setRotationMat(float *rt, float *up, float *lk);
int readNextTrData(FILE* ft, int nf, bool& bTrdata);
void setIntrMat(int dw, int dh, int sw, int sh, float fx, float fy);
int invert3x3(float K[3][3]);
float3 sphToCart(float2 sph);

// ==============================================================================
//                               MAPPINGS
// ==============================================================================
/*
 * sph2[mapping]
 * ---------------------------------------------
 * Descrp: Given 3D coordinate on sphere, calculates corresponding 
 *         point on mapped projection.
 *         This is generally used to get a destination value given a sph point
 * Inputs: (int*)          fill value address (increment to get multiple points)
 *         (float*)        mapped row
 *         (float*)        mapped column
 *         (int)           map height
 *         (int)           map width
 *         (const float*)  3d point
 *         (int)           indicator of which value to generate (based on blending)
 * Return: (int) success flag
 */
int sph2rect(int *f, float* i, float* j, const image* img, const float* v, int b);
int sph2eqar(int *f, float *i, float *j, const image* img, const float *v, int b);
int sph2dyad(int *f, float *i, float *j, const image* img, const float *v, int b);
int sph2cube(int *f, float *i, float *j, const image* img, const float *v, int b);
int sph2mult(int *f, float *i, float *j, const image* img, const float *v, int b);
int sph2bmul(int *f, float *i, float *j, const image* img, const float *v, int b);
int sph2trec(int *f, float *i, float *j, const image* img, const float *v, int b);
int sph2brec(int *f, float *i, float *j, const image* img, const float *v, int b);
int sph2grid(int *f, float *i, float *j, const image* img, const float *v, int b);
int sph2beqr(int *f, float *i, float *j, const image* img, const float *v, int b);
int sph2teqr(int *f, float *i, float *j, const image* img, const float *v, int b);
int sph2merc(int *f, float *i, float *j, const image* img, const float *v, int b);
int sph2cos2(int *f, float *i, float *j, const image* img, const float *v, int b);

/*
 * [mapping]2sph
 * ---------------------------------------------
 * Descrp: Given 2d coordinates on a mapped projection, calculates corresponding
 *         point on sphere.
 *         This is generally used to get the desired point from the destination
 * Inputs: (int*)          fill value
 *         (float*)        mapped row
 *         (float*)        mapped column
 *         (int)           map height
 *         (int)           map width
 *         (const float*)  3d point
 * Return: (int) success flag
 */
int eqar2sph(int f, float i, float j, int h, int w, float *v);
int rect2sph(int f, float i, float j, int h, int w, float *v);
int dyad2sph(int f, float i, float j, int h, int w, float *v);
int merc2sph(int f, float i, float j, int h, int w, float *v);
int cube2sph(int f, float i, float j, int h, int w, float *v);
int cos22sph(int f, float i, float j, int h, int w, float *v);
int view2sph(int f, float i, float j, int h, int w, float *v);
int trec2sph(int f, float i, float j, int h, int w, float *v);
int brec2sph(int f, float i, float j, int h, int w, float *v);
int grid2sph(int f, float i, float j, int h, int w, float *v);
int beqr2sph(int f, float i, float j, int h, int w, float *v);
int teqr2sph(int f, float i, float j, int h, int w, float *v);

// ==============================================================================
//                               INTERPOLATION
// ==============================================================================
/* 
 * aBlend[mapping]
 * -----------------------------------------------
 * Descrp: Generates the alpha blending factor used between two images
 *         Returns 1 if there should be no alpha blending
 * Inputs: (const float*) sphere location
 * Return: (float) alpha blend factor
 */
float aBlendBrec(const float *v);
float aBlendGrid(const float *v);
float aBlendBeqr(const float *v);
float aBlendBmul(const float *v);

// ==============================================================================
//                               INTERPOLATION
// ==============================================================================
/*
 * filter_[filtertype]
 * -----------------------------------------------
 * Descrp: Sample an image at row i column j using [filtertype] interpolation
 * Inputs: (const image*) mapped img
 *         (const image*) access img
 *         (float)        row
 *         (float)        column
 *         (float)        val
 * Return: (none)
 */
void filter_nearest(const image *img, const image *acs, float i, float j, float *p);
void filter_linear( const image *img, const image *acs, float i, float j, float *p);
void filter_bicubic(const image *img, const image *acs, float i, float j, float *p);

#endif
