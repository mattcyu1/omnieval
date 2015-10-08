#ifndef _YUV_HELPER_H_
#define _YUV_HELPER_H_

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <string>
#include <algorithm>

struct image
{
    float *p;  // data
    int    h;  // height
    int    w;  // width
    int    c;  // sample count
    int    b;  // sample depth
    int    s;  // sample format
};
typedef struct image image;

struct imageYUV
{
    image* Y;
    image* U;
    image* V;
};
typedef struct imageYUV imageYUV;

struct imgdim{
    int w;
    int h;
    int p;
};
typedef struct imgdim imgdim;

/* cYuv
 * ------------------------------------
 * Descrp: Base class for YUV reader/writer
 */
class cYuv
{
protected:
    std::vector<int> m_Ws; // width
    std::vector<int> m_Hs; // height
    int       m_F; // number of faces (rect:1, cube:6)
    int       nFiles; // number of files needed to read

    image    *m_pY, *m_pU, *m_pV;               // y,u,v pointers
    image    *m_pTmpY, *m_pTmpU, *m_pTmpV;      // temp y,u,v pointers

    std::vector<FILE*>    m_pFiles;
    image*    img_alloc(int w, int h, int n);
    void      img_free(image* img);
    bool      multFlag;
    void      setDimensions(image* img, std::vector<int> ws, std::vector<int> hs, float sf, int pad=0);

public:
    cYuv()
    {
	m_pY = NULL; m_pU = NULL; m_pV = NULL;
        m_pTmpY = NULL; m_pTmpU = NULL; m_pTmpV = NULL;
        m_Ws.clear(); m_Hs.clear(); m_F = 0; nFiles = 0;
    };
    ~cYuv()
    {
        if (m_pY) img_free(m_pY);
        if (m_pU) img_free(m_pU);
        if (m_pV) img_free(m_pV);

        if (m_pTmpY) img_free(m_pTmpY);
        if (m_pTmpU) img_free(m_pTmpU);
        if (m_pTmpV) img_free(m_pTmpV);
    }; 
    void init(std::vector<int> ws,std::vector<int> hs, int n);
    void clear();
    void setY(float* ydata);
    void setU(float* udata);
    void setV(float* vdata);
    void yuvFileOpen(const char* fbase, int n, bool readFlag);
    std::string toNthYuv(const char* fbase, int n);

    image* getY() { return m_pY; }
    image* getU() { return m_pU; }
    image* getV() { return m_pV; }

    std::vector<int> getW() { return m_Ws; }
    std::vector<int> getH() { return m_Hs; }
    int getF() { return m_F; }

};

/* cYuvReader
 * ------------------------------------
 * Descrp: Reads an input YUV420 file and generates 3 images representing
 *         the respective Y,U,V channels
 */
class cYuvReader: public cYuv
{
private:
    uint8_t *m_pBuffer;
    void buffer2Img(image* img, uint8_t* bufferIn, int w, int h, int offset, int n);
    void buffer2ImgMRQ(image& img, uint8_t* bufferIn, long nRead, long offset);

public:
    cYuvReader()  { m_pBuffer = NULL; multFlag = false;};
    ~cYuvReader() { if (m_pBuffer) free(m_pBuffer); };
    void init(const char* fname, std::vector<int> ws, std::vector<int> hs, int n, bool mFlag);
    bool readNextFrame();
    bool readNextFrame(int pID);
};

/* cYuvWriter
 * ------------------------------------
 * Descrp: Writes yuv data to file in append mode
 */
class cYuvWriter: public cYuv
{
public:
    cYuvWriter() {multFlag=false;};
    void init(const char* fname, int w, int h, int n);
    bool writeNextFrame(bool color=true);
};


/* cCubeHelper
 * ------------------------------------
 * Descrp: Helper class for cube map
 */
typedef float *(*rot)(image *, int, int);

class cCubeHelper
{
private:

    /* border
     * ------------------------------------------------------
     * Descrp: Adds a border from other pages to the current page
     * Inputs: (image*) img 1
     *         (rot)    rotation location 1
     *         (image*) img 2
     *         (rot)    rotation location 2
     *         (int)    n/a (adjustment to do with border width?)
     * Return:
     */
    void border(image *a, rot rota, image *b, rot rotb, int d);
    
    /* blit
     * -------------------------------------------------------
     * Descrp: copies img over from source to destination (while doing something?)
     * Inputs: (void*) destination data
     *         (int)   d width
     *         (int)   d row index
     *         (int)   d col index
     *         (void*) source data
     *         (int)   s width
     *         (int)   s row index
     *         (int)   s col index
     *         (int)   img width
     *         (int)   img height
     *         (int)   number of pages
     *         (int)   number of?
     * Return: (none)
     */
    void blit(void *dst, int dw, int dx, int dy,
              void *src, int sw, int sx, int sy,
              int w,     int h,  int c,  int b);

public:
    void image_border(image *src, image *dst);
};
#endif
