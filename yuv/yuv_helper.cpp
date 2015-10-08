#include "yuv_helper.h"

/* init
 * ------------------------------------
 * Descrp: initializes the yuv base class
 */
void cYuv::init(std::vector<int> ws, std::vector<int> hs, int n){
   // alloc memory
    m_Ws=ws; m_Hs=hs; m_F=n;
    int w = *std::max_element(m_Ws.begin(),m_Ws.end());
    int h = *std::max_element(m_Hs.begin(),m_Hs.end());
    if ( n == 1 ){
        m_pY = img_alloc(  w,  h,n);
        m_pU = img_alloc(w/2,h/2,n);
        m_pV = img_alloc(w/2,h/2,n);
    }
    else if ( n == 6 ){
        // alloc with border
        m_pY = img_alloc(  w+2,  h+2,6);
        m_pU = img_alloc(w/2+2,h/2+2,6);
        m_pV = img_alloc(w/2+2,h/2+2,6);
    }
    else {
        fputs ("Num faces error", stderr);
        return;
    }
}

/* toNthYuv
 * -------------------------------------------
 * Descrp: Generates a numbered yuv file given a base yuv file
 *         E.g., test.yuv -> test-001.yuv
 */
std::string cYuv::toNthYuv(const char* fbase, int n){
    std::string s = std::string(fbase);
    int idx = s.find(".");
    if (idx <0)
	return std::string("");

    s.resize(idx);
    char buf[10];
    sprintf(buf,"-%03d.yuv",n);
    s+= std::string(buf);
    return s;
}

/* yuvFileOpen
 * ----------------------------------
 * Descrp: Use to potentially open 1 or more yuv files for reading or writing
 */
void cYuv::yuvFileOpen(const char* fname, int n, bool readFlag){
    // just one file to read/write
    if(!multFlag){
	if(readFlag)
	    m_pFiles.push_back(fopen(fname,"rb"));
	else
	    m_pFiles.push_back(fopen(fname,"w+"));
	if(m_pFiles[0] == NULL)
	    fputs ("File read error\n", stderr);
    }
    // multiple files
    else{
	if(readFlag){
	    for(int i=0;i<n;i++){
		m_pFiles.push_back(fopen(toNthYuv(fname,i).c_str(),"rb"));
		if(m_pFiles[i]==NULL)
		    fputs ("File read error\n", stderr);
	    }
	}
    }
}

/* setDimensions
 * -------------------------------------
 * Descrp: Sets the img dimensions of the preallocated imgs
 */
void cYuv::setDimensions(image* img, std::vector<int> ws, std::vector<int>hs, float sf, int pad){
    for(int i=0;i<ws.size();i++){
	img[i].w = ws[i]*sf+pad;
	img[i].h = hs[i]*sf+pad;
    }
}

/* init
 * ------------------------------------
 * Descrp: initializes the yuv reader
 */
void cYuvReader::init(const char *fname, std::vector<int> ws, std::vector<int> hs, int n, bool mFlag)
{
    // set parameters
    multFlag = mFlag;
    m_Ws=ws; m_Hs=hs; m_F=n;
    int w = *std::max_element(m_Ws.begin(),m_Ws.end());
    int h = *std::max_element(m_Hs.begin(),m_Hs.end());
    yuvFileOpen(fname,n,true);
    nFiles = 1;
    if(multFlag)
	nFiles = n;
    int pad = 0;

    // alloc memory
    if ( n == 1 || multFlag){
        m_pY = img_alloc(  w,  h,n);
        m_pU = img_alloc(w/2,h/2,n);
        m_pV = img_alloc(w/2,h/2,n);
    }
    // cube special case
    else if ( n == 6 ) {
        // temp for reading
        m_pTmpY = img_alloc(  w,  h,6);
        m_pTmpU = img_alloc(w/2,h/2,6);
        m_pTmpV = img_alloc(w/2,h/2,6);

        // alloc with border
        m_pY = img_alloc(  w+2,  h+2,6);
        m_pU = img_alloc(w/2+2,h/2+2,6);
        m_pV = img_alloc(w/2+2,h/2+2,6);
	
	pad = 2;
    }
    else {
        fputs ("Num faces error", stderr);
        return;
    }

    // alloc buffer for fread
    long lSize = n*(w*h + w*h/2);
    m_pBuffer = (uint8_t*) malloc (sizeof(uint8_t)*lSize);

    // set image sizes
    setDimensions(m_pY,m_Ws,m_Hs,  1,pad);
    setDimensions(m_pU,m_Ws,m_Hs,0.5,pad);
    setDimensions(m_pV,m_Ws,m_Hs,0.5,pad);
}


/* init
 * ------------------------------------
 * Descrp: initializes the yuv writer
 */
void cYuvWriter::init(const char *fname, int w, int h, int n)
{
    m_Ws.push_back(w); m_Hs.push_back(h); m_F=n;
    yuvFileOpen(fname,n,false);

    // alloc memory
    m_pY = img_alloc(  w,  h,n);
    m_pU = img_alloc(w/2,h/2,n);
    m_pV = img_alloc(w/2,h/2,n);
}

/* readNextFrame
 * -----------------------------------------
 * Descrp: Wrapper to read frames from all sub files
 */
bool cYuvReader::readNextFrame(){
    if (m_pFiles[0] == NULL) 
	return false;

    // read next frame for each file
    for(int i=0; i < nFiles; i++){
	if(!readNextFrame(i))
	    return false;
    }
    return true;
}

/* readNextFrame
 * ------------------------------------
 * Descrp: Reads an input YUV file and generates 3 images representing
 *         the respective Y,U,V channels
 * Inputs: (char*) YUV file
 *         (int)   image width
 *         (int)   image height
 *         (int)   num subimages (e.g., 6 in cube)
 * Return: (image*) Y image
 */
bool cYuvReader::readNextFrame(int pID){
    // how many bytes do we need to read per file
    long wh    = m_Ws[pID]*m_Hs[pID];
    long lSize = wh + wh/2;
    if(nFiles ==1)
	lSize*=m_F;

    // legacy
    int m_W = m_Ws[0];
    int m_H = m_Hs[0];
    
    // read data into the buffer
    size_t result = fread(m_pBuffer,1,(size_t)lSize,m_pFiles[pID]);
    if(result != lSize)
	return false;

    // generate img from buffer data (assuming YUV420)    
    if (nFiles > 1 || m_F==1){
	buffer2ImgMRQ(m_pY[pID], m_pBuffer,   wh,       0);
	buffer2ImgMRQ(m_pU[pID], m_pBuffer, wh/4,      wh);
	buffer2ImgMRQ(m_pV[pID], m_pBuffer, wh/4, wh+wh/4);
    }
    if (m_F == 6 && nFiles == 1) {
        buffer2Img(m_pTmpY, m_pBuffer,  m_W,  m_H,                      0,m_F);
        buffer2Img(m_pTmpU, m_pBuffer,m_W/2,m_H/2,            m_W*m_H*m_F,m_F);
        buffer2Img(m_pTmpV, m_pBuffer,m_W/2,m_H/2,(m_W*m_H+m_W*m_H/4)*m_F,m_F);

        cCubeHelper cubeHelper;
        cubeHelper.image_border(m_pTmpY, m_pY);
        cubeHelper.image_border(m_pTmpU, m_pU);
        cubeHelper.image_border(m_pTmpV, m_pV);
    }
    return true;
}

/* writeNextFrame
 * ------------------------------------
 * Descrp: Writes yuv data to file in append mode
 */
bool cYuvWriter::writeNextFrame(bool color)
{
    if (m_pFiles[0] == NULL) return false;

    if (color)
    {
        image *Y = getY();
        for(int f=0; f<m_F; f++)
            for(int i=0; i<Y->w*Y->h; i++)
                fprintf(m_pFiles[0],"%c",(unsigned int)(Y[f].p[i]*255 + 0.499));

        image *U = getU();
        for(int f=0; f<m_F; f++)
            for(int i=0; i<U->w*U->h; i++)
                fprintf(m_pFiles[0],"%c",(unsigned int)(U[f].p[i]*255 + 0.499));

        image *V = getV();
        for(int f=0; f<m_F; f++)
            for(int i=0; i<V->w*V->h; i++)
                fprintf(m_pFiles[0],"%c",(unsigned int)(V[f].p[i]*255 + 0.499));
    }
    else
    {
        image *Y = getY();
        for(int f=0; f<m_F; f++)
            for(int i=0; i<Y->w*Y->h; i++)
                fprintf(m_pFiles[0],"%c",(unsigned int)(Y[f].p[i]*16));

        image *U = getU();
        for(int f=0; f<m_F; f++)
            for(int i=0; i<U->w*U->h; i++)
                fprintf(m_pFiles[0],"%c",(unsigned int)(128));

        image *V = getV();
        for(int f=0; f<m_F; f++)
            for(int i=0; i<V->w*V->h; i++)
                fprintf(m_pFiles[0],"%c",(unsigned int)(128));
    }

    return true;
}


/* buffer2ImgMRQ
 * -------------------------------
 * Descrp: Version2 reads a certain amount from the buffer and puts it in the data of the resulting img
 */
void cYuvReader::buffer2ImgMRQ(image& img, uint8_t* bufferIn, long nRead, long offset){
    float* bufferOut = img.p;
    for(size_t i=0; i<nRead;i++)
	bufferOut[i] = 1/255.0f*(float)bufferIn[offset+i];
}

/* buffer2Img
 * ----------------------------------
 * Descrp:
 * Inputs:
 * Return:
 */
void cYuvReader::buffer2Img(image* img, uint8_t* bufferIn, int w, int h, int offset, int n)
{
    for(int f=0;f<n;f++){
        float* bufferOut = img[f].p;
        for(size_t i=0; i < w*h; i++)
            bufferOut[i] =1/255.0f*(float)bufferIn[offset+f*w*h+i];
    }
}

/*
 * img_alloc
 * ------------------------------------
 * Descrp:
 * Inputs:
 * Return:
 */
image* cYuv::img_alloc(int w, int h, int n)
{
    image* img = 0;
    img = (image *) calloc(n, sizeof (image));
    for (int i=0; i <n; i++)
    {
        img[i].p = (float *) calloc(w * h, sizeof (float));
        img[i].w = w;
        img[i].h = h;
        img[i].c = 1;
        img[i].b = 0;
        img[i].s = 0;
    }
    return img;
}

/*
 * img_free
 * ------------------------------------
 * Descrp:
 * Inputs:
 * Return:
 */
void cYuv::img_free(image *img)
{
    for (int i=0; i <m_F; i++)
    {
        if (img[i].p) free(img[i].p);
    }
    if (img) free(img);
}

/*
 * clear
 * ------------------------------------
 * Descrp:
 * Inputs:
 * Return:
 */
void cYuv::clear()
{
    for (int i=0; i <m_F; i++)
    {
        memset(m_pY[i].p, 0, m_pY[i].w*m_pY[i].h*sizeof(float));
        memset(m_pU[i].p, 0, m_pU[i].w*m_pU[i].h*sizeof(float));
        memset(m_pV[i].p, 0, m_pV[i].w*m_pV[i].h*sizeof(float));
    }
}

/*
 * setY
 * -------------------------------
 * Descrp:
 * Inputs:
 * Return:
 */
void cYuv::setY(float* ydata)
{
    for(int i=0; i<m_pY->w*m_pY->h;i++)
	m_pY->p[i]=ydata[i];
}



/*
 * setU
 * -------------------------------
 * Descrp:
 * Inputs:
 * Return:
 */
void cYuv::setU(float* udata)
{
    for(int i=0; i<m_pU->w*m_pU->h;i++){
	m_pU->p[i]=udata[i];
    }
}

/*
 * setV
 * -------------------------------
 * Descrp:
 * Inputs:
 * Return:
 */
void cYuv::setV(float* vdata)
{
    for(int i=0; i<m_pV->w*m_pV->h;i++){
	m_pV->p[i]=vdata[i];
    }
}

//--------------------------------------------
float* rotN(image *img, int i, int j)
{
    const int ii = i;
    const int jj = j;
    return img->p + img->c * (img->w * ii + jj);
}

float* rotL(image *img, int i, int j)
{
    const int ii = j;
    const int jj = img->h - i - 1;
    return img->p + img->c * (img->w * ii + jj);
}

float* rotR(image *img, int i, int j)
{
    const int ii = img->w - j - 1;
    const int jj = i;
    return img->p + img->c * (img->w * ii + jj);
}
//--------------------------------------------

/* Add borders to a cubemap image. Assume the given image pointer is an array */
/* of six images. Copy each to a new sef of six images, each two pixels wider */
/* and higher. Also copy the borders. This is necessary for correct cubemap   */
/* sampling.                                                                  */
void cCubeHelper::image_border(image *src, image *dst)
{
    const int d = 1;

    if ((src) && (dst))
    {
        const int n = src[0].w;
        const int c = src[0].c;
        const int b = 4;

        const int N = n + 2 * d;

        /* Copy all page data. */
        blit(dst[0].p, N, d, d, src[0].p, n, 0, 0, n, n, c, b);
        blit(dst[1].p, N, d, d, src[1].p, n, 0, 0, n, n, c, b);
        blit(dst[2].p, N, d, d, src[2].p, n, 0, 0, n, n, c, b);
        blit(dst[3].p, N, d, d, src[3].p, n, 0, 0, n, n, c, b);
        blit(dst[4].p, N, d, d, src[4].p, n, 0, 0, n, n, c, b);
        blit(dst[5].p, N, d, d, src[5].p, n, 0, 0, n, n, c, b);
	
        border(dst + 0, rotN, dst + 5, rotN, d);
        border(dst + 5, rotN, dst + 1, rotN, d);
        border(dst + 1, rotN, dst + 4, rotN, d);
        border(dst + 4, rotN, dst + 0, rotN, d);

        border(dst + 1, rotR, dst + 2, rotN, d);
        border(dst + 1, rotL, dst + 3, rotN, d);

        border(dst + 2, rotN, dst + 0, rotL, d);
        border(dst + 3, rotN, dst + 0, rotR, d);

        border(dst + 2, rotL, dst + 4, rotL, d);
        border(dst + 2, rotR, dst + 5, rotL, d);
        border(dst + 3, rotL, dst + 5, rotR, d);
        border(dst + 3, rotR, dst + 4, rotR, d);
    }
    return;
}

void cCubeHelper::border(image *a, rot rota, image *b, rot rotb, int d)
{
    const size_t s = b->c * sizeof (float);
    const int    n = b->h;

    for     (int i = d; i < n - d; i++)
        for (int j = 0; j <     d; j++)
        {
            memcpy(rota(a, i, n - d + j), rotb(b, i,         d + j), s);
            memcpy(rotb(b, i,         j), rota(a, i, n - d - d + j), s);
        }
}

void cCubeHelper::blit(void *dst, int dw, int dx, int dy,
		       void *src, int sw, int sx, int sy,
		       int w, int h, int c, int b)
{
    int i;
    for (i = 0; i < h; i++)
        memcpy( (char *) dst + ((i + dy) * dw + dx) * c * b,
                (char *) src + ((i + sy) * sw + sx) * c * b,
                w * c * b);
}

