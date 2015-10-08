#include <getopt.h>
#include "yuv_helper.h"
#include "map_utils.h"
#include "panomapper.h"

static int usage(const char *exe){
    fprintf(stderr,
            "%s [-i input1] [-o input2] [-f filter] [-m m] [-n n] [-z z] [-w w] [-s] src1 src2 sph\n"
            "\t-i ... Input file type: cube, rect  eqarea, merc, dyad          [rect]\n"
            "\t-o ... Input file type: cube, rect, eqarea, merc, dyad          [rect]\n"
            "\t-f ... Filter type: nearest, linear, bicubic                 [bicubic]\n"
            "\t-w ... Latitude weighting function                                 [1]\n"
            "\t-s ... Sphere weighting indicator                              [false]\n"
            "\t-m ... Src1 height                                               [500]\n"
            "\t-b ... Src1 width                                                 [2m]\n"
            "\t-n ... Src2 height                                               [500]\n"
            "\t-v ... Src1 width                                                 [2n]\n"
            "\t-z ... Num frames                                             [INTMAX]\n",
             exe);
    return 0;
}

int main(int argc, char **argv){
    // check cmd inputs
    int c,z=__INT_MAX__;
    bool swFlag = false, mserFlag=false;
    const char *i = NULL, *o=NULL, *f=NULL, *w=NULL, *m=NULL,*n=NULL, *b=NULL,*v=NULL;
    while ((c = getopt(argc, argv, "i:o:m:n:w:f:z:spv:b:")) != -1){
	switch (c){
	case 'i': i        = optarg;                    break;
	case 'o': o        = optarg;                    break;
	case 'f': f        = optarg;                    break;
	case 'w': w        = optarg;                    break;
	case 's': swFlag   = true;                      break;
	case 'p': mserFlag = true;                      break;
	case 'm': m        = optarg;                    break;
	case 'n': n        = optarg;                    break;
	case 'b': b        = optarg;                    break;
	case 'v': v        = optarg;                    break;
	case 'z': z        = (int)strtol(optarg, 0, 0); break;
	default: return usage(argv[0]);
        }
    }
    if(argc<=3)
	return usage(argv[0]);

    sphcomparer sc;
    sc.init(i,o,f,m,n,z,b,v,w,swFlag,argv[optind],argv[optind+1],argv[optind+2]);
    double t = sc.sphcomp(mserFlag);
    printf("PSNR: %.15f\n",t);

    return 0;
}
