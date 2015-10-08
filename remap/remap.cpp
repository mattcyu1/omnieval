#include <getopt.h>
#include "yuv_helper.h"
#include "map_utils.h"
#include "panomapper.h"

static int usage(const char *exe){
    fprintf(stderr,
            "%s [-i input] [-o output] [-f filter] [-m m] [-n n] [-w w] [-h h] [-t tf] [-y] src dst\n"
            "\t-i ... Input  file type: cube, rect, eqar, merc                   [rect]\n"
            "\t-o ... Output file type: cube, rect, eqar, merc, view             [rect]\n"
            "\t-f ... Filter type: nearest, linear, bicubic                   [bicubic]\n"
            "\t-m ... Input  height list                                          [500]\n"
            "\t-b ... Input  width                                                 [2m]\n"
            "\t-n ... Output height                                               [500]\n"
            "\t-v ... Output width                                                 [2n]\n"
            "\t-w ... Viewport width                                              [200]\n"
            "\t-h ... Viewport height                                             [200]\n"
            "\t-x ... Viewport fov x in degree                                     [90]\n"
            "\t-y ... Viewport fov y in degree                                     [90]\n"
            "\t-p ... Viewport center position phi (degrees)                        [0]\n"
            "\t-l ... Viewport center position tht (degrees)                        [0]\n"
	    "\t-t ... Tracking data file                                         [none]\n"
	    "\t-y ... Blend data together (only works with orec, etc ...)         [off]\n"
            "\t-z ... Number of frames                                            [MAX]\n",
            exe);
    return 0;
}

int main(int argc, char **argv){
    // check command line inputs
    int c,n=500,z=__INT_MAX__,v=-1;
    float x=90.,y=90.,w=200,h=200;
    float p=0.,l=0.;
    bool  blendFlag = false;
    const char *i=NULL, *o=NULL, *f=NULL, *t=NULL, *a=NULL,*m=NULL, *b=NULL;
    while ((c = getopt(argc, argv, "i:o:m:n:z:t:x:y:w:h:f:a:b:v:p:l:r")) != -1){
        switch (c){
            case 'i': i = optarg;                      break;
	    case 'o': o = optarg;                      break;
            case 't': t = optarg;                      break;
            case 'f': f = optarg;                      break;
            case 'a': a = optarg;                      break;
            case 'm': m = optarg;                      break;
            case 'n': n = (int)strtol(optarg, 0, 0);   break;
            case 'z': z = (int)strtol(optarg, 0, 0);   break;
            case 'x': x = (float)strtof(optarg, 0);    break;
            case 'y': y = (float)strtof(optarg, 0);    break;
            case 'w': w = (float)strtof(optarg, 0);    break;
            case 'h': h = (float)strtof(optarg, 0);    break;
            case 'b': b = optarg;                      break;
            case 'v': v = (int)strtol(optarg, 0, 0);   break;
            case 'p': p = (float)strtol(optarg, 0, 0); break;
	    case 'l': l = (float)strtol(optarg, 0, 0); break;
	    case 'r': blendFlag = true;                break;
            default : return usage(argv[0]);
        }
    }
    if(argc<=2)
	return usage(argv[0]);

    // mapping / yuv read write init
    remapper pm;    
    pm.init(i,o,f,m,blendFlag,n,x,y,w,h,z,b,v,p,l,t,a,argv[optind],argv[optind+1]);
    pm.remapFrames();

    // map all frames
    return 0;
}
