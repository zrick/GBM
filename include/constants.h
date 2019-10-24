//
//  constants.h
//  GBM
//
//  Created by Cedrick Ansorge on 30.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#ifndef constants_h
#define constants_h

#define GBM_FILE_ERR "gbm.err"
#define GBM_FILE_LOG "gbm.log"
#define GBM_FILE_OUT "gbm.out"

// Hard-coded array sizes for memory compactness of triangulation
#define MAX_CHAR_LEN 512

#define MAX_VRT_TTR 64
#define MAX_VRT_TRI 64
#define MAX_VRT_EDG 64

#define MAX_EDG_TRI 32
#define MAX_EDG_TTR 32

#define MAX_TTR_ADJ 64

#define GBM_NBOX 5

// OUTPUT FORMATS
#define FMT_ASC_BASIC 0
#define FMT_ASC_EXT   1
#define FMT_ASC_FULL  2
#define FMT_XML_VTK   3

#define PI 3.141592653589793

#define PERIODIC_TOLERANCE 1e-16

#define PERIODIC_NONE -1
#define PERIODIC_X   1
#define PERIODIC_Y   2
#define PERIODIC_Z   4
#define PERIODIC_XY  3
#define PERIODIC_XZ  5
#define PERIODIC_YZ  6
#define PERIODIC_XYZ 7



#define SPLINE_NONE     -1
#define SPLINE_O1       0
#define SPLINE_O2_LOCAL 1
#define SPLINE_O2_EDGES 2
#define SPLINE_O2_CLSTR 3
#define MAX_SPLINES     4

#define SPLINE_O2_LOCAL_POINTS  { {0,0},{1,1},{2,2},{3,3},{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
#define SPLINE_O2_EDGES_POINTS  { {4,4},{5,5},{6,6},{7,7},{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};

#define GBMDER_NONE 0
#define GBMDER_X    1
#define GBMDER_Y    2
#define GBMDER_Z    3
#define GBMDER_XX   4
#define GBMDER_YY   5
#define GBMDER_ZZ   6
#define GBMDER_XY   7
#define GBMDER_XZ   8
#define GBMDER_YZ   9

#endif /* constants_h */
