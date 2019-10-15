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

#define MAX_VRT_TTR 128
#define MAX_VRT_TRI 128
#define MAX_VRT_EDG 128

#define MAX_EDG_TRI 32
#define MAX_EDG_TTR 32

#define MAX_TTR_ADJ 64

// OUTPUT FORMATS
#define FMT_ASC_BASIC 0
#define FMT_ASC_EXT   1
#define FMT_ASC_FULL  2
#define FMT_XML_VTK   3

#define PI 3.141592653589793

#define PERIODIC_TOLERANCE 1e-16

#define PERIODIC_X   1
#define PERIODIC_Y   2
#define PERIODIC_Z   4
#define PERIODIC_XY  3
#define PERIODIC_XZ  5
#define PERIODIC_YZ  6
#define PERIODIC_XYZ 7

#endif /* constants_h */
