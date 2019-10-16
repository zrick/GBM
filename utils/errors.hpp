//
//  errors.hpp
//  GBM
//
//  Created by Cedrick Ansorge on 02.08.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#ifndef errors_hpp
#define errors_hpp

#include <exception>
#include "utils.hpp"

// any kind of error
#define GBMERR              1

// general issues
#define GBMERR_GENERAL   10
#define GBMERR_INDEX     11
#define GBMERR_DIM       12
#define GBMERR_NAMELIST  13
#define GBMERR_RESIDUAL  14

// grid-related problems
#define GBMERR_GRID      20
#define GBMERR_MAXEDGTRI 21
#define GBMERR_MAXEDGTTR 23
#define GBMERR_MAXEDGVRT 24
#define GBMERR_MAXVRTTRI 25
#define GBMERR_MAXTTRVRT 26
#define GBMERR_HALO      27
#define GBMERR_TOPOLOGY  28

// numerical problems / linear systems
#define GBMERR_NUMERIC     30
#define GBMERR_SINGULARITY 31

// misuse
#define GBMERR_MISUSE      90
#define GBMERR_PARAMETER   98
#define GBMERR_UNDEVELOPED 99
#endif /* errors_hpp */
