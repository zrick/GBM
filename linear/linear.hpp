//
//  linear.hpp
//  linear
//
//  Created by Cedrick Ansorge on 23.07.19.
//  Copyright © 2019 Cedrick Ansorge. All rights reserved.
//

#ifndef LINEAR_
#define LINEAR_

int qrdcmp(int nn);
int qrslv(int nn);

void crs3(double a[3], double b[3], double *c);
void dot3(double a[3], double b[3], double *c);
void renorm3(double *v);
#endif
