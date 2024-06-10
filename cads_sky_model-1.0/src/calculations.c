/*
  calculations.c
  sky_model

  Created by Jayant Murthy on 07/10/12.
  Copyright (c) 2012 Jayant Murthy. All rights reserved.
*/
#include "sky_model.h"

float CALC_SCALE_FACTOR(struct SPECTRA x, struct SPECTRA y)
{
    float scale_factor, int_lambda;
    float new_spec;
    int ix, iy;
    
    iy = 0;
    scale_factor = 0;
    int_lambda   = 0;
    
    /*Scale to the filter; Assume that they are all monotonic increasing.
      Counts = int(flux*effarea)/(int(effarea)*/
    for (ix = 0; ix < x.nelements; ++ix) {
        while (y.wavelength[iy] < x.wavelength[ix])
            ++iy;
        if (iy == 0)
            new_spec = 0;
        else
            new_spec = y.spectrum[iy - 1] + (y.spectrum[iy] - y.spectrum[iy - 1])/
            (y.wavelength[iy] - y.wavelength[iy - 1])*
            (x.wavelength[ix] - y.wavelength[iy - 1]);
        if (ix > 0)
            scale_factor += new_spec*x.spectrum[ix]*(x.wavelength[ix] -
                                                     x.wavelength[ix - 1]);
    }
    return (scale_factor);
}

