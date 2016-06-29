////////////////////////////////////////////////////////////////////////////////  
// Copyright (c) 2013, Lawrence Livermore National Security, LLC. 
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory. 
// Written by Erik Draeger (draeger1@llnl.gov) and Francois Gygi (fgygi@ucdavis.edu).
// Based on the Qbox code by Francois Gygi Copyright (c) 2008 
// LLNL-CODE-635376. All rights reserved. 
//
// qb@ll is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details, in the file COPYING in the
// root directory of this distribution or <http://www.gnu.org/licenses/>.
//
/*******************************************************************************
 *
 * spline.h
 *
 ******************************************************************************/

#define SPLINE_FLAT_BC 0.0       /* Flat boundary condition (y'=0) */
#define SPLINE_NATURAL_BC 1.e31  /* Natural boundary condition (Y"=0) */

void spline(double *x, double *y, int n, double yp1, double ypn, double *y2);
void splint (double *xa, double *ya, double *y2a, int n, double x, double *y);
void splintd (double *xa, double *ya, double *y2a, 
 int n, double x, double *y, double *dy);
