/*=================================================================
 *
 * Backproject.c
 * 
 * This .mex file performs the linear and nearest-neighbour
 * interpolation versions of the backprojection in iradon.m.
 * It does the same operation as the following piece of
 * code taken from iradon.m
 *
 * For nearest-neighbour interpolation
 *         for i=1:length(theta)   
 *            proj = p(:,i);
 *            t = round(x*costheta(i) + y*sintheta(i));
 *            img = img + proj(t+ctrIdx);
 *         end
 *
 * For linear interpolation
 *         for i=1:length(theta)
 *             proj = p(:,i);
 *             t = x.*costheta(i) + y.*sintheta(i);
 *             a = floor(t);
 *             img = img + (t-a).*proj(a+1+ctrIdx) + (a+1-t).*proj(a+ctrIdx);
 *         end
 *
 * Each loop can be replaced with one of
 *
 *         img = Backproject(p, costheta, sintheta, N, 0); % NN-interp
 *         img = Backproject(p, costheta, sintheta, N, 1); % linear interp
 *
 * The input parameters are:
 *
 *         p          Radon transform, with projections stored in columns
 *                    (as output by radon.m, for example)
 *         costheta   a vector (column or row) holding the cosine of the
 *                    projections angles
 *         sintheta   a vector (column or row) holding the sine of the
 *                    projections angles
 *         N          the output image with be NxN (this is a double float)
 *         interp     0 => nearest-neighbour interpolation
 *                    1 => linear interpolation
 *
 * The output is:
 *
 *         img        an NxN (double float) matrix
 *
 * WARNING: This code has not been exhaustively tested.  Use at your
 *          own risk. (standard disclaimer)
 *
 * Written by:
 *   Jeff Orchard
 *   Assistant Professor, David R. Cheriton School of Computer Science
 *   University of Waterloo
 *   Waterloo, Ontario, CANADA
 *   e-mail: jorchard@cs.uwaterloo.ca
 *   Web: http://www.cs.uwaterloo.ca/~jorchard/
 *
 * November 1, 2006
 *
 * This is a MEX-file for MATLAB.  
 * Copyright 1984-2004 The MathWorks, Inc.
 *
 *=================================================================*/
/* $Revision: 1.0 $ */
#include <math.h>
#include "mex.h"



void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{
	/* INPUT PARAMETERS */
    double *p;			/* filtered projections (one per col) */
	double *costheta;	/* pre-computed cos of proj angles */
	double *sintheta;	/* pre-computed sin of proj angles */
	double *Nptr;		/* size of reconstructed image */
	double *interp_ptr;	/* 0 for NN interp, 1 for linear interp */
	
	/* OUTPUT PARAMETERS */
	double *img;		/* output image */
	
	
	/* Other variables */
	int N;				/* integer copy of Nptr (above) */
	int interp_flag;	/* integer copy of interp_ptr (above) */
	int x_idx, y_idx, theta_idx;	/* loop indecies */
	double x;			/* stores x-coordinate of pixel */
	/*double y;*/		/* stores y-coordinate of pixel (not needed) */
	double ctr, xleft, ytop;	/* used to tranform from matrix indecies */
								/* to (x,y) coords (see iradon.m) */
	int len;		/* length of each projection (spatial dimension) */
	int n_proj;		/* number of projection angles (length of costheta vector, eg) */
	int ctr_idx;	/* centre index for projections */
	
	/* Temporary variable used for code optimization */
	double cos_theta, sin_theta, t;
	int ctr_idx_1;	/* = ctr_idx + 1 (for optimization purposes) */
	double x_cos_theta;	/* = x * cos_theta */
	int a;
	double *proj;	/* points at the start of a projection (a column) */
	double *img_ptr;
	
    
    /* Check for proper number of arguments */
    
    if (nrhs != 5) { 
		mexErrMsgTxt("Five input arguments required."); 
    } else if (nlhs > 1) {
		mexErrMsgTxt("Too many output arguments."); 
    } 
	
    p = mxGetPr(prhs[0]);
    costheta = mxGetPr(prhs[1]);
    sintheta = mxGetPr(prhs[2]);
    Nptr = mxGetPr(prhs[3]);
	interp_ptr = mxGetPr(prhs[4]);
	
	N = (int)(Nptr[0]);
	interp_flag = (int)(interp_ptr[0]);
	len = mxGetM(prhs[0]);
	n_proj = mxGetM(prhs[1]) * mxGetN(prhs[1]); /* use product to account for row- or col-vectors */
	
	/* The line below takes care of the difference between Matlab matrix indexing
	 * (base-1) and C-style matrix indexing (base-0).  The corresponding Matlab code is
	 *       ctr = floor((N + 1)/2);
	 * eg. if N=8
	 *          Matlab -> ctr = index of 4
	 *          C      -> ctr = index of 3
	 */
	ctr = floor((N-1) / 2);
	
	xleft = -ctr;
	ytop = ctr;
	
	ctr_idx = (int)floor(len/2);  /* centre index for projections */
	ctr_idx_1 = ctr_idx + 1;
    
    
    /* Create a matrix for the return argument */
    plhs[0] = mxCreateDoubleMatrix(N, N, mxREAL);
    img = mxGetPr(plhs[0]);  /* Get pointer to the data array */
	
	
    
    for (theta_idx=0 ; theta_idx<n_proj ; theta_idx++)
	{
	
		cos_theta = costheta[theta_idx];
		sin_theta = sintheta[theta_idx];
		proj = (p+theta_idx*len);  /* point at proper collumn */
		img_ptr = img;
		
		/* The loops below are designed to fly through the elements of the output image (img) in
		 * linear-memory sequence... makes it faster.  The matricies are stored column-by-column.
		 * The switch statement is outside the loops for optimization purposes.
		 */		
		
		switch (interp_flag)
		{
		
		case 0: /* nearest-neighbour interpolation */
		
			for (x_idx=0 ; x_idx<N ; x_idx++)
			{
				x = xleft + x_idx;  /* x-coord */
				x_cos_theta = x*cos_theta;  /* perform this multiplication outside the y-loop */
				t = x_cos_theta + ytop*sin_theta;  /* After this, t can simply be decremented by sin_theta each y-iter */
			
				for (y_idx=0 ; y_idx<N ; y_idx++)
				{
					/*y = ytop - y_idx; */
					/*t = x_cos_theta + y*sin_theta;*/  /* removed for optimization purposes */
				
					//a = round(t);
                    a = (int)(t);
				
					*img_ptr++ += proj[a+ctr_idx];
				
					t -= sin_theta;  /* decrement t */
				} /* end of y-loop */
			} /* end of x-loop */
			
			break;
			
		case 1: /* linear interpolation */
		
			for (x_idx=0 ; x_idx<N ; x_idx++)
			{
				x = xleft + x_idx;  /* x-coord */
				x_cos_theta = x*cos_theta;  /* perform this multiplication outside the y-loop */
				t = x_cos_theta + ytop*sin_theta;  /* After this, t can simply be decremented by sin_theta each y-iter */
			
				for (y_idx=0 ; y_idx<N ; y_idx++)
				{
					/*y = ytop - y_idx;
					/*t = x_cos_theta + y*sin_theta;*/  /* removed for optimization purposes */
				
					a = floor(t);
                    //a = int(t);
                    //if((t-a)>0.5 || (t-a)<(-0.5))
                      //  a = a+1;
				
					/* The following line is the slower version (2 multiplications) */
					/**img_ptr++ += (t-a)*proj[a+ctr_idx_1] + (1-t+a)*proj[a+ctr_idx]; */
				
					/* This is the speedy version (1 multiplication) */
					*img_ptr++ += (t-a)*( proj[a+ctr_idx_1] - proj[a+ctr_idx] ) + proj[a+ctr_idx];
				
					t -= sin_theta;  /* decrement t */
				} /* end of y-loop */
			} /* end of x-loop */
			
			break;
			
		} /* end of switch statement */
	} /* end of theta loop */
	
  	
    return;
    
}


