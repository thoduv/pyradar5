/*
 * This file implements a Python extension, wrapping the capabilities of the "RADAR5" FORTRAN code.
 */

#include <stdio.h>
#include <assert.h> 
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <Python.h>
#include "capsulethunk.h"

#define WITH_TCC

#ifdef WITH_TCC
#include "tcc/libtcc.h"
#endif

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

#ifdef _MSC_VER 
#	define GSL_DLL
#	define WIN32
#endif

#ifdef WITH_GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#else
#include "cspline.c"
#endif

#include "radar5.h"

#ifndef MIN
#define MIN(a,b) (((a)>(b))?(b):(a))
#endif  
#ifndef MAX
#define MAX(a,b) (((a)>(b))?(a):(b))
#endif

typedef void (*equ_c_t)(double *F, const double *Y, double X, double (*lag)(int, int, void*), void *ctx, void *params);
typedef double (*lag_c_t)(double X, const double *Y);

typedef struct
{
	double *X,*Y;
	ARGLAG_t ARGLAG;
	PHI_t PHI;
	double *RPAR;
	int *IPAR;
	double *PAST;
	int *IPAST, *NRDS;
	int lag_flag[10];
	int ipos[10];
	double alphas[10];
}  _lag_ctx;

typedef struct params_ params_t;
struct params_
{
	PyObject_HEAD
	
	int full_output;
	int xpos;
	double *xvalues;
	int xvalues_len;
	
	double *array;
	unsigned int array_size;
	unsigned int array_pos;
	unsigned int array_incr;
	
	
	PyObject *py_fcn;
	equ_c_t c_fcn;
	
	PyObject *py_y0[10];
	double constant_y0[10];
	PyObject *py_lagfuns[10];
	lag_c_t c_lagfuns[10];
	double constant_lags[10];
	PyThreadState *_thread_save;
	int interactive;
	int fail;
	void *user_params;
	
	PyObject *py_rpar;
	PyObject *py_lag_callback;

#if WITH_GSL
	gsl_interp_accel *spline_acc;
	gsl_spline **y0_spline;
#else
	cspline_t *y0_spline;
#endif
	double y0_t0;
	
	_lag_ctx *lag_ctx;
	
	params_t *next, *prev;
};

static double _lag_c(int il, int ip, _lag_ctx *ctx)
{
	int IL = il + 1, IP = ip + 1;
	
	/* Store the result of LAGR5 for a little speedup if a delay is used for more than one variable... */
	if(ctx->lag_flag[il] == 0)
	{
		lagr5_(&IL,ctx->X,ctx->Y,ctx->ARGLAG,ctx->PAST,&ctx->alphas[il],&ctx->ipos[il],ctx->RPAR,ctx->IPAR,ctx->PHI,ctx->IPAST,ctx->NRDS);   
		ctx->lag_flag[il] = 1;
	}
	
	/* Effective computation of the delayed value for this variable. */
	return ylagr5_(&IP,&ctx->alphas[il],&ctx->ipos[il],ctx->PHI,ctx->RPAR,ctx->IPAR,ctx->PAST,ctx->IPAST,ctx->NRDS);
}

static PyObject *_lag_py(PyObject *self, PyObject *args)
{
	int il, ip;
	params_t *p = PyCapsule_GetPointer(self, NULL);
	PyArg_ParseTuple(args, "ii", &il, &ip);
	return PyFloat_FromDouble(_lag_c(il, ip, p->lag_ctx));
}

static void FCN(int *N, double *X, double *Y, double *F, ARGLAG_t ARGLAG, PHI_t PHI, double *RPAR, int *IPAR, double *PAST,int *IPAST, int *NRDS)
{
	int i, len;
	PyObject *seq;
	PyObject *ret;
	PyObject *y_array;
	npy_intp dims[1] = {*N};
	params_t *p = (params_t*)IPAR;
	_lag_ctx ctx = {X,Y,ARGLAG,PHI,RPAR,IPAR,PAST,IPAST,NRDS,{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0}};
	
	/* If it's a C function then it's easy, just run it. */
	if(p->c_fcn)
	{
		p->c_fcn(F,Y,*X,(void*)_lag_c,&ctx,RPAR);
		return;
	}
	
	if(p->fail)
		return;
	
	/* If it's a Python function, then we have to tell Python we're going to maybe modify things. */
	if(!p->interactive)
	{
		PyEval_RestoreThread(p->_thread_save);
		p->_thread_save = NULL;
	}
	
	p->lag_ctx = &ctx;
	
	/* Create an array from Y values */
	y_array = PyArray_SimpleNewFromData(1,dims,NPY_DOUBLE,Y);
		
	/* Call the function, with following signature : func(Y, X, lagfun, params if any) */
	if(p->py_rpar)
		ret = PyObject_CallFunctionObjArgs(p->py_fcn, y_array, PyFloat_FromDouble(*X), p->py_lag_callback, p->py_rpar, NULL);
	else
		ret = PyObject_CallFunctionObjArgs(p->py_fcn, y_array, PyFloat_FromDouble(*X), p->py_lag_callback, NULL);
	
	if(!ret)
	{
		PyErr_Print();
		Py_DECREF(y_array);
		if(!p->interactive) p->_thread_save = PyEval_SaveThread();
		p->fail = 1;
		return;
	}
	
	Py_DECREF(y_array);
	
	seq = PySequence_Fast(ret, "Equation must return a sequence !");
	if(!seq)
	{
		PyObject *a, *b, *c;
		double dy;
		
		PyErr_Fetch(&a, &b, &c);
		PyErr_Clear();
		
		dy = PyFloat_AsDouble(ret);
		
		if(!PyErr_Occurred() && *N == 1)
		{
			if(a) Py_DECREF(a);
			if(b) Py_DECREF(b);
			if(c) Py_DECREF(c);
			F[0] = dy;
			goto cleanup;
		}
		else
		{
			PyErr_Restore(a, b, c);
			if(!p->interactive) p->_thread_save = PyEval_SaveThread();
			p->fail = 4;
			return;
		}
	}
	
	len = PySequence_Size(ret);
	
	if(len != *N)
	{
		PySys_WriteStderr("Equation returned a sequence of wrong size (%d, expected %d)...\n",len,*N);
		PyErr_SetString(PyExc_RuntimeError, "Equation function returned bad F(Y) length.");
		p->fail = 5;
		if(!p->interactive) p->_thread_save = PyEval_SaveThread();
		Py_DECREF(ret);
		return;
	}
	
	/* Give the output values to the integrator */
	for(i = 0; i < *N; i++)
	{
		PyObject *item = PySequence_Fast_GET_ITEM(seq, i);
		F[i] = PyFloat_AsDouble(item);
	}
	
	Py_DECREF(seq); 
	
cleanup:	
	Py_DECREF(ret);
	
	if(!p->interactive) p->_thread_save = PyEval_SaveThread();
}

static double PHI(int *I, double *X, double *RPAR, int *IPAR)
{
	PyObject *ret;
	int i = *I - 1;
	params_t *p = (params_t*)IPAR;
	PyObject *fun = p->py_y0[i];
	double y0;
	
	if(fun)
	{	
		if(!p->interactive && p->_thread_save) {
			PyEval_RestoreThread(p->_thread_save);
		}
		
		ret = PyObject_CallFunction(fun, "d", *X);
		
		if(!ret || PyErr_Occurred())
		{
			PyErr_Print();
			p->fail = 1;
			if(!p->interactive && p->_thread_save) p->_thread_save = PyEval_SaveThread();
			return 0;
		}
		y0 = PyFloat_AsDouble(ret);
		if(PyErr_Occurred())
		{
			PyErr_Print();
			p->fail = 4;
		}
		Py_DECREF(ret);
		
		if(!p->interactive && p->_thread_save) p->_thread_save = PyEval_SaveThread();
		return y0;
	}
	else if(p->y0_spline)
	{
#ifdef WITH_GSL
		return gsl_spline_eval (p->y0_spline[i], *X + p->y0_t0, p->spline_acc);
#else
		double y = 0;
		cspline_eval(&p->y0_spline[i], *X + p->y0_t0, &y);
		return y;
#endif
//		fprintf(stderr,"PHI(%g->%g)[%d]=%g\n",*X,*X+p->y0_t0,i,y);
	}
	else
	{
		return p->constant_y0[i];
	}
}

static double ARGLAG(int *IL, double *X, double *Y,double *RPAR, int *IPAR, PHI_t PHI, double *PAST, int *IPAST, int *NRDS)
{
	PyObject *ret;
    int ilag = *IL-1;
    params_t *p = (params_t*)IPAR;
    PyObject *pyfun = p->py_lagfuns[ilag];
	lag_c_t cfun = p->c_lagfuns[ilag];
	double r;
	
    if(!pyfun)
    {
        if(!cfun)
        {
            return *X-p->constant_lags[ilag];
        }
        else
        {
            return cfun(*X,Y);
        }
    }
    else
    {
        if(!p->interactive && p->_thread_save) PyEval_RestoreThread(p->_thread_save);
        ret = PyObject_CallFunction(pyfun, "dd", *X, *Y);
        if(!ret || PyErr_Occurred())
        {
            PyErr_Print();
        }
        
		r = PyFloat_AsDouble(ret);
        Py_DECREF(ret);
//          fprintf(stderr,"Lag %d: %p -> %f\n", ilag, fun, r);
	if(!p->interactive && p->_thread_save) p->_thread_save = PyEval_SaveThread();
        return r;
    }
}

static void JFCN(int *N, double *X, double *Y, double *DFY, int LDFY, ARGLAG_t ARGLAG, PHI_t PHI, double *RPAR, int *IPAR, double *PAST,int *IPAST, int *NRDS)
{
    fprintf(stderr,"JFCN\n");
}
static void JACLAG(int *N, double *X, double *Y, int *DFYL, ARGLAG_t ARGLAG, PHI_t PHI, int *IVE, int *IVC, int *IVL, double *RPAR, int *IPAR, double *PAST, int *IPAST, int *NRDS)
{
    fprintf(stderr,"JACLAG\n");
}


static void SOLOUT(int *NR, double *XOLD, double *X, double *HSOL, double *Y, double *CONT, int *LRC, int *N, double *RPAR, int *IPAR, int *IRTRN)
{   
    int i;
    params_t *p = (params_t*)IPAR;
	double *ptr;
	double *new_array;
	double next_x;
    
    if(p->interactive)
    {
        if(PyErr_CheckSignals())
        {
            p->fail = 6;
        }
    }
      
    if(p->fail)
    {
        *IRTRN = -1;
        return;
    } 
    
//   fprintf(stderr,"h=%g, dx=%g\n",*HSOL,*X-*XOLD);
//     fprintf(stderr,"SOLOUT %f %f / %p %d %d %d\n",*X,*Y, p->array,p->array_size,p->array_pos,p->array_incr);
    
    if(p->xvalues)
    {    
	    if(*NR == 1)
	    {
			ptr = p->array + p->array_pos * (*N+1);
			ptr[0] = *X;
		
			for(i = 0; i < *N; i++)
			{
				ptr[i + 1] = Y[i]; 
			}
			
			p->array_pos += 1;
			p->xpos += 1;
	    }
	    
	    next_x = p->xvalues[p->xpos];
	    
	    while(*X >= next_x && p->xpos < p->xvalues_len)
	    {
			if(p->array_size == p->array_pos)
			{
				p->array_size += p->array_incr;
			//         fprintf(stderr,"realloc %d->%d (%i)", p->array_size - p->array_incr, p->array_size , p->array_size*(*N+1)*sizeof(double));
				new_array = PyMem_Realloc(p->array, p->array_size*(*N+1)*sizeof(double));
				if(new_array == NULL)
				{
					p->fail = 2;
					*IRTRN = -1;
					return;
				}
			//         fprintf(stderr," -> %p\n", new_array);
				p->array = new_array;
			}    
			
			ptr = p->array + p->array_pos * (*N+1);
			ptr[0] = next_x;
		
			for(i = 1; i < *N + 1; i++)
			{
				ptr[i] = contr5_(&i, N, &next_x, CONT, X, HSOL);
			}
			
			p->array_pos += 1;
			p->xpos += 1;
			next_x = p->xvalues[p->xpos];
	    }
    }
    
        
    if(p->full_output)
    {
		if(p->array_size == p->array_pos)
		{
			p->array_size += p->array_incr;
		//         fprintf(stderr,"realloc %d->%d (%i)", p->array_size - p->array_incr, p->array_size , p->array_size*(*N+1)*sizeof(double));
			new_array = PyMem_Realloc(p->array, p->array_size*(*N+1)*sizeof(double));
			if(new_array == NULL)
			{
				p->fail = 2;
				*IRTRN = -1;
				return;
			}
		//         fprintf(stderr," -> %p\n", new_array);
			p->array = new_array;
		}	    
			
		ptr = p->array + p->array_pos * (*N+1);
		ptr[0] = *X;
		
		for(i = 0; i < *N; i++)
		{
			ptr[i+1] = Y[i]; 
		}
		
		p->array_pos += 1;
    }
}

void _tcc_error(void *opaque, const char *msg)
{
	PyErr_SetString(PyExc_RuntimeError, msg);
	//PySys_WriteStderr("! %s\n", msg);
}


// radar5(equ,y0,xend,[rpar],grid=grid)  

static params_t *current_p = NULL;

static PyMethodDef py_lag_def = {"lag", _lag_py, METH_VARARGS, "Returns delayed values of variables inside right-hand-side function."};

static PyObject *radar5_radar5(PyObject *self, PyObject *args,PyObject *keywds)
{    
	PyObject *output_array = NULL;
	int verbose = 0;
	int NGRID = 200;
	double XEND = 100;
	int N = 0;
        int MXST = 4*16384;
	params_t par;
	params_t *p = &par;
	int *IPAR = (void*)p;
	
	double X = 0;
	double *Y = NULL;

	double H = 1e-6;
	double RTOL = 1e-6;
	double ATOL = RTOL;
	double ITOL = 0;
	int IJAC = 0;
	int MLJAC = N;
	int MUJAC = 0; // not needed
	int NLAGS = 0;
	int NJACL = 0;
	int IMAS = 0;
	int IOUT = verbose ? 1 : 3;
	double WORK[30];
	int IWORK[30];
	double *RPAR = NULL;

	int IDID = 0;

	int *IPAST = NULL;
	int DUMMY;
	int MLMAS = N;
	int MUMAS = 0; // not needed
	
	int RPAR_len = 0;
	
	double *GRID = NULL;

	int i;
	
	int prev_interactive;
	
	PyObject* rpar_obj=NULL;
	PyObject* lagvars_obj=NULL;
	PyObject* lagfuns_obj=NULL; 
	PyObject* equ_obj=NULL;
	PyObject* y0_obj=NULL;
	PyObject* xend_obj=NULL;
	
	
	int full_output_flag = 0;
	
	PyObject *py_p = NULL;
	
	int len;

#ifdef WITH_TCC
	TCCState *tcc = NULL;
#endif	
	
	/* Parse arguments */
	static char *kwlist[] = {"equ", "y0", "xend", "rpar", "ngrid", "lagvars", "lags", "verbose", "mxst", "full_output", "rtol", "atol", "initial_stepsize", NULL};
	
	memset(p,0,sizeof(params_t));
	
	if (!PyArg_ParseTupleAndKeywords(args, keywds, "OOO|OiOOiiifff", kwlist, &equ_obj, &y0_obj, &xend_obj, &rpar_obj, &NGRID, &lagvars_obj, &lagfuns_obj, &verbose, &MXST,&full_output_flag,&RTOL,&ATOL,&H))
		return NULL;
        
	if(verbose) PySys_WriteStderr("RADAR5 wrapper (verbose = %d)\n", verbose);

	
	/* SPECIFIED X VALUES OR XEND */
	if(xend_obj)
	{
		PyArrayObject *xvalues = (PyArrayObject*)PyArray_ContiguousFromAny(xend_obj, NPY_DOUBLE, 0, 1);
		if(xvalues == NULL)
		{
			PyErr_SetString(PyExc_ValueError, "Cannot parse x values.");
			goto cleanup;
		}
		
		p->full_output = full_output_flag;	/* Full output only if specified... */
		
		if(PyArray_NDIM(xvalues) == 0)	/* if a single value if given for X, force full output */
		{
			p->xvalues = NULL;
			p->full_output = 1;
			XEND = *(double*)PyArray_GETPTR1(xvalues, 0);
			if(verbose) PySys_WriteStderr("Full output up to X=%g (all points)\n", XEND);
		}
		else
		{
			p->xvalues = (double*)PyArray_GETPTR1(xvalues, 0);
			p->xvalues_len = PyArray_DIM(xvalues, 0);
			p->xpos = 0;
			XEND = p->xvalues[p->xvalues_len - 1];
			if(verbose) PySys_WriteStderr("%s output up to X=%g (%d points)\n", p->full_output ? "Full" : "User-specified", XEND, p->xvalues_len);
		}
	}
	
	/* LOADING INITIAL VALUES : */
	
	if(y0_obj)
	{
		PyObject *seq = PySequence_Fast(y0_obj, "Initial values must be provided as a sequence.");
		if(!seq)
			return NULL;
		len = PySequence_Size(y0_obj);

		if(len == 2)
		{
			/* It might be a tuple (time, values) */
			PyObject *time = PySequence_Fast_GET_ITEM(seq, 0);
			PyObject *values = PySequence_Fast_GET_ITEM(seq, 1);
			
			/* Check each of (time, values) is iterables. If it is, we suppose we can use them as array. */
			if(PySequence_Check(time) && PySequence_Check(values))
			{
				int y0_len;
				
				PyArrayObject *time_array = (PyArrayObject*)PyArray_ContiguousFromAny(time, NPY_DOUBLE, 1, 1);
				PyArrayObject *values_array = (PyArrayObject*)PyArray_ContiguousFromAny(values, NPY_DOUBLE, 1, 2);
				
				/* Account for format errors */
				if(!time_array || !values_array)
				{
					PyErr_Print();
					return NULL;
				}
				
				y0_len = PyArray_DIM(time_array, 0);	/* Get initial sequence length */
				
				if(PyArray_NDIM(values_array) == 1) N = 1;
				else N = PyArray_DIM(values_array, 1);	/* Get system dimension from initial data */
				
				p->y0_t0 = *(double*)PyArray_GETPTR1(time_array, y0_len - 1);	/* Get the final time of the initial sequence */
				if(verbose) PySys_WriteStderr("Initial values array as (time->%g, [i0:%d])\n", p->y0_t0, y0_len);
#ifdef WITH_GSL
				p->y0_spline = (gsl_spline**)malloc(sizeof(void*) * N);
#else
				p->y0_spline = (cspline_t*)malloc(sizeof(cspline_t) * N);
#endif				
				
				for(i = 0; i < N; i++)
				{
					int j;
#ifdef WITH_GSL
					gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, y0_len);
#else
					cspline_t *spline = &p->y0_spline[i];
					cspline_init(spline);
					cspline_alloc_xy(spline, y0_len);
#endif

					/* Copy time data */
					memcpy (spline->x, (double*)PyArray_DATA(time_array), y0_len * sizeof(double));
					
					/* Copy variable data */
					for(j = 0; j < y0_len; j++)
					{
						spline->y[j] = *(double*)PyArray_GETPTR2(values_array, j, i);
					}
					
					/* Initialize interpolator */
#ifdef WITH_GSL					
					gsl_interp_init (spline->interp, spline->x, spline->y, y0_len);
					p->y0_spline[i] = spline;
#else
					cspline_compute_nat (spline);
#endif
				}
				
				Py_DECREF(time_array);
				Py_DECREF(values_array);
			}
		}
		else if(len > 0 && PySequence_Check(PySequence_Fast_GET_ITEM(seq, 0)))
		{
			int y0_len;
			
			/* Looks like a 2D array : time data might be its first row. This is useful to "continue" an integration, because it's the integrator output format. */
			PyArrayObject *full_array = (PyArrayObject*)PyArray_ContiguousFromAny(y0_obj, NPY_DOUBLE, 2, 2);
			if(PyErr_Occurred())
			{
				Py_DECREF(seq);
				goto cleanup;
			}
			
			y0_len = PyArray_DIM(full_array, 0);
			N = PyArray_DIM(full_array, 1) - 1;
					
			p->y0_t0 = *(double*)PyArray_GETPTR2(full_array, y0_len - 1, 0);
#ifdef WITH_GSL
			p->y0_spline = (gsl_spline**)malloc(sizeof(void*) * N);
#else
			p->y0_spline = (cspline_t*)malloc(sizeof(cspline_t) * N);
#endif	
			
			if(verbose) PySys_WriteStderr("Initial values array as [time->%g, y0:%d]\n", p->y0_t0, y0_len);
				
			for(i = 0; i < N; i++)
			{
				int j;
#ifdef WITH_GSL
				gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, y0_len);
#else
				cspline_t *spline = &p->y0_spline[i];
				cspline_init(spline);
				cspline_alloc_xy(spline, y0_len);
#endif

				for(j = 0; j < y0_len; j++)
				{
					spline->x[j] = *(double*)PyArray_GETPTR2(full_array, j, 0);
					spline->y[j] = *(double*)PyArray_GETPTR2(full_array, j, i + 1);
				}
					
#ifdef WITH_GSL					
				gsl_interp_init (spline->interp, spline->x, spline->y, y0_len);
				p->y0_spline[i] = spline;
#else
				cspline_compute_nat (spline);
#endif
			}
			
			Py_DECREF(full_array);
		}
		
		if(!p->y0_spline)
		{
			/* Each variable is specified on its own. */
			for (i = 0; i < len; i++)
			{
				PyObject *item = PySequence_Fast_GET_ITEM(seq, i);
				
				/* Check whether it's a Python function or a constant value. */
				if(PyCallable_Check(item))
				{
					p->py_y0[i] = item;
					if(verbose)
					{
						PyObject *repr = PyObject_Repr(item);
#if PY_MAJOR_VERSION >= 3
						PySys_WriteStderr("Initial value: y[%d](t<0) = %s\n", i, PyUnicode_AsUTF8(repr));
#else			
						PySys_WriteStderr("Initial value: y[%d](t<0) = %s\n", i, PyString_AsString(repr));
#endif
						Py_DECREF(repr);
					}
				}
				else
				{
					double y0 = PyFloat_AsDouble(item);
					if(PyErr_Occurred())
					{
						Py_DECREF(seq);
						goto cleanup;
					}
					p->constant_y0[i] = y0;
					p->py_y0[i] = NULL;
					if(verbose) PySys_WriteStderr("Initial value: y[%d](t<0) = %g (constant)\n", i, y0);
				}
			}
			
			N = len;
		}
		
		Py_DECREF(seq);
	}
	
	if(N == 0)
	{
		PyErr_SetString(PyExc_ValueError, "Could not guess system dimension from initial values.");
		goto cleanup;
	}
	
	if(verbose) PySys_WriteStderr("System dimension: N = %d\n", N);
	
	Y = (double*)malloc(N * sizeof(double));
	IPAST = (int*)malloc(N * sizeof(int));

	/* Set all integrator parameters to default values. */
	for(i = 0; i < 30; i++)
	{
		WORK[i] = 0;
		IWORK[i] = 0;
	}
	
	/* Load system numerical parameters. */
	if(rpar_obj)
	{
		PyObject *seq = PySequence_Fast(rpar_obj, "Constant parameters must be provided as a sequence.");
		if(!seq) goto cleanup;
		RPAR_len = PySequence_Size(rpar_obj);
		
		RPAR = malloc(sizeof(double) * RPAR_len);
		
		for (i = 0; i < RPAR_len; i++) {
			PyObject *item = PySequence_Fast_GET_ITEM(seq, i);
			RPAR[i] = PyFloat_AsDouble(item);
			
			if(verbose) PySys_WriteStderr("Constant parameter %d : %g\n", i,RPAR[i]);
		}
		Py_DECREF(seq); 
	}
        
	if(lagvars_obj)
	{
		PyObject *seq = PySequence_Fast(lagvars_obj, "Delayed variables must be provided as a sequence.");
		if(!seq) goto cleanup;
		len = PySequence_Size(lagvars_obj);

		for (i = 0; i < len; i++) {
			PyObject *item = PySequence_Fast_GET_ITEM(seq, i);
			int lagparam = PyLong_AsLong(item);
			
			if(lagparam < 0 || lagparam >= N)
			{
				PyErr_SetString(PyExc_ValueError, "Invalid lag parameter");
				Py_DECREF(seq);
				goto cleanup;
			}
			
			IPAST[i]=lagparam+1;
			
			if(verbose) PySys_WriteStderr("Variable %d is delayed (#%d)\n", lagparam, i);
		}
		
		IWORK[14] = i; // number of dense solution (lags)
		Py_DECREF(seq); 
	}
	else
	{
		/* By default, all variables are delayed */
		for (i = 0; i < N; i++)
		{
			IPAST[i]=i+1;
		}
		IWORK[14] = N;
		if(verbose) PySys_WriteStderr("Delayed variables unspecified: by default all %d variables are selected.\n", N);
	}
	
			
	if(lagfuns_obj)
	{
		PyObject *seq = PySequence_Fast(lagfuns_obj, "Lags must be provided as a sequence.");
		if(!seq) goto cleanup;
		len = PySequence_Size(lagfuns_obj);

		for (i = 0; i < len; i++) {
			PyObject *item = PySequence_Fast_GET_ITEM(seq, i);
			
			if(PyCallable_Check(item))
			{
				p->py_lagfuns[i]=item;
				p->c_lagfuns[i]=NULL;
				
				if(verbose)
				{
					PyObject *repr = PyObject_Repr(item);
#if PY_MAJOR_VERSION >= 3
					PySys_WriteStderr("Lag %d is a Python function: %s\n", i, PyUnicode_AsUTF8(repr));
#else			
					PySys_WriteStderr("Lag %d is a Python function: %s\n", i, PyString_AsString(repr));
#endif
					Py_DECREF(repr);
				}
			}
			else if(PyNumber_Check(item))
			{
				double delay = PyFloat_AsDouble(item);
				p->constant_lags[i]=delay;
				p->py_lagfuns[i]=NULL;
				p->c_lagfuns[i]=NULL;
				if(verbose) PySys_WriteStderr("Lag %d is constant: %g\n", i, delay);
			}
			else
			{
				p->c_lagfuns[i]=PyCapsule_GetPointer(item, NULL);
				p->py_lagfuns[i]=NULL;
				if(verbose) PySys_WriteStderr("Lag %d is a raw C function: %p\n", i, p->c_lagfuns[i]);
			}
		}
		
		Py_DECREF(seq); 
	}
	
	/* Set RHS function, either Python or C. */
	if(PyCallable_Check(equ_obj))
	{
		p->py_fcn = equ_obj;
		p->c_fcn = NULL;
		
		if(verbose)
		{
			PyObject *repr = PyObject_Repr(equ_obj);
#if PY_MAJOR_VERSION >= 3
			PySys_WriteStderr("Equation is a Python function:  %s\n", PyUnicode_AsUTF8(repr));
#else			
			PySys_WriteStderr("Equation is a Python function:  %s\n", PyString_AsString(repr));
#endif
			Py_DECREF(repr);
		}
		
		/* Create a Python wrapper for the function that will be used to get delayed values in RHS */
		py_p = PyCapsule_New(p, NULL, NULL);
		p->py_lag_callback = PyCFunction_New(&py_lag_def, py_p);
		
		/* We need a Python version of RPAR */
		if(RPAR_len > 0)
		{
			npy_intp dims[1] = {RPAR_len};
			p->py_rpar = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, RPAR);
		}
	}
#if PY_MAJOR_VERSION >= 3	
	else if(PyUnicode_Check(equ_obj))
	{
#ifdef WITH_TCC		
		PyObject *str = PyUnicode_AsASCIIString(equ_obj);
		if(!str)
			goto cleanup;
#endif
#else
	else if(PyString_Check(equ_obj))
	{
		PyObject *str = equ_obj;
#endif
#ifdef WITH_TCC

		/* We need to provide some math functions in order to do useful stuff. */
		const struct
		{
			const char *name;
			void *ptr;
		} symbols[] = {{"cos", cos},	{"sin", sin},	{"tan", tan},
					   {"acos", acos},	{"asin", asin},	{"atan", atan},
					   {"log", log},	{"pow", pow},	{"log10", log10},
					   {"cosh", cosh},	{"sinh", sinh},	{"tanh", tanh}};
		const int nsymbols = sizeof(symbols)/sizeof(symbols[0]);
		char *source = malloc(PyBytes_Size(str) + nsymbols * 128 + 128);
		
		/* This one can be useful too. */
		int pos = sprintf(source, "#define M_PI 3.14159265358979323846\n");
		
		for(i = 0; i < nsymbols; i++)
			pos += sprintf(source + pos, "double %s(double);", symbols[i].name);

		sprintf(source + pos, "\n%s", PyBytes_AsString(str));
		
		if(verbose) PySys_WriteStderr("Equation is C code, try to compile it...\n");
		
		/* Create a minimal TCC context */
		tcc = tcc_new();
		tcc_set_options(tcc, "-nostdlib");
		tcc_set_error_func(tcc, NULL, _tcc_error);
		tcc_set_output_type(tcc, TCC_OUTPUT_MEMORY);
		
		/* Compile the code snippet */
		if(tcc_compile_string(tcc, source) < 0)
		{
			free(source);
			goto cleanup;
		}
		
		//Py_DECREF(str);
		free(source);
		
		/* Add math symbols for the linking step. */
		for(i = 0; i < nsymbols; i++)
			tcc_add_symbol(tcc, symbols[i].name, symbols[i].ptr);
		
		/* Linking */
		if(tcc_relocate(tcc, TCC_RELOCATE_AUTO) < 0)
		{
// 			PyErr_SetString(PyExc_RuntimeError, "Equation could not be compiled (relocation error).");
			goto cleanup;
		}
		
		p->py_fcn = NULL;
		p->c_fcn = tcc_get_symbol(tcc, "equation");		/* Get pointer to the compiled function */
		
		if(!p->c_fcn)
		{
			p->c_fcn = tcc_get_symbol(tcc, "y");
			if(!p->c_fcn)
			{
				p->c_fcn = tcc_get_symbol(tcc, "rhs");
				if(!p->c_fcn)
				{
					PyErr_SetString(PyExc_ValueError, "No function named 'equation', 'y' or 'rhs' in provided C code.");
					goto cleanup;
				}
			}
		}
		
		if(verbose) PySys_WriteStderr("Equation successfully compiled : %p\n", p->c_fcn);
#else
		PyErr_SetString(PyExc_RuntimeError, "Equation provided as string, but we were compiled without TCC.");
		goto cleanup;
#endif
	}
	else
	{
		p->py_fcn = NULL;
		p->c_fcn = PyCapsule_GetPointer(equ_obj, NULL);
		if(verbose) PySys_WriteStderr("Equation is a raw C function: %p\n", p->c_fcn);
	}
	
	GRID = malloc(sizeof(double) * (NGRID + 1));
			
	
	if(p->xvalues)
	{
		p->array_size = p->xvalues_len;
	}
	else
	{
		p->array_size = NGRID;    // taille initiale
	}
	
	p->array = PyMem_Malloc(p->array_size*(N+1)*sizeof(double));
	if(verbose) PySys_WriteStderr("Allocate initial array (size %gkB) at %p\n", p->array_size*(N+1.)*sizeof(double)/1024, p->array);
	
	p->array_pos = 0;
	p->array_incr = p->array_size / 4;    // increment taille
	
	p->user_params = RPAR;
        
	p->fail = 0;
	p->interactive = 0;
	
	IWORK[13] = 1;      //  ne pas utiliser JACLAG
	
	/* NOTE: Black magic... */
	IWORK[12] = NGRID;      //  WORKSPACE FOR GRID
	
//	MXST = 1024*1024*32; //4*16384*IWORK[14]; //MAX(NGRID*IWORK[14],8192); // pourquoi c'était divisé par /N ?
	IWORK[11] = MXST;       // WORKSPACE FOR PAST 
	
	IWORK[1] = INT_MAX;
	
	for(i = 0; i < NGRID; i++)
	{
		GRID[i] = XEND/(NGRID+1)*i;
	} 
        
        /* Load initial values */
	prev_interactive = p->interactive;
	p->interactive = 1;
	for(i=1;i<N+1;i++)
	{
		Y[i-1] = PHI(&i, &X, RPAR, IPAR);
		if(verbose) PySys_WriteStderr("Initial value: Y[%d]=%g\n", i-1, Y[i-1]);
	}
	p->interactive = prev_interactive;
	
        
	/* MAIN CALL TO THE INTEGRATOR AHEAD */
	
	if(!p->interactive)
	{
		p->_thread_save = PyEval_SaveThread();	/* Promise Python that we are not going to mess with any of its data. */
							/* This allows other threads/event loops to run during the integration. */
	}

	current_p = p;
	
	radar5_(&N,FCN,PHI,ARGLAG,&X,Y,&XEND,&H, \
					&RTOL,&ATOL,&ITOL, \
					JFCN,&IJAC,&MLJAC,&MUJAC, \
					JACLAG,&NLAGS,&NJACL, \
					&IMAS,SOLOUT,&IOUT, \
					WORK,IWORK,RPAR,IPAR,&IDID, \
					GRID,IPAST,&DUMMY,&MLMAS,&MUMAS);

	current_p = NULL;
	
	if(!p->interactive)
	{
		PyEval_RestoreThread(p->_thread_save);
	}
	
	/* DONE ! */
	
	if(verbose)
	{
		switch(p->fail)
		{
			case 1: PySys_WriteStderr("Python function produced an error.\n"); break;
			case 3: PySys_WriteStderr("Computation halted by stop().\n"); break;
			case 4: 
			case 5: PySys_WriteStderr("Python function provided wrong return type.\n"); break;
			case 6: PySys_WriteStderr("Computation halted by an external signal.\n"); break;
		}
	}
	
	if(p->fail == 2)
	{
		PyErr_SetString(PyExc_MemoryError, "Cannot grow output array.");
	}
	
	if(verbose)
	{
		PySys_WriteStderr("Solver output: IDID = %d (", IDID);
		switch(IDID)
		{
			case  1: PySys_WriteStderr("Successful computation"); break;
			case  2: PySys_WriteStderr("Successful computation interrupted by routine SOLOUT"); break;
			case -1: PySys_WriteStderr("Non-consistent input values"); break;
			case -2: PySys_WriteStderr("Too many stepsizes required (>NMAX)"); break;
			case -3: PySys_WriteStderr("Stepsize becomes too small"); break;
			case -4: PySys_WriteStderr("Jacobian matrix repeatedly singular"); break;
			case -5: PySys_WriteStderr("Computation interrupted by routine YLAGR5"); break;
			case -6: PySys_WriteStderr("The equation makes use of advanced arguments"); break;
			default: PySys_WriteStderr("Unexpected return code"); break;
		}
		PySys_WriteStderr(")\nComputed points: %d\n", p->array_pos);
	}
	
	/* Create numpy array based on integrator output : */
	{
		npy_intp dims[2] = {p->array_pos, N+1};
		output_array = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, p->array);
	}
        
#if NPY_API_VERSION >= 9
	PyArray_ENABLEFLAGS((PyArrayObject*)output_array, NPY_ARRAY_OWNDATA);
#else
	((PyArrayObject*)output_array)->flags |= NPY_OWNDATA;
#endif
 
	
	/* CLEANUP : */
cleanup:

	if(RPAR) free(RPAR);
	if(Y) free(Y);
	if(IPAST) free(IPAST);
	if(GRID) free(GRID);
	
	if(p->y0_spline)
	{
		if(verbose) PySys_WriteStderr("Clearing initial values...\n");
		for(i = 0; i < N; i++)
		{
#ifdef WITH_GSL			
			if(p->y0_spline[i])
				gsl_spline_free(p->y0_spline[i]);
#else
			cspline_free(&p->y0_spline[i]);
#endif
		}
		free(p->y0_spline);
	}
	
#ifdef WITH_GSL	
	if(p->spline_acc)
	{
		gsl_interp_accel_free (p->spline_acc);
	}
#endif
	
	if(p->py_lag_callback)
		Py_DECREF(p->py_lag_callback);
	
	if(py_p)
		Py_DECREF(py_p);
	
	if(p->py_rpar)
		Py_DECREF(p->py_rpar);
	
	if(p->xvalues)
		Py_DECREF(p->xvalues);
	
#ifdef WITH_TCC
	if(tcc)
		tcc_delete(tcc);
#endif		

	return output_array;
}

static PyObject *radar5_stop(PyObject *self, PyObject *args)
{    
	if(current_p)
	{
		current_p->fail = 3;
	}
	Py_RETURN_NONE;
}


PyDoc_STRVAR(radar5_doc,
"Call the RADAR5 solver, which numerically solves delayed differential equations (DDE) using an implicit collocation Runge-Kutta method based on Radau nodes.\n"
"\n"
"Beware : some functionalities are not implemented yet :\n"
" - providing custom jacobians ;\n"
" - mass matrix ;\n"
" - lots of options...\n"
"\n"
"For more extensive details, see the documentation [1]_ of the original FORTRAN code, which is the work of Nicola Guglielmi and Ernst Hairer.\n"
"More details can be found on : http://www.unige.ch/~hairer/software.html\n"
"\n"
":param equation : callable, pointer to c function, or string\n"
"    Right-hand side of the equation, computing the time derivatives.\n"
"    It is fundamentally a function that takes three or four arguments:\n"
"        * y[i], an array containing the current values:\n"
"        * t, the current time;\n" 
"        * lag(i,j), a function that can be called compute delayed variable i with the lag j.\n"
"        * params, an extra array of user-defined parameters\n\n"
"   This function can be a Python callable, but this can be quite slow when integrating large systems. For better performances, a pointer to a C function can be provided. This function can be loaded by ctypes, or be part of a C extension, etc. It must have the following signature:\n"
"     void equation(double *dy, const double *y, double t, double (*lag)(int,int,void*), void* ctx, void *params)\n"
"   The output derivatives should be written in dy, and the extra 'ctx' argument has to be provided as 'lag' third parameters.\n"
"   For more convenience, a string containing C code can also be provided. It must contain a function called 'equation', 'equ', or 'rhs' with the same signature. It will be internally compiled using the 'tcc' compiler. Standard math function (sin, cos, ...) are made available in the global scope.\n"
"\n"
":param y0 : array, or tuple of arrays\n"
"    Provides initial and past values. It can be either a constant, and in that case all t<0 values are the same. If non-constant past values are needed, they must be provided as a tuple (time, values). Missing points will be computed using a cubic spline interpolation. For a more precise control, a list of Python callable computing the past values for each variable can also be used. Finally, a single array can be used, and its first row will be considered as the time values. This is particularly useful for resuming integration from the output of another call.\n"
"\n"
":param time : array, or float\n"
"    Specifies the points on which the result of the integration has to be stored. Note that this does not correspond to the integrating grid, as the solver uses an adaptative step. If a single number is provided, full output is turned on, and all output steps are stored.\n"
"\n"
":param params : array, optional\n"
"    Extra parameters to be provided to the integrating function.\n"
"\n"
":param lagvars : array, optional\n"
"    List of variables subjected to delay. If it is not specified, all variables are considered delayed variable. You can specifiy an empty list or tuple if no delayed variable are needed.\n"
"\n"
":param lags : array, optional\n"
"    List of delays used in the equations. They can be either constant, or functions of the other variables or current time.\n"
"\n"
":param ngrid : int, optional\n"
"    Changes the size of the internal grid. Might need some tuning if integration fails unexpectedly."
":param mxst : int, optional\n"
"    Changes the size of the internal past storage. Might need some tuning if integration fails unexpectedly with 'MXST' related errors."
":param atol : float, optional\n"
"    Changes the absolute tolerance (default is 1e-6)\n"
":param rtol : float, optional\n"
"    Changes the relative tolerance (default is 1e-6)\n"
":param initial_stepsize : float, optional\n"
"    Initial guess for the stepsize; for stiff equations with initial transient, H=1/(norm of F'); usually 1e-3 or 1e-5 is good. (default is 1e-6)\n" 
":param verbose : bool, optional\n"
"    Enable printing of troubleshooting informations."
":param full_output : bool, optional\n"
"    Turns on the output of every integration step. This can make the output array much bigger than the specified times.\n"
"\n"
":return: The result of the integration, the first row of which being the \"time\" variable, and the latter the integrated components.\n"
"\n"
".. References :\n"
"[1] N. Guglielmi and E. Hairer, \"Users' Guide for the code RADAR5 - Version 2.1\", Technical Report, July 2005.\n"
);

static PyMethodDef Radar5Methods[] = {
    {"radar5",  (PyCFunction)radar5_radar5, METH_VARARGS|METH_KEYWORDS, radar5_doc},
    {"stop",  (PyCFunction)radar5_stop, METH_VARARGS, "Stop running computation."},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef Radar5Def = {
        PyModuleDef_HEAD_INIT,
        "radar5",
        NULL,
        -1,
		Radar5Methods
};

#define INITERROR return NULL

#else
#define INITERROR return
#endif

#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC PyInit_radar5(void)
{
    PyObject *module = PyModule_Create(&Radar5Def);
#else
PyMODINIT_FUNC initradar5(void)
{	
    PyObject *module = Py_InitModule("radar5", Radar5Methods);
#endif

    if (module == NULL)
        INITERROR;
	
	import_array();

	// TODO : add dummy integration WITH lag so integration without lag works properly...
	
#if PY_MAJOR_VERSION >= 3
    return module;
#endif
}

