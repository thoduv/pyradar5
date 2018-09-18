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

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

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
	PyThreadState *_thead_save;
	int interactive;
	int fail;
	void *user_params;
	
	PyObject *py_rpar;
	PyObject *py_lag_callback;

	gsl_interp_accel *spline_acc;
	gsl_spline **y0_spline;
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
	PyEval_RestoreThread(p->_thead_save);
	
	p->lag_ctx = &ctx;
	
	/* Create an array from Y values */
	npy_intp dims[1] = {*N};
	PyObject *y_array = PyArray_SimpleNewFromData(1,dims,NPY_DOUBLE,Y);
		
	/* Call the function, with following signature : func(Y, X, lagfun, params if any) */
	PyObject *ret;
	if(p->py_rpar)
		ret = PyObject_CallFunctionObjArgs(p->py_fcn, y_array, PyFloat_FromDouble(*X), p->py_lag_callback, p->py_rpar, NULL);
	else
		ret = PyObject_CallFunctionObjArgs(p->py_fcn, y_array, PyFloat_FromDouble(*X), p->py_lag_callback, NULL);
	
	if(!ret)
	{
		PyErr_Print();
		Py_DECREF(y_array);
		p->_thead_save = PyEval_SaveThread();
		p->fail = 1;
		return;
	}
	
	Py_DECREF(y_array);
	
	PyObject *seq = PySequence_Fast(ret, "Equation must return a sequence !");
	if(!seq)
	{
		p->_thead_save = PyEval_SaveThread();
		p->fail = 1;
		return;
	}
	
	int len = PySequence_Size(ret);
	
	if(len != *N)
	{
		PySys_WriteStderr("Equation returned a sequence of wrong size (%d, expected %d)...\n",len,*N);
		PyErr_SetString(PyExc_RuntimeError, "Equation function returned bad F(Y) length.");
		p->fail = 1;
		p->_thead_save = PyEval_SaveThread();
		Py_DECREF(ret);
		return;
	}
	
	/* Give the output values to the integrator */
	int i;
	for(i = 0; i < *N; i++)
	{
		PyObject *item = PySequence_Fast_GET_ITEM(seq, i);
		F[i] = PyFloat_AsDouble(item);
	}
	
	Py_DECREF(ret);
	Py_DECREF(seq); 
	
	p->_thead_save = PyEval_SaveThread();
}

static double PHI(int *I, double *X, double *RPAR, int *IPAR)
{
	int i = *I - 1;
	params_t *p = (params_t*)IPAR;
	PyObject *fun = p->py_y0[i];
	
	if(fun)
	{
		PyEval_RestoreThread(p->_thead_save);
		PyObject *ret = PyObject_CallFunction(fun, "d", *X);
		p->_thead_save = PyEval_SaveThread();
		if(!ret || PyErr_Occurred())
		{
			PyErr_Print();
		}
		double y0 = PyFloat_AsDouble(ret);
		Py_DECREF(ret);
//         fprintf(stderr,"Phi %d: %p -> %f\n", i, fun, y0);
		return y0;
	}
	else if(p->y0_spline)
	{
		return gsl_spline_eval (p->y0_spline[i], *X + p->y0_t0, p->spline_acc);
//		fprintf(stderr,"PHI(%g->%g)[%d]=%g\n",*X,*X+p->y0_t0,i,y);
	}
	else
	{
		return p->constant_y0[i];
	}
}

static double ARGLAG(int *IL, double *X, double *Y,double *RPAR, int *IPAR, PHI_t PHI, double *PAST, int *IPAST, int *NRDS)
{
    int ilag = *IL-1;
    params_t *p = (params_t*)IPAR;
    PyObject *pyfun = p->py_lagfuns[ilag];
    if(!pyfun)
    {
        lag_c_t cfun = p->c_lagfuns[ilag];
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
        PyEval_RestoreThread(p->_thead_save);
        PyObject *ret = PyObject_CallFunction(pyfun, "dd", *X, *Y);
        p->_thead_save = PyEval_SaveThread();
        if(!ret || PyErr_Occurred())
        {
            PyErr_Print();
        }
        double r = PyFloat_AsDouble(ret);
        Py_DECREF(ret);
//          fprintf(stderr,"Lag %d: %p -> %f\n", ilag, fun, r);
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
    
    if(p->interactive)
    {
        if(PyErr_CheckSignals())
        {
            p->fail = 1;
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
		double *ptr = p->array + p->array_pos * (*N+1);
		ptr[0] = *X;
	
		for(i = 0; i < *N; i++)
		{
			ptr[i + 1] = Y[i]; 
		}
		
		p->array_pos += 1;
		p->xpos += 1;
	    }
	    
	    double next_x = p->xvalues[p->xpos];
	    
	    while(*X >= next_x && p->xpos < p->xvalues_len)
	    {
		if(p->array_size == p->array_pos)
		{
			p->array_size += p->array_incr;
		//         fprintf(stderr,"realloc %d->%d (%i)", p->array_size - p->array_incr, p->array_size , p->array_size*(*N+1)*sizeof(double));
			double *new_array = realloc(p->array, p->array_size*(*N+1)*sizeof(double));
			if(new_array == NULL)
			{
				p->fail = 2;
				*IRTRN = -1;
				return;
			}
		//         fprintf(stderr," -> %p\n", new_array);
			p->array = new_array;
		}    
		
		double *ptr = p->array + p->array_pos * (*N+1);
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
		double *new_array = realloc(p->array, p->array_size*(*N+1)*sizeof(double));
		if(new_array == NULL)
		{
			p->fail = 2;
			*IRTRN = -1;
			return;
		}
	//         fprintf(stderr," -> %p\n", new_array);
		p->array = new_array;
	}	    
	    
	double *ptr = p->array + p->array_pos * (*N+1);
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
	
	memset(p,0,sizeof(params_t));
	
	int i;
	
	PyObject* rpar_obj=NULL;
	PyObject* lagvars_obj=NULL;
	PyObject* lagfuns_obj=NULL; 
	PyObject* equ_obj=NULL;
	PyObject* y0_obj=NULL;
	PyObject* xend_obj=NULL;
	
	
	int full_output_flag = -1;
	
	PyObject *py_p = NULL;
	
	PyObject* seq;
	int len;

#ifdef WITH_TCC
	TCCState *tcc = NULL;
#endif	
	
	/* Parse arguments */
	static char *kwlist[] = {   "equ", "y0", "xend", "rpar", "ngrid", "lagvars", "lags", "verbose", "mxst", "full_output", NULL};
	
	if (!PyArg_ParseTupleAndKeywords(args, keywds, "OOO|OiOOiii", kwlist, &equ_obj, &y0_obj, &xend_obj, &rpar_obj, &NGRID, &lagvars_obj, &lagfuns_obj, &verbose, &MXST,&full_output_flag))
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
		
		p->full_output = full_output_flag == 1 ? 1 : 0;	/* Full output only if specified... */
		
		if(full_output_flag != 1 && PyArray_SIZE(xvalues) == 1)	/* ...or if a single value if given for X */
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
		seq = PySequence_Fast(y0_obj, "Initial values must be provided as a sequence.");
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
				// better use PyArray_AsCArray ?
				PyArrayObject *time_array = (PyArrayObject*)PyArray_ContiguousFromAny(time, NPY_DOUBLE, 1, 1);
				PyArrayObject *values_array = (PyArrayObject*)PyArray_ContiguousFromAny(values, NPY_DOUBLE, 2, 2);
				
				/* Account for format errors */
				if(!time_array || !values_array)
				{
					PyErr_Print();
					return NULL;
				}
				
				int y0_len = PyArray_DIM(time_array, 0);	/* Get initial sequence length */
				N = PyArray_DIM(values_array, 1);	/* Get system dimension from initial data */
				
				p->y0_t0 = *(double*)PyArray_GETPTR1(time_array, y0_len - 1);	/* Get the final time of the initial sequence */
				p->y0_spline = (gsl_spline**)malloc(sizeof(void*) * N);
				
				if(verbose) PySys_WriteStderr("Initial values array as (time->%g, [y0:%d])\n", p->y0_t0, y0_len);
				
				for(i = 0; i < N; i++)
				{
					int j;
					gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, y0_len);

					/* Copy time data */
					memcpy (spline->x, (double*)PyArray_DATA(time_array), y0_len * sizeof(double));
					
					/* Copy variable data */
					for(j = 0; j < y0_len; j++)
					{
						spline->y[j] = *(double*)PyArray_GETPTR2(values_array, j, i);
					}
					
					/* Initialize interpolator */
					gsl_interp_init (spline->interp, spline->x, spline->y, y0_len);
					p->y0_spline[i] = spline;
				}
				
				Py_DECREF(time_array);
				Py_DECREF(values_array);
			}
		}
		else if(len > 0 && PySequence_Check(PySequence_Fast_GET_ITEM(seq, 0)))
		{
			/* Looks like a 2D array : time data might be its first row. This is useful to "continue" an integration, because it's the integrator output format. */
			PyArrayObject *full_array = (PyArrayObject*)PyArray_ContiguousFromAny(y0_obj, NPY_DOUBLE, 2, 2);
			if(PyErr_Occurred())
			{
				Py_DECREF(seq);
				goto cleanup;
			}
			
			int y0_len = PyArray_DIM(full_array, 0);
			N = PyArray_DIM(full_array, 1) - 1;
					
			p->y0_t0 = *(double*)PyArray_GETPTR2(full_array, y0_len - 1, 0);
			p->y0_spline = (gsl_spline**)malloc(sizeof(void*) * N);
			
			if(verbose) PySys_WriteStderr("Initial values array as [time->%g, y0:%d]\n", p->y0_t0, y0_len);
				
			for(i = 0; i < N; i++)
			{
				int j;
				gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, y0_len);

				for(j = 0; j < y0_len; j++)
				{
					spline->x[j] = *(double*)PyArray_GETPTR2(full_array, j, 0);
					spline->y[j] = *(double*)PyArray_GETPTR2(full_array, j, i + 1);
				}
					
				gsl_interp_init (spline->interp, spline->x, spline->y, y0_len);
				p->y0_spline[i] = spline;
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
						PySys_WriteStderr("Initial value: y[%d](t<0) = ", i);
						PyObject_Print(item, stderr, 0);
						PySys_WriteStderr("\n");
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
					if(verbose) PySys_WriteStderr("Initial value: y[%d](t<0) = %f (constant)\n", i, y0);
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
		seq = PySequence_Fast(rpar_obj, "Constant parameters must be provided as a sequence.");
		if(!seq) goto cleanup;
		RPAR_len = PySequence_Size(rpar_obj);
		
		RPAR = malloc(sizeof(double) * RPAR_len);
		
		for (i = 0; i < RPAR_len; i++) {
			PyObject *item = PySequence_Fast_GET_ITEM(seq, i);
			RPAR[i] = PyFloat_AsDouble(item);
			
			if(verbose) PySys_WriteStderr("Constant parameter %d : %f\n", i,RPAR[i]);
		}
		Py_DECREF(seq); 
	}
        
	if(lagvars_obj)
	{
		seq = PySequence_Fast(lagvars_obj, "Delayed variables must be provided as a sequence.");
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
	
			
	if(lagfuns_obj)
	{
		seq = PySequence_Fast(lagfuns_obj, "Lags must be provided as a sequence.");
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
					PySys_WriteStderr("Lag %d is a Python function: ", i);
					PyObject_Print(item, stderr, 0);
					PySys_WriteStderr("\n");
				}
			}
			else if(PyNumber_Check(item))
			{
				double delay = PyFloat_AsDouble(item);
				p->constant_lags[i]=delay;
				p->py_lagfuns[i]=NULL;
				p->c_lagfuns[i]=NULL;
				if(verbose) PySys_WriteStderr("Lag %d is constant: %f\n", i, delay);
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
	
	
	GRID = malloc(sizeof(double) * (NGRID + 1));
			
	
	if(p->xvalues)
	{
		p->array_size = p->xvalues_len;
	}
	else
	{
		p->array_size = NGRID;    // taille initiale
	}
	
	double *initial_array = malloc(p->array_size*(N+1)*sizeof(double));
	if(verbose) PySys_WriteStderr("Allocate initial array (size %gkB) at %p\n", p->array_size*(N+1.)*sizeof(double)/1024, initial_array);
	
	p->array =  initial_array;
	p->array_pos = 0;
	p->array_incr = p->array_size / 4;    // increment taille
	
	p->user_params = RPAR;
	
	/* Set RHS function, either Python or C. */
	if(PyCallable_Check(equ_obj))
	{
		p->py_fcn = equ_obj;
		p->c_fcn = NULL;
		
		if(verbose)
		{
			PySys_WriteStderr("Equation is a Python function: ");
			PyObject_Print(equ_obj, stderr, 0);
			PySys_WriteStderr("\n");
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
		
		if(verbose) PySys_WriteStderr("Equation is C code, try to compile it...\n");
		
		/* Create a minimal TCC context */
		tcc = tcc_new();
		tcc_set_options(tcc, "-nostdlib");
		tcc_set_error_func(tcc, NULL, _tcc_error);
		
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
		
		/* Compile the code snippet */
		int result = tcc_compile_string(tcc, source);
		if(result < 0)
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
		result = tcc_relocate(tcc, TCC_RELOCATE_AUTO);
		if(result < 0)
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
	for(i=1;i<N+1;i++)
	{
		Y[i-1] = PHI(&i, &X, RPAR, IPAR);
		if(verbose) PySys_WriteStderr("Initial value Y[%d]=%g\n", i, Y[i-1]);
	}
        

	/* MAIN CALL TO THE INTEGRATOR */
	
	if(!p->interactive)
	{
		p->_thead_save = PyEval_SaveThread();	/* Promise Python that we are not going to mess with any of its data. */
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
		PyEval_RestoreThread(p->_thead_save);
	}
	
	/* DONE ! */
	
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
	npy_intp dims[2] = {p->array_pos, N+1};
	output_array = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, p->array);
        
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
			if(p->y0_spline[i])
				gsl_spline_free(p->y0_spline[i]);
		}
		free(p->y0_spline);
	}
		
	if(p->spline_acc)
	{
		gsl_interp_accel_free (p->spline_acc);
	}
	
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
		current_p->fail = 1;
	}
	Py_RETURN_NONE;
}


PyDoc_STRVAR(radar5_doc,
    "radar5(equ, y0, xend, rpar=[], ngrid=200, lagvars=[], lags=[], verbose=False)\n"
	"\n"
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
    "Parameters\n"
	"----------\n"
    "equ : function, pointer to a raw c function or string\n"
	"y0 : iterable or callable, initial and past values: either constant, or a (time, values) array, or a functions.\n"
	"xvalues/xend : real or iterable : specified points where the integration should take place\n"
	"rpar : iterable\n"
	"ngrid : int\n"
	"lagvars : iterable\n"
	"lags : iterable\n"
	"verbose : bool\n"
	"full_output : bool\n"
	"\n"
    "Returns\n"
	"-------\n"
	"result : The result of the integration, the first row of which being the \"time\" variable. Please note that the length of the array is mostly unpredictable.\n"
	"\n"
	"References\n"
	"----------\n"
	"[1] N. Guglielmi and E. Hairer, \"Users' Guide for the code RADAR5 - Version 2.1\", Technical Report, July 2005."
    "\n");

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

