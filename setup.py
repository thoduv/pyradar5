# -*- coding: utf-8 -*-
from distutils.util import get_platform
import numpy
from numpy.distutils.core import setup
from distutils.sysconfig import get_python_inc, get_python_lib
from numpy.distutils.misc_util import Configuration
from os import chdir, getcwd, system

def configuration(parent_package='', top_path=None):	
		
	config = Configuration(None, parent_package='', top_path=None)

	config.add_library('radar5f', sources=['radar5/radar5.f', 'radar5/dc_decdel.f','radar5/decsol.f' ,'radar5/contr5.f'],
						extra_f77_compile_args=['-static'],)

	config.add_extension('radar5',
						sources = ['wrapper.c','tcc/libtcc.c'],
						include_dirs=[numpy.get_include(), get_python_inc(), 'C:\\ProgramData\\Anaconda3\\Library\\include'],
						library_dirs=[get_python_lib(), 'C:\\ProgramData\\Anaconda3\\Library\\lib'],
	#                    libraries = ['gfortran', 'gcc', 'gsl', 'gslcblas', 'radar5f', 'winpthread','pthread','mingw32','quadmath','mingwex','win32'])
						libraries = ['radar5f', 'gsl'] + ['blas'] if not 'win' in get_platform() else [])
						
	return config

metadata = dict(version= '0.1',
		    author= u'Aurélien Thorette',
		    author_email='thoduv@free.fr',
		    description = 'Python wrapper of the RADAR5 solver for delayed differential equations (DDE).',
		    long_description= '''The radar5 module is a numerical solver for delay differential equations. It is a wrapper around an original FORTRAN code, that uses a colocation method on Radau nodes for an efficient integration of stiff problems. Please note that all the capabilities of the FORTRAN code are not wrapped yet.
	The original FORTRAN code is the work of Nicola Guglielmi and Ernst Hairer, and can be found on: http://www.unige.ch/~hairer/software.html''',
			install_requires=['numpy', 'gsl'],
		    classifiers=[
	"Topic :: Scientific/Engineering",
	"Programming Language :: Python :: 3",
	"Programming Language :: Python :: 2",
	"Programming Language :: Python :: Implementation :: CPython",
        "License :: OSI Approved :: LGPL License",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX"
    ],
	configuration = configuration)

setup(**metadata)
