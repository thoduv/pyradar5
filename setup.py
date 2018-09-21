# -*- coding: utf-8 -*-
try:
    import numpy
except ImportError:
    try:
       import pip 
       pip.main(['install', 'numpy'])
    except AttributeError:
       import pip._internal
       pip._internal.main(['install', 'numpy'])   
    import numpy

import platform
import setuptools
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration
#from setuptools.distutils.sysconfig import get_python_inc, get_python_lib

def configuration(parent_package='', top_path=None):	
		
	config = Configuration('radar5', parent_package, top_path)

	config.add_library('radar5f', sources=['radar5/radar5.f', 'radar5/dc_decdel.f','radar5/decsol.f' ,'radar5/contr5.f'],
						extra_f77_compile_args=['-static', '-w'],)

	config.add_extension('radar5',
						sources = ['wrapper.c','tcc/libtcc.c'],
						include_dirs=[numpy.get_include()],
						libraries = ['radar5f'])
						
	return config

metadata = dict(version= '0.6',
		    author= u'Aurélien Thorette',
		    author_email='thoduv@free.fr',
		    url='https://github.com/thoduv/pyradar5',
		    license='LGPL',
		    description = 'Python wrapper of the RADAR5 solver for delayed differential equations (DDE).',
		    long_description= '''The radar5 module is a numerical solver for delay differential equations. It is a wrapper around an original FORTRAN code, that uses a colocation method on Radau nodes for an efficient integration of stiff problems. Please note that all the capabilities of the FORTRAN code are not wrapped yet.
	The original FORTRAN code is the work of Nicola Guglielmi and Ernst Hairer, and can be found on: http://www.unige.ch/~hairer/software.html''',
			install_requires=['numpy'],
			setup_requires=['numpy'],
		    classifiers=[
	"Topic :: Scientific/Engineering",
	"Programming Language :: Python :: 3",
	"Programming Language :: Python :: 2",
	"Programming Language :: Python :: Implementation :: CPython",
        "License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX"
    ],
	configuration = configuration)

setup(**metadata)
