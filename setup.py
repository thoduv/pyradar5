# -*- coding: utf-8 -*-
from distutils.util import get_platform
from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration
from os import chdir, getcwd, system

with open("README.md", "r") as fh:
    long_description = fh.read()

    
config = Configuration('radar5', '', None)

prev = getcwd()
if 'win' in get_platform():
	chdir('tcc\\win32\\')
	system('build.bat')
else:
	chdir('tcc')
	system('./configure')
chdir(prev)

config.add_library('radar5f', sources=['radar5/radar5.f', 'radar5/dc_decdel.f','radar5/decsol.f' ,'radar5/contr5.f'])

config.add_extension('radar5',
                    sources = ['wrapper.c',
			       'tcc/libtcc.c'],
                    extra_compile_args = ['-Ofast'],
                    libraries = ['m','gsl', 'blas', 'lapack', 'radar5f'])

config.dict_append(version= '1.0',
		    author= u'Aurélien Thorette',
		    author_email='thoduv@free.fr',
		    description = 'RADAR5 DDE Solver Wrapper',
		    long_description=long_description,
		    long_description_content_type="text/markdown",
		    classifiers=[
	"Topic :: Scientific/Engineering",
	"Programming Language :: Python :: 3",
	"Programming Language :: Python :: 2",
	"Programming Language :: Python :: Implementation :: CPython",
        "License :: OSI Approved :: LGPL License",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX"
    ])

setup(**config.todict())
