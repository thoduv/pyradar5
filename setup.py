from distutils.core import setup, Extension

module1 = Extension('radar5',
                    sources = ['wrapper.c'],
                    extra_objects = ['radar5/radar5.so', 'radar5/dc_decdel.so','radar5/decsol.so' ,'radar5/contr5.so'] +
                    ['tcc/'+x for x in ['libtcc.o', 'tccpp.o', 'tccgen.o', 'tccelf.o', 'tccasm.o', 'tccrun.o', 'x86_64-gen.o', 'x86_64-link.o', 'i386-asm.o']],
                    libraries = ['gfortran','m','gsl', 'blas', 'lapack'])

setup (name = 'radar5',
       version = '1.0',
       description = 'RADAR5 DDE Solver Wrapper',
       ext_modules = [module1])
