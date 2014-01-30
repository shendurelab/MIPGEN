import sys
# ensure Cython is installed
# sys.path.append Cython build if necessary
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
  cmdclass = { 'build_ext': build_ext},
  ext_modules = [Extension("genome_sam_collapser", ["genome_sam_collapser.pyx"])]
)
