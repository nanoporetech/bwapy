import os
import sys
import shutil
import re
import shutil
import platform
from glob import glob
from setuptools import setup, find_packages, Extension
from setuptools import Distribution, Command
from setuptools.command.install import install
import pkg_resources


__pkg_name__ = 'bwapy'
__author__ = 'cwright'
__description__ = 'Bindings to bwa alinger.'

__path__ = os.path.dirname(__file__)
__pkg_path__ = os.path.join(os.path.join(__path__, __pkg_name__))

# Get the version number from __init__.py
verstrline = open(os.path.join(__pkg_name__, '__init__.py'), 'r').read()
vsre = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(vsre, verstrline, re.M)
if mo:
    __version__ = mo.group(1)
else:
    raise RuntimeError('Unable to find version string in "{}/__init__.py".'.format(__pkg_name__))


# Get requirements from file, we prefer to have
#   preinstalled these with pip to make use of wheels.
dir_path = os.path.dirname(__file__)
install_requires = []
with open(os.path.join(dir_path, 'requirements.txt')) as fh:
    reqs = (
        r.split('#')[0].strip()
        for r in fh.read().splitlines() if not r.startswith('#')
    )
    for req in reqs:
        if req == '':
            continue
        if req.startswith('git+https'):
            req = req.split('/')[-1].split('@')[0]
        install_requires.append(req)

extra_requires = {}

extensions = []
extensions.append(Extension(
    'bwalib',
    sources=[os.path.join('bwapy', 'libbwapy.c'), os.path.join('bwapy', 'memopts.c')],
    include_dirs=['bwa'],
    extra_compile_args=['-pedantic', '-Wall', '-std=c99', '-march=native', '-ffast-math', '-DUSE_SSE2', '-DNDEBUG'],
    libraries=['z'],
    extra_objects=[os.path.join('bwa','libbwa.a')]
))


setup(
    name=__pkg_name__,
    version=__version__,
    url='https://github.com/nanoporetech/{}'.format(__pkg_name__),
    author=__author__,
    author_email='{}@nanoporetech.com'.format(__author__),
    description=__description__,
    dependency_links=[],
    ext_modules=extensions,
    install_requires=install_requires,
    tests_require=[].extend(install_requires),
    extras_require=extra_requires,
    # don't include any testing subpackages in dist
    packages=find_packages(exclude=['*.test', '*.test.*', 'test.*', 'test']),
    package_data={},
    zip_safe=False,
    data_files=[],
    entry_points={
        'console_scripts': [
            'bwamempy = {}.libbwa:main'.format(__pkg_name__)
        ]
    },
    scripts=[]
)
