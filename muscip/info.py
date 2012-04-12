"""This file contains defines parameters for muscip that we use to
fill settings in setup.py, the muscip top-level docstring, and for
building the docs.  In setup.py in particular, we exec this file, so
it cannot import muscip

( This info script has been adapted from the one found in nipype )

"""

_version_major = 0
_version_minor = 1
_version_micro = 0
_version_extra = '.dev'

def get_muscip_gitversion():
    """Muscip version as reported by the last commit in git"""
    import os
    import subprocess
    try:
        import muscip
        gitpath = os.path.realpath(os.path.join(os.path.dirname(muscip.__file__),
                                                os.path.pardir))
    except:
        gitpath = os.getcwd()
    gitpathgit = os.path.join(gitpath, '.git')
    if not os.path.exists(gitpathgit):
        return None
    ver = None
    try:
        o, _ = subprocess.Popen('git describe', shell=True, cwd=gitpath,
                                stdout=subprocess.PIPE).communicate()
    except Exception:
        pass
    else:
        ver = o.strip().split('-')[-1]
    return ver

if '.dev' in _version_extra:
    gitversion = get_muscip_gitversion()
    if gitversion:
        _version_extra = '.' + gitversion + '-' + 'dev'

# Format expected by setup.py and doc/source/conf.py: string of form "X.Y.Z"
__version__ = "%s.%s.%s%s" % (_version_major,
                              _version_minor,
                              _version_micro,
                              _version_extra)

CLASSIFIERS = ["Development Status :: 2 - Pre-Alpha",
               "Environment :: Console",
               "Intended Audience :: Science/Research",
               "License :: OSI Approved :: Apache Software License",
               "Operating System :: OS Independent",
               "Programming Language :: Python",
               "Topic :: Scientific/Engineering"]

description  = 'Utilities to aid in processing neuroimaging data.'

# Note: this long_description is actually a copy/paste from the top-level
# README.txt, so that it shows up nicely on PyPI.  So please remember to edit
# it only in one place and sync it correctly.
long_description = \
"""
========================================================
MUSCIP: MUSC Image Processing
========================================================

TODO: write long description
"""

# versions
NIBABEL_MIN_VERSION = '1.0'
NETWORKX_MIN_VERSION = '1.6'
NUMPY_MIN_VERSION = '1.6'
SCIPY_MIN_VERSION = '0.10'
NIPYPE_MIN_VERSION = '0.5.3'

NAME                = 'muscip'
MAINTAINER          = "Travis Nesland"
MAINTAINER_EMAIL    = "nesland@musc.edu"
DESCRIPTION         = description
LONG_DESCRIPTION    = long_description
URL                 = "http://github.com/tnez/muscip"
DOWNLOAD_URL        = "http://github.com/tnez/muscip/archives/master"
LICENSE             = "BSD License"
CLASSIFIERS         = CLASSIFIERS
AUTHOR              = "Travis Nesland"
AUTHOR_EMAIL        = "nesland@musc.edu"
PLATFORMS           = "OS Independent"
MAJOR               = _version_major
MINOR               = _version_minor
MICRO               = _version_micro
ISRELEASE           = _version_extra == ''
VERSION             = __version__
REQUIRES            = ["nibabel (>=1.0)", "networkx (>=1.6)", "numpy (>=1.6)",
                       "nipype (>=0.5)"]
STATUS              = 'pre-alpha'

# TODO: need to add pointer to my fork of nipype to in order to
# incorporate bug-fixes
