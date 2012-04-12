"""muscip : Image processing library developed to facilitate daily-dos.

( This setup script has been adapted from the one found in nipype)

"""

import sys
from glob import glob
from muscip import info

packages = ["muscip",
            "muscip.interfaces",
            "muscip.workflows",
            "muscip.workflows.connectivity"]

_info_fname = 'muscip/info.py'
INFO_VARS = {}
exec(open(_info_fname, 'rt').read(), {}, INFO_VARS)

################################################################################
# For some commands, use setuptools

if len(set(('develop', 'bdist_egg', 'bdist_rpm', 'bdist', 'bdist_dumb',
            'bdist_wininst', 'install_egg_info', 'egg_info', 'easy_install',
            )).intersection(sys.argv)) > 0:
    from setup_egg import extra_setuptools_args

# extra_setuptools_args can be defined from the line above, but it can
# also be defined here because setup.py has been exec'ed from
# setup_egg.py.
if not 'extra_setuptools_args' in globals():
    extra_setuptools_args = dict()


################################################################################
    
def main(**extra_args):
    from numpy.distutils.core import setup

    setup(name=INFO_VARS['NAME'],
          maintainer=INFO_VARS['MAINTAINER'],
          maintainer_email=INFO_VARS['MAINTAINER_EMAIL'],
          description=INFO_VARS['DESCRIPTION'],
          long_description=INFO_VARS['LONG_DESCRIPTION'],
          url=INFO_VARS['URL'],
          download_url=INFO_VARS['DOWNLOAD_URL'],
          license=INFO_VARS['LICENSE'],
          classifiers=INFO_VARS['CLASSIFIERS'],
          author=INFO_VARS['AUTHOR'],
          author_email=INFO_VARS['AUTHOR_EMAIL'],
          platforms=INFO_VARS['PLATFORMS'],
          version=INFO_VARS['VERSION'],
          requires=INFO_VARS['REQUIRES'],
          scripts = glob('bin/*'),
          packages = packages,
          **extra_args)
    
if __name__ == '__main__':
    main(**extra_setuptools_args)
