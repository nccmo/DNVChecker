#!/usr/bin/env python

from distutils.core import setup

setup(name='dnvchecker',
      version='0.1.0',
      description='Python tools to curate DNP, TNP from MAF files.',
      author='Yuki Saito',
      author_email='js3050ys1990@gmail.com',
      url='---',
      package_dir = {'': 'lib'},
      packages=['dnvchecker'],
      scripts=['dnvchecker'],
      license='GPL-3'
     )
