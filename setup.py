from setuptools import setup 
import os
import glob

setup(name='pointLISA',
      version='0.1',
      description='Obtaining static and dynamic pointing vectors, angles and control',
      url='https://github.com/eabram/LISA/tree/master/pointLISA',
      author='Ester Abram',
      author_email='esterabram@hotmail.com',
      license='Nikhef/TNO',
      packages=['pointLISA'],
      package_dir={'pointLISA': 'pointLISA'},
      package_data={'pointLISA': ['parameters/Abram/*.txt','parameters/Waluschka/*.txt']},
      zip_safe=False)


