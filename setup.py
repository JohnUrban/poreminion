import os
import sys
from setuptools import setup

version_py = os.path.join(os.path.dirname(__file__), 'poreminion', 'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"','')
long_description = """
``poreminion`` provides additional tools on top of poretools for working with nanopore sequencing data'
"""

with open("requirements.txt", "r") as f:
    install_requires = [x.strip() for x in f.readlines()]

setup(
        name="poreminion",
        version=version,
        install_requires=install_requires,
        requires = ['python (>=2.7, <3.0)'],
        packages=['poreminion',
                  'poreminion.scripts'],
        author="John Urban (poreminion); Nick Loman and Aaron Quinlan (poretools)",
        description='Additional tools on top of poretools for working with nanopore sequencing data',
        long_description=long_description,
        url="https://github.com/JohnUrban/poreminion",
        package_dir = {'poreminion': "poreminion"},
        package_data = {'poreminion': []},
        zip_safe = False,
        include_package_data=True,
        #scripts = ['poreminion/scripts/poreminion'],
        scripts = ['poreminion/scripts/getAttributes.sh'],
        entry_points = {
            'console_scripts' : [
                 'poreminion = poreminion.poreminion_main:main', 
            ],
        },  
        author_email="mr.john.urban@gmail.com",
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License (GPL)',
            'Topic :: Scientific/Engineering :: Bio-Informatics']
    )
