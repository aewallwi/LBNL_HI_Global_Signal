"""A setuptools based setup module.
See:
https://packaging.python.org/guides/distributing-packages-using-setuptools/
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup, find_packages
import pathlib
import os
import sys

sys.path.append("LBNL_HI_Global_Signal")

def package_files(package_dir, subdirectory):
    # walk the input package_dir/subdirectory
    # return a package_data list
    paths = []
    directory = os.path.join(package_dir, subdirectory)
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            path = path.replace(package_dir + '/', '')
            paths.append(os.path.join(path, filename))
    return paths

here = pathlib.Path(__file__).parent.resolve()
long_description = (here / 'README.md').read_text(encoding='utf-8')
data_files = package_files('LBNL_HI_Global_Signal', 'data')

setup(
    name='LBNL_HI_Global_Signal',  # Required
    version='0.0.1',  # Required
    description='Global Signal Simulation Tools',  # Optional
    long_description=long_description,  # Optional
    long_description_content_type='text/markdown',  # Optional (see note above)
    url='https://github.com/aewallwi/garray21cm',  # Optional
    author='A. Ewall-Wice',  # Optional
    author_email='aaronew@berkeley.edu',  # Optional
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3 :: Only',
    ],
    keywords='21cm, cosmology, foregrounds, radio astronomy, cosmic dawn',
    package_dir={'LBNL_HI_Global_Signal': 'LBNL_HI_Global_Signal'},
    packages=['LBNL_HI_Global_Signal'],
    python_requires='>=3.6, <4',
    install_requires=['pyuvdata==2.1.5',
                      'numpy',
                      'pygdsm @ git+git://github.com/telegraphic/pygdsm',
                      'uvtools @ git+git://github.com/HERA-Team/uvtools',
                      'pyuvsim @ git+git://github.com/RadioAstronomySoftwareGroup/pyuvsim',
                      'hera_sim @ git+git://github.com/HERA-Team/hera_sim',
                      ],
    include_package_data=True,
    package_data={'LBNL_HI_Global_Signal': data_files},
    exclude = ['tests'],
    zip_safe = False,
    )
