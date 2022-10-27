import setuptools
from setuptools import setup


with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='minedatabase',
      version='2.2.0',
      description='Metabolic In silico Network Expansions',
      long_description=long_description,
      long_description_content_type="text/markdown",
      url='https://github.com/tyo-nu/MINE-Database',
      author='Kevin Shebek, Jonathan Strutz',
      author_email='jonstrutz11@gmail.com',
      license='MIT',
      packages=setuptools.find_packages(exclude=["docs", "tests"]),
    #   install_requires=['pymongo'],#, 'rdkit', 'scikit-learn'],
      package_data={'minedatabase': ['data/*'],
                    'minedatabase.NP_Score': ['*.gz'],
                    'minedatabase.tests': ['data/*'],
                    },
      include_package_data=True,
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: MIT License',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
          'Topic :: Scientific/Engineering :: Chemistry',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.9',
          'Programming Language :: Python :: 3.8',
      ],
      )
