from setuptools import setup

setup(name='minedatabase',
      version='0.1',
      description='Metabolic In silico Network Expansions',
      url='http://github.com/JamesJeffryes/mine-database',
      author='James Jeffryes',
      author_email='jamesgjeffryes@gmail.com',
      license='MIT',
      packages=['minedatabase',
                'minedatabase.NP_Score'],
      install_requires=['pymongo'],
      extras_require={},
      )
