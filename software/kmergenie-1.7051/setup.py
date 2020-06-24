import subprocess
from setuptools import setup, find_packages, Extension

#http://stackoverflow.com/questions/15440115/how-would-i-run-a-script-file-as-part-of-the-python-setup-py-install
from setuptools.command.install import install

class MyInstall(install):

    def run(self):
        import os
        os.system('make')
        install.run(self)

setup(
    name='kmergenie',
    version='1.7051',
    author='kmergenie',
    license='Free Software License',
    packages=['kmergenie'],
    package_dir={'kmergenie': '.'},
    package_data={ 'kmergenie': ['scripts/test_install', 'scripts/decide', 'scripts/*.py', 'scripts/*.r', 'scripts/*.R', 'third_party/*.py', #'specialk', 
                                 'ntCard/ntcard', 'readfq.py']},
    zip_safe=False, # so that package data doesn't get zipped
    scripts=['kmergenie'],
    #include_package_data=True, 
    cmdclass={'install': MyInstall}
)
