from setuptools import setup
import numpy,os
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()
setup(
    name = "pytdt",
    version = "0.0.7",
    author = "Danny Antaki",
    author_email = "dantaki@ucsd.edu",
    description = ("python tdt"),
    license = "GPLv2",
    keywords = "tdt",
    url = "https://github.com/dantaki/pytdt/",
    packages=['pytdt'],
    package_dir={'pytdt':'pytdt'},
    long_description=read('README.md'),
    include_package_data=True,
    install_requires=['numpy','pandas','scipy','tqdm','statsmodels'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
    ],
)
