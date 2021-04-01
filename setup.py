from setuptools import setup, find_packages
from pathlib import Path

setup(
    name='ddqc',
    version='0.1.0',
    description='Description',
    url='https://github.com/ayshwaryas/ddqc',
    author='Author',
    author_email='email@example.gov',
    license='BSD 2-clause',
    packages=find_packages(),
    install_requires=[r.strip() for r in Path("requirements.txt").read_text("utf-8").splitlines()],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: POSIX :: Linux',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
    ],

    keywords="single cell/nucleus quality control",
)