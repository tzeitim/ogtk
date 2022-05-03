import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ogtk",
    version="0.0.1",
    author="polivares",
    author_email="pedroy.final@gmail.com",
    description="general tools for genomics",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tzeitim/ogtk",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires = ['fastcluster',
                        'colorhash',
                        'pysam', 
                        'pyfaidx', 
                        'regex', 
                        'pyfasta', 
                        'fasta', 
                        'pyaml'],
    python_requires='>=3.6',
)
# python3 setup.py sdist bdist_wheel

#from distutils.core import setup
#
#setup(name='ogtk',
#      version='0.1',
#      description='general tools',
#      author='Pedro Olivares',
#      author_email='pedro.olivares@mdc-berlin.de',
#      #packages=['UM'],
#     )
