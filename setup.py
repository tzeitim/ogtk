import setuptools
#from setuptools_rust import Binding, RustExtensioni

with open("README.md", "r") as fh:
    long_description = fh.read()

print('install pysam via conda since zlib-devel cannot be installed using pip')
packages = setuptools.find_packages()
packages.append('rogtk')
setuptools.setup(
    name="ogtk",
    version="0.0.1",
    author="polivares",
    author_email="pedroy.final@gmail.com",
    description="general tools for genomics",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tzeitim/ogtk",
    #packages=setuptools.find_packages(),
    packages=packages,
    #rust_extensions=[RustExtension("rogtk.rogtk", binding=Binding.PyO3)],
    # rust extensions are not zip safe, just like C-extensions.
    #zip_safe=False,
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires = [
        'bgzip',
        'colorhash',
        #'fasta',
        'fastcluster',
        'matplotlib',
        'pandas',
        'polars',
        'pyarrow',
        'pyaml', 
        'pysam', 
        'pyfaidx', 
        'pyfasta', 
        'regex', 
        'rich',
        'scipy',
#        'metacells',
#        'scanpy',
        'seaborn',
        'setuptools_rust',
        'hdbscan',
        'ngs_tools',
        'pyseq_align',
        'tables'
    ], 
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
