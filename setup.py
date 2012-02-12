from distutils.core import setup
setup(name='PyPLIF',
    version='0.1-alpha',
    author='Muhammad Radifar',
    author_email='m.radifar05@gmail.com',
    packages=['pyplif']
    py_modules=['PyPLIF', 'tanimoto_coef', 'ring', 'atom_property', 'interactions'],
    url=['http://code.google.com/p/pyplif/']
    license=['LICENSE.txt']
    description=['PyPLIF is a program/script written in Python to analyze protein-ligand interaction from the molecular docking result']
    long_description=open('README.txt').read(),
    install_requires=[
    "openbabel >= 0.6",
    "numpy >= 1.5.0",
    "bitarray >= 0.3.0",
    ],
)
