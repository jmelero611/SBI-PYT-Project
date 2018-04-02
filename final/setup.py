from distutils.core import setup

setup(name='Complex_builder',
	version='1.0',
	description='Reconstruction of a macro-comex',
	author='Lydia Fortea, Juan Luis Melero',
	packages=['modules'],
	scripts=['complex_reconstruction.py'],
	requires=['Bio', 're', 'sys', 'os', 'subprocess', 'argparse', 'string', 'copy', 'numpy'])
