import setuptools


setuptools.setup(

	name="pathogen-profiler",
	version="1.5.0",
	packages=["pathogenprofiler",],
	license="GPL3",
	long_description="Pathogen profiling tool",
	scripts=[
		'scripts/splitchr.py',
		'scripts/add_dummy_AD.py',
		'scripts/pathogen-profiler-get-mutations.py',
		'scripts/pathogen-profiler-gather.py',
		'scripts/delly-overlaps.py'
	]
)
