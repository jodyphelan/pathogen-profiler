import setuptools


setuptools.setup(

	name="pathogen-profiler",
	version="0.1dev",
	packages=["pathogenprofiler",],
	license="MIT",
	long_description="TBProfiler command line tool",
	scripts=['scripts/splitchr.py']
)
