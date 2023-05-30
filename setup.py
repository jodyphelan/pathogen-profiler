import setuptools
import glob

version = [l.strip() for l in open("pathogenprofiler/__init__.py") if "version" in l][0].split('"')[1]

setuptools.setup(

	name="pathogen-profiler",
	version=version,
	packages=["pathogenprofiler",],
	license="GPL3",
	long_description="Pathogen profiling tool",
	scripts= glob.glob("scripts/*.py") 
)
