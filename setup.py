import setuptools

version = [l.strip() for l in open("pathogenprofiler/__init__.py") if "version" in l][0].split('"')[1]

setuptools.setup(

	name="pathogen-profiler",
	version=version,
	packages=["pathogenprofiler",],
	license="GPL3",
	long_description="Pathogen profiling tool",
	scripts=[
		'scripts/combine_vcf_variants.py',
		'scripts/rename_vcf_chrom.py',
		'scripts/add_dummy_AD.py'
	]
)
