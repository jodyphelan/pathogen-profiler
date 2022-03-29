import setuptools


setuptools.setup(

	name="pathogen-profiler",
	version="2.0.3",
	packages=["pathogenprofiler",],
	license="GPL3",
	long_description="Pathogen profiling tool",
	scripts=[
		'scripts/combine_vcf_variants.py',
		'scripts/rename_vcf_chrom.py',
		'scripts/add_dummy_AD.py'
	]
)
