# Pathogen-Profiler
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/pathogen-profiler/README.html) [![Anaconda-Server Badge](https://anaconda.org/bioconda/pathogen-profiler/badges/license.svg)](https://anaconda.org/bioconda/pathogen-profiler) [![Anaconda-Server Badge](https://anaconda.org/bioconda/pathogen-profiler/badges/latest_release_date.svg)](https://anaconda.org/bioconda/pathogen-profiler)

 This library gives acces to classes and functions create a profiling tool to look for mutations from NGS data. This library is used as the scaffolding for [TBProfiler](https://github.com/jodyphelan/TBProfiler).

### Installation

#### Conda

```
conda install -c bioconda pathogen-profiler
```
#### Manual
The library is written in pure python and has no dependancies to install:
```
git clone https://github.com/jodyphelan/pathogen-profiler.git
python setup.py install
```
At runtime, some of the functions will need tools involved with trimming, mapping and variant calling. Please install the following tools: *trimmomatic, bwa, minimap2, bowtie2, samtools, bcftools and parallel*.
