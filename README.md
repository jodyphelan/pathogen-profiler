# Pathogen-Profiler
 This library gives acces to classes and functions create a profiling tool to look for mutations from NGS data. This library is used as the scaffolding for [TBProfiler2](https://github.com/jodyphelan/TBProfiler2).

### Installation
The library is written in pure python and has no dependancies to install:
```
git clone https://github.com/jodyphelan/pathogen-profiler.git
python setup.py install
```
At runtime, some of the functions will need tools involved with trimming, mapping and variant calling. Please install the following tools: *trimmomatic, bwa, minimap2, bowtie2, samtools, bcftools and parallel*.
