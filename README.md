# Compare Assemblies
Python script for comparing two genome assemblies

* get scaffolding mapping using minimap2 (done)
* write out mapping results (done)
* get scaffold statistics: length, gaps, repeat content, GC content (done)
* support cases with partial synteny between assemblies (to do)
* get contig metrics (to do)


```
usage: compareAssemblies.py [-h] -r ASM1FILENAME -q ASM2FILENAME
                            [-m MINQUERYLEN] [-o OUTPUTPREFIX] [-d]

compare assemblies with minimap2

optional arguments:
  -h, --help       show this help message and exit
  -r ASM1FILENAME  input assembly 1 (ref) fasta
  -q ASM2FILENAME  input assembly 2 (query) fasta
  -m MINQUERYLEN   minimum sequence length in assembly 2 (query), default=20Mb
  -o OUTPUTPREFIX  output prefix, default=output
  -d               turn off saving minimap2 results

```