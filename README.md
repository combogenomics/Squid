![squidd](https://user-images.githubusercontent.com/69002653/153582716-f780928e-b976-4da8-a814-82e1f563de44.svg)

# Introduction  
Squid is a sequencing reads ungapped mapper and splitter, designed for high performance and customizability. It is highly intuitive, supports multithreading and it does not need a separate step for index building. Here are some of its key features:
- 6.4x faster than Bowtie2;
- outputs BED / BEDPE for fast downstream analyses;
- outputs FASTQ for reads that map or do not map to an input database;
- compatible with the most commonly used sequencing library preparation protocols;

## Usage  
```console
Usage: squid -i <str> -R1 <str> [-R2 <str>] -l <str> -o <str> [Options]

Mandatory arguments:
   -i         <str>         input database in FASTA format (can be gzipp'd)
   -R1        <str>         read in forward direction (R1) (can be gzipp'd)
   -R2        <str>         read in reverse direction (R2) (can be gzipp'd)
   -o         <str>         basename, Squid will add "_R1.fastq", "_R2.fastq" and/or ".bed"

   "-l" argument is also mandatory, with one of the following format strings:
   -l SF                    Stranded Forward. R1 comes from the forward strand, R2 from the reverse strand
   -l SR                    Stranded Reverse. R1 comes from the reverse strand, R2 from the forward strand
   -l U                     Unstranded. R1 or R2 can derive from both strands.
   -l ISF                   Inward Stranded Forward. R1 and R2 behave as in SF. R1 must map upstream to R2
   -l ISR                   Inward Stranded Reverse. R1 and R2 behave as in SR. R1 must map downstream to R2
   -l IU                    Inward Unstranded. R1 and R2 behave as in U. With this option Squid tries ISF and ISR
   -l OSF                   Outward Stranded Forward. R1 and R2 behave as in SF. R1 must map downstream to R2
   -l OSR                   Outard Stranded Reverse. R1 and R2 behave as in SR. R1 must map upstream to R2
   -l OU                    Outard Unstranded. R1 and R2 behave as in U. With this option Squid tries OSF and OSR

Squid also provides a number of additional arguments for a more flexible mapping.

Boolean arguments:
   --diff                   when FASTQ(s) output is enabled, return reads that do not map to database.
                            By default this is switched off, meaning that only mapping reads will be written.
   --disjoin                when database is a multi-FASTA, allow R1 and R2 to map to different sequences.
                            Default is to coerce R1 and R2 to map to the the same seqid. When on,
                            disjoined read pairs will switch the score field in the BEDPE to 1 instead of 0.
   --ignore_N               do not treat Ns as mismatches, simply ignore them (default: OFF)
   --mask-lower             do not capitalize lowercase letters in database (default is to make them uppercase)
   --no-bed                 do not produce BED/BEDPE output file (default is to write it)
   --no-fastq               do not produce FASTQ output file(s) (default is to write them)
   --quiet                  do not print log to stderr (default is to be verbose)

Scanning and performance arguments:
   -e         <int>         evaluate <int> number of alternative positionings of R1 and R2, looking for a better match.
                            Default is to break as soon as a suitable match is found (-e 0). This option is meaningful
                            when BED/BEDPE output is enabled. Greater values of <int> affect performance but could
                            report a more accurate mapping when higly similar sequences are present in the database.
   -k         <int>         kmer size: 9, 11, 13 or 15 (default: 15)
   -m         <int>         max % of mismatches allowed during ungapped extension
                            Default is to force 85% sequence identity, hence -m 15.
   -s         <int>         step size while sliding over the sequencing reads
                            for a perfect match of length k.
                            Lower s increases sensitivity but decreases speed.
                            Min=1 (sliding window of 1), default: 17.
   -t         <int>         number of threads (default: 1)

```
