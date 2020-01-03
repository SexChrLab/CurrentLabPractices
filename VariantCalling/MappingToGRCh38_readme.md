**This file documents mapping to GRCh38 that takes into account alternative contigs**
- Read this tutorial: https://software.broadinstitute.org/gatk/documentation/article?id=8017

1. 1000G version of GRCh38
    1. The files were downloaded by Angela and are stored on the cluster here: `/data/storage/SAYRES/REFERENCE_GENOMES/GRCh38_reference_genome/`
    2. Download location: `ftp://ftp-trace.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/`.
  
2. Create a mask to mask:
    1. From `https://gatkforums.broadinstitute.org/gatk/discussion/7857/reference-genome-components`, it says: `The chrY location of PAR1 and PAR2 on GRCh38 are chrY:10,000-2,781,479 and chrY:56,887,902-57,217,415.`
    2. Working direcoty is: `/agavescratch/tphung3/SmallProjects/MappingToGRCh38/reference/`.
    3. `touch mask.bed`
    4. `echo -e chrY"\t"10000"\t"2781479 > mask.bed`
    5. `echo -e chrY"\t"56887902"\t"57217415 >> mask.bed`
  
3. Rename the alt index file
    - `bwa` requires that the basename of the alternative contigs have to match with the reference file.
    ```
    cp GRCh38_full_analysis_set_plus_decoy_hla.fa.alt GRCh38_full_analysis_set_plus_decoy_hla.XXonly.fa.alt
    ```
  
4. The snakefile `mapping_to_GRCh38.snakefile` has the rules for mapping to the 1000 Genome version of GRCh38 with and without taking into account alternative contig. Briefly, to map with taking into account alternative contig, just add the flag `-j` to the `bwa` command.

5. Compare the number of variants between mapping with and without taking into account alternative contigs
    1. I used fastq files from three exome from the placenta project.
    2. After mapping, I did joint-genotyping across the three placenta individuals
    3. Results:
      1. On the autosomes:
        1. There are **251,030** variants that were genotyped when mapping **without** taking into account alternative contigs
        2. There are **265,658** variants that were genotyped when mapping **with** taking into account alternative contigs
        3. There are **250,956** variants that are common between the two mapping methods. **250,530** of these variants match exactly in terms of genotypes
        4. There are **74** variants that were genotyped **only** when mapping **without** taking into account alternative contigs
        5. There are **14,702** variants that were genotyped **only** when mapping **with** taking into account alternative contigs
      2. On the X x_chromosome:
        1. There are **5,741** variants that were genotyped when mapping **without** taking into account alternative contigs
        2. There are **5,836** variants that were genotyped when mapping **with** taking into account alternative contigs
        3. There are **5,741** variants that are common between the two mapping methods. **5,737** of these variants match exactly in terms of genotypes
        5. There are **95** variants that were genotyped **only** when mapping **with** taking into account alternative contigs
