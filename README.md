# HMM Gene Search Tools
This is a tool for automating HMM-based gene searches using profile hidden markov models. This program takes an input folder of HMMs & and an input folder of Prokka annotated genomes. It searches each genome with each HMM, recording the HMMER output in an intermediate directory (hmm_out). Then, it parses the HMMER output using SimpleHMMER and records all hits, their locus ID, and score/e-value for each pHMM in an hmm_hits directory, where each genome's hits are recorded in a genome-specific ouput csv. Hits beyond the TC score threshold (designated in the profile HMM) are counted and recorded in an output presence/absence matrix, where each genome is a row, and each input HMM is a column. This code all allows you to pull sequences (amino acid) from the hmm hits, recording them in an output .fa or .fna file for downstream phylogenetic analyses.
Some aspects of this work are based on the structure put forth by Gray Chadwick, which can be found at: [Phylogenomic Gene Cluster Display](https://github.com/gchadwick/phylogenomic_gene_cluster_display) 

## Summary Workflow: 
![hmm-tools_summary](https://github.com/kshalv/hmm_tools/assets/143134539/97743fe1-b064-44fd-a21e-361bf07759a1)


## Step 1: 
You'll need a set of genomes (organized in their own directory) and a set of HMMs of interest (organized in their own directory). For the genomes, make sure they have a .fna extension (see step #2 for caveat). The output will use whatever name the HMM file is labelled as (i.e. 'mcrA' for the 'mcrA.HMM' file). For my work, I have used genomes downloaded from GTDB (latest release [here](https://data.gtdb.ecogenomic.org/releases/) & annotated with [prokka](https://github.com/tseemann/prokka), as well as [Biopython](https://biopython.org/), [HMMER](http://hmmer.org/documentation.html) (installs with Prokka, otherwise will have to independently install), and [SimpleHMMER](https://github.com/minillinim/SimpleHMMER/tree/master).

## Step 2: 
The HMMER code runs on any set of genomes, so long as they have been annotated and use the ".fna" extension (though the script can be modified to take in other nucleotide file extensions). If you are using genomes supplied by GTDB, these should download with the "_genomic.fna" extension. This works best when you first annotate your genomes using Prokka; instructions and usage can be found [here](https://github.com/tseemann/prokka/)). If you would like to loop Prokka through a directory of genomes, use the following command (put forth by the prokka dev's [here](https://github.com/tseemann/prokka/issues/187)): 



    for F in ./*_genomic.fna; do  
    N=$(basename $F _genomic.fna) ;   
    prokka --gcode 15 --locustag $N --outdir $N --prefix $N  $F ; 
    done  
    


### Step 3: 
Running the code. This script has two options: 1. 'hmm' which loops through all hmm's, thresholds, & counts hits, & 2. 'seq' which uses the thresholded hmm hits to pull sequences from the prokka annotated .faa files for downstream phylogenetic analyses. You can customize an e-value for use here, otherwise it uses a default value of 1e-3.

To run '**hmm**':

``` 
usage: python hmm_tools.py hmm [-h] [-g GENOMES] [-hm HMMS] [-e EVALUE]

options:
  -h, --help            show this help message and exit
  -g GENOMES, --genomes GENOMES
                        Path to genome directory
  -hm HMMS, --hmms HMMS
                        Path to hmm directory. All files should have a .HMM extension and named
                        according to gene
  -e EVALUE, --evalue EVALUE
                        Specify e-value for parsing hmm output.

```



To run '**seq**':

```
usage: python hmm_tools.py seq [-h] [-hi HITS] [-hm HMMS] [-g GENOMES] [-s {nt,aa}] [-f {all,tc}]

options:
  -h, --help            show this help message and exit
  -hi HITS, --hits HITS
                        Path to hmm_hits. This is an intermediate directory produced by hmm searcher.
  -hm HMMS, --hmms HMMS
                        Path to hmm directory
  -g GENOMES, --genomes GENOMES
                        Path to genome directory
  -s {nt,aa}, --sequence {nt,aa}
                        Specify sequence type. Nucleotide (-nt) or amino acid (-aa)
  -f {all,tc}, --filter {all,tc}
                        Specify filtering threshold. "All" removes score threshold, whereas "tc" imposes
                        the TC cutoff from the seed HMM.
```
