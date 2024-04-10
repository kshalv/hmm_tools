# HMM Gene Search Tools
This is a tool for automating HMM-based gene searches using profile hidden markov models. This program takes an input folder of HMMs & and an input folder of Prokka annotated genomes. It searches each genome with each HMM, recording the HMMER output in an intermediate directory (hmm_out). Then, it parses the HMMER output using SimpleHMMER and records all hits, their locus ID, and score/e-value for each pHMM in an hmm_hits directory, where each genome's hits are recorded in a genome-specific ouput csv. Then, parsed hits beyond the TC score threshold designated in the profile HMM are counted and recorded in an output presence/absence matrix, where each genome is a row, and each input HMM is a column. This code all allows you to pull sequences (amino acid) from the hmm hits, recording them in an output .fa file for downstream phylogenetic analyses.
Some aspects of this work are loosely based on the structure put forth by Gray Chadwick, which can be found at: [Phylogenomic Gene Cluster Display](https://github.com/gchadwick/phylogenomic_gene_cluster_display) 

## Step 1: 
You'll need a set of genomes (organized in their own directory) and a set of HMMs of interest (organized in their own directory). For the genomes, make sure they have a .fna extension (see step #2 for caveat). The output will use whatever name the HMM file is labelled as (i.e. 'mcrA' for the 'mcrA.HMM' file). For my work, I have used genomes downloaded from GTDB (latest release [here](https://data.gtdb.ecogenomic.org/releases/) & annotated with [prokka](https://github.com/tseemann/prokka), as well as [Biopython](https://biopython.org/), [HMMER](http://hmmer.org/documentation.html) (installs with Prokka, otherwise will have to independently install), and [SimpleHMMER](https://github.com/minillinim/SimpleHMMER/tree/master).

## Step 2: 
The HMMER code runs on any set of genomes, so long as they use the ".fna" extension (though the script can be modified to take in other nucleotide file extensions). If you are using genomes supplied by GTDB, these should download with the "_genomic.fna" extension. This works best when you first annotate your genomes using Prokka; instructions and usage can be found here. If you would like to loop Prokka through a directory of genomes, use the following command (put forth by the prokka dev's [here](https://github.com/tseemann/prokka/issues/187)): 



    for F in ./*_genomic.fna; do  
    N=$(basename $F _genomic.fna) ;   
    prokka --gcode 15 --locustag $N --outdir $N --prefix $N  $F ; 
    done  
    


### Step 3: 
Run the code! This code has two major options: 1. 'hmm' which loops through all hmm's, thresholds, & counts hits, & 2. 'seq' which uses the thresholded hmm hits to pull sequences from the prokka annotated .faa files for downstream phylogenetic analyses. 

To run '**hmm**':

`` 
python hmm_tools.py hmm [-g genome directory path] [-hm HMM directory path] [-f {all, threshold}]
``



To run '**seq**':

``
python hmm_tools.py seq [-hi hits directory path] [-g genome directory path] [-hm HMM directory path]
``
