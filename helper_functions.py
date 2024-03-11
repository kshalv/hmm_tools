#!/usr/bin/env python
import pandas as pd
import os
import shutil
import csv
import numpy as np
from simplehmmer.simplehmmer import HMMERRunner
from simplehmmer.simplehmmer import HMMERParser
from Bio import SeqIO

# loop through a directory of HMM's, against your directory of genomes
def hmmer_run(root_dir, hmm_dir):
    hmms = []
    archaea_genomes = []
    source = os.getcwd()
    output_dir = source+'/hmm_out/'
    root_dir = root_dir+'/'
    hmm_dir = hmm_dir+'/'
    os.makedirs(output_dir, exist_ok=True)

    for root, dirs, files in os.walk(root_dir):
        for file in files:
            if file.endswith(".fna"):
                genome_name = os.path.splitext(file)[0]
                archaea_genomes.append(genome_name)

    for file in os.listdir(hmm_dir):
        if not file.startswith('.DS_Store'):
            hmms.append(file.split('.HMM')[0])

    print('HMM Input List:\t'+str(hmms))
    print('Number of Genomes:  '+str(len(archaea_genomes))+'\n')

    for hmm in hmms: 
        HR = HMMERRunner(prefix=hmm)
        i = 0
        for archaea in archaea_genomes: 
            i += 1
            HR.search(hmm_dir+hmm+'.HMM',root_dir+archaea+'/'+archaea+'.faa',output_dir+archaea+'_hmm')  
        print('HMM Name: '+hmm+'\t\t'+'Genome Count: '+str(i))

    return output_dir
    
# parse the hmm output from hmmer
def hmmer_parser(root_dir, hmm_dir):
    dirs = []
    hmms = []
    source = os.getcwd()
    output_dir = source+'/hmm_hits/'
    root_dir = root_dir+'/'
    hmm_dir = hmm_dir+'/'
    os.makedirs(output_dir, exist_ok=True)

    # make a list of archaeal genomes from the folder of hmm output  
    for name in os.listdir(root_dir):
        if name != '.DS_Store': 
            dirs.append(name)
    # make a list of hmms
    for name in os.listdir(hmm_dir): 
        if name != '.DS_Store':
            name = name.split('.')[0]
            hmms.append(name)

    # print the number of archaeal genomes recorded in the hmm_dirs list. this should be the same as your total # of genomes
    print('number of genomes (input): ' + str(len(dirs)))

    i = 0
    for archaea in dirs:
        i += 1
        data = np.empty((0, 4), dtype='str')
        
        output_file = output_dir + archaea + '_hits.csv'
        besthmmhit = {}  
        bestevalhit = {}
        bestscore = {}
        # initialize count of hmm's
        h_count = 0
        
        for hmm in hmms:
            h_count += 1
            with open(root_dir + archaea + '/' + hmm + '_out.txt', 'r') as f:
                HP = HMMERParser(f)
                while True:
                    result = HP.next()
                    if result is None:
                        break       

                    gene = str(result).split('\t')[0]
                    if gene in besthmmhit:
                        if bestevalhit[gene] > float(str(result).split('\t')[6]):
                            bestevalhit[gene] = float(str(result).split('\t')[6])
                            bestscore[gene] = float(str(result).split('\t')[7])
                            besthmmhit[gene] = hmm
                    else:
                        if float(str(result).split('\t')[6]) < 1e-3:
                            bestevalhit[gene] = float(str(result).split('\t')[6])
                            bestscore[gene] = float(str(result).split('\t')[7])
                            besthmmhit[gene] = hmm
        
        with open(output_file, 'w') as f: 
            f.write('geneID,hmm,eval,score\n')
            for gene in besthmmhit:
                f.write(f"{gene},{besthmmhit[gene]},{bestevalhit[gene]},{bestscore[gene]}\n")
    
    print('number hmm count files (output): '+str(i))
    return output_dir


# filter hits by the score threshold designated in the profile hmm
def filt_count(hmm_dir, hit_dir, option): 
    hmm_dir = hmm_dir+'/'
    hit_dir = hit_dir+'/'
    hmm_paths = []
    score_dict = {}
    dict_list = []
    col_list = []

    output_file = './hmm_counts-filtering.csv'

    # records all hmm file paths in a list (hmm_paths), hmm name & threshold scores a dictionary (score_dict), & hmm names for output columns
    for file in os.listdir(hmm_dir): 
        if file != '.DS_Store':
            filepath=os.path.join(hmm_dir, file)
            hmm_paths.append(filepath)

            with open(filepath, 'r') as f: 
                hmm_name = filepath.split('/')[-1].split('.')[0]
                for line in f: 
                    if line.startswith('TC'):
                        score = float(line.strip().split(' ')[4])
                        score_dict.update({hmm_name:score})
                col_list.append(hmm_name)

    # initializes empty dataframe to record counts
    result_df = pd.DataFrame(columns=col_list)
    

    # loop through all files in the hmm hit directories and creates a genome specific dictionary that contains the hmm name & count
    # for all hmm hits that exceed the score threshold denoted by the .HMM file 
    i=0
    for file in os.listdir(hit_dir):
        if file.endswith('.csv'):
            i+=1
            filepath = os.path.join(hit_dir, file)
            df = pd.read_csv(filepath, sep=',')  
            
            part = filepath.split('/')[-1].split('_')
            genome_name = part[0]+'_'+part[1]
            genome_dict = {'genome_ID':genome_name}
    
            if option == 'all':
                for gene in score_dict:
                    df_filt = df[(df['score']>=score_dict[gene]) & (df['hmm']==gene)]
                                    
                    if df_filt.shape[0]>=1: 
                        grouped = df_filt.groupby('hmm').size().reset_index(name='Count')
                        ser = grouped['Count']
                        genome_dict.update({gene:ser[0]})
                    else: 
                        genome_dict.update({gene:0})      
                dict_list.append(genome_dict)
            
            else:
                has_hits = False
                
                for gene in score_dict:
                    df_filt = df[(df['score']>=score_dict[gene]) & (df['hmm']==gene)]
                    
                    if df_filt.shape[0]>=1: 
                        grouped = df_filt.groupby('hmm').size().reset_index(name='Count')
                        ser = grouped['Count']
                        genome_dict.update({gene:ser[0]})
                        has_hits = True
                    else: 
                        genome_dict.update({gene:0})

                if has_hits:         
                    dict_list.append(genome_dict)
        
    result_df = pd.DataFrame(dict_list)
    result_df.fillna(0, inplace=True)
    result_df.to_csv(output_file, index=False)
    print(score_dict)
    print('number of hmm hit files: '+str(i))
    print('number of output rows: '+str(result_df.shape[0]))


# optional function to pull sequences by hits
def seq_puller(hit_dir, hmm_dir, genome_dir):
    arch_list = []
    score_dict = {}
    hit_dict = {}
    hmms = []
    output_dir = os.getcwd()+'/seq_out/'
    hit_dir = hit_dir+'/'
    genome_dir = genome_dir+'/'

    #check to see if output directory exists, makes it if not
    os.makedirs(output_dir, exist_ok=True)
    
    #create file list of all the hit .csv's to loop through
    for file in os.listdir(hit_dir):
        if file!='.DS_Store':
            filepath = os.path.join(hit_dir, file)
            arch_list.append(filepath)

    #create score dictionary, we'll use this to threshold the hit dataframes
    for file in os.listdir(hmm_dir): 
        if file != '.DS_Store':
            filepath=os.path.join(hmm_dir, file)
            with open(filepath, 'r') as f: 
                hmm_name = filepath.split('/')[-1].split('.')[0]
                hmms.append(hmm_name)
                score = float(f.readlines()[15].split(' ')[4])
                score_dict.update({hmm_name:score})

    # loop through each hmm in the score dictionary
    for hmm in score_dict.keys():
        #count genome (i), number of sequnces there should be (based on the geneIDs pulled), & number of sequences written to your output file
        i=0
        seq_in = 0
        seq_writ = 0

        #intialize empty output file
        with open(output_dir+hmm+'_seq.fa', 'a+') as f:
            for archaea in arch_list: 
                i+=1
                genome = archaea.split('/')[-1].split('_h')[0]
                faa = genome_dir+genome+'/'+genome+'.faa'

                #create list of loci that fit the score cutoff for each hmm
                df = pd.read_csv(archaea, sep=',')
                df_filt = df.loc[np.where((df['hmm']==hmm) & (df['score']>=score_dict[hmm]))]
                loci = df_filt['geneID'].tolist()
                seq_in = seq_in + int(len(loci))

                #parse the fasta file for each genome and write the sequences that are present in the thresholded loci to a new fasta
                for locus in loci: 
                    for record in SeqIO.parse(faa, 'fasta'):
                        head = record.description
                        gene = head.split()[0][-6:]
                        if record.id == locus: 
                            f.write('>'+head+'\n')
                            f.write(str(record.seq)+'\n')
                            seq_writ += 1
                            
            print(hmm+'\t genomes:\t'+str(i))
            print('\t seq in:\t'+str(seq_in))
            print('\t seq out:\t'+str(seq_writ))
