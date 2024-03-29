#!/usr/bin/env python

# Code to 1) filter errors in oligo library synthesis or sequencing
# 2) extract barcode counts and 3) other info for qc and further analysis

import os, math, sys
from multiprocessing import Pool as ThreadPool
from datetime import datetime

from Bio import pairwise2
from Bio.Seq import Seq
from Bio import SeqIO
from collections import Counter
from collections import defaultdict

import pandas as pd
from natsort import natsorted

### helper functions ###
def exp_err(phred_arr):
    E = 0.0
    for q in phred_arr:
        E += math.pow(10,(-q/10))
    return E   

def tabulate(key, counts):
    counts[key] += 1

def tabulate_read(key,value,read):
    read[key].append(value)

def write_dict(dict, outputfile):
    # write sequences and counts to a text file
    print (outputfile)
    file = open(outputfile, 'w+')
    for w in sorted(dict, key=dict.get, reverse=True):
        file.write('{seq}, {num}\n'.format(seq=w, num=dict[w]))
    file.close()

def write_list(list, outputfile):
    # write sequences and counts to a text file
    print (outputfile)
    file = open(outputfile, 'w+')
    for w in list:
        file.write('{seq}\n'.format(seq=w))
    file.close()

def parse_dna(file):
    
#     fadapter = "GACCGGACCAC" #RC
#     radapter = "GACTCGCAGCA" #RC
    
    fadapter = 'TGCTGCGAGTC' #DNAseq_primer
    radapter = 'GTGGTCCGGTCT' #UTRseq

    #alignment scoring
    match = 0
    mismatch = -1
    gap_open = -1
    gap_extend = -1
    
    score_dist = []
    exp_err_dist = []
    read_trim_len_dist = []
    prom_len_dist = []
    
    count = 0    
    lowq_count = 0
    missingadapter_count = 0
    frag_count = 0
    badbc_count = 0
    goodbc_perfectalign_count = 0
    goodbc_goodalign_count = 0
    goodbc_badalign_count = 0
    
    lowq = defaultdict(int)
    missingadapter = defaultdict(int)
    frag = defaultdict(int)
    badbc = defaultdict(int)
    goodbc_perfectalign = defaultdict(list)
    goodbc_goodalign = defaultdict(list)
    goodbc_badalign = defaultdict(list)
    goodbc_perfectalign_bccounts = defaultdict(int)
    goodbc_goodalign_bccounts = defaultdict(int)

    print("Parsing of " + file + " started at:" + str(datetime.now()))

    for rec in SeqIO.parse(file, 'fastq'):
        count +=1
        qscore = rec.letter_annotations["phred_quality"]
        read = str(rec.seq)
        error = exp_err(qscore)
        exp_err_dist.append(error)
        #print('{}, r{}'.format(read,len(read)))

#        if count > 100:
#            break

        if error < 3:
            if fadapter in read and radapter in read:
                fpos = read.find(fadapter)
                rpos = read.find(radapter)
                read_trim = read[fpos+len(fadapter):rpos]
                read_trim_len_dist.append(len(read_trim))
                #print('{}, rt{}'.format(read_trim,len(read_trim)))
                if 165 >= len(read_trim) >= 70:
                    prom = read_trim[:-12]
                    prom_len_dist.append(len(prom))
                    bc = read_trim[-12:]
                    sbc = str(bc)
                    #print('{}, {}, p{}'.format(prom,bc,len(prom)))
                    if 155>=len(prom)>=55:
                        if sbc in sref_bc:
                            if prom == prom_lib[sbc]:
                                goodbc_perfectalign_count +=1
                                tabulate(sbc,goodbc_perfectalign_bccounts)
                                tabulate_read(sbc,prom,goodbc_perfectalign)
                                score_dist.append(0)
                            else:
                                score = pairwise2.align.globalms(prom,prom_lib[sbc],match,mismatch,gap_open,gap_extend,score_only = True)
                                score_dist.append(score)

                                if score >= (math.ceil(len(prom)*0.05)*-1):
                                    goodbc_goodalign_count +=1
                                    tabulate(sbc,goodbc_goodalign_bccounts)
                                    tabulate_read(sbc,prom,goodbc_goodalign)
                                else:
                                    goodbc_badalign_count +=1
                                    tabulate_read(sbc,prom,goodbc_badalign)
                        else:
                            badbc_count +=1
                            tabulate(read_trim,badbc)      
                    else:
                        badbc_count +=1
                        tabulate(read_trim,badbc)
                else:
                    frag_count +=1  
                    tabulate(read_trim,frag)
            else:
                missingadapter_count += 1
                tabulate(read,missingadapter)
        else:
            tabulate(read,lowq)
            lowq_count += 1

    combined = {}
    combined.update(goodbc_perfectalign_bccounts)
    combined.update(goodbc_goodalign_bccounts)
        
    combined_bccounts = {}
    for key in combined:
        if key in goodbc_goodalign_bccounts or key in goodbc_perfectalign_bccounts:
            combined_bccounts[key] = goodbc_goodalign_bccounts[key] + goodbc_perfectalign_bccounts[key]
    
    print(file + " log:")
    print("total count: " + str(count))
    print("lowq count: " +str(lowq_count) + " = " + "{0:.2f}".format(float(lowq_count/count)))
    print("missing_adapter count: "+ str(missingadapter_count) + " = " + "{0:.2f}".format(float(missingadapter_count/count)))
    print("frag count: " + str(frag_count) + " = " + "{0:.2f}".format(float(frag_count/count)))
    print("badbc count: " + str(badbc_count) + " = " + "{0:.2f}".format(float(badbc_count/count)))
    print("goodbc_badalignment count: " + str(goodbc_badalign_count) + " = " + "{0:.2f}".format(float(goodbc_badalign_count/count)))
    print("goodbc_goodalignment count: " + str(goodbc_goodalign_count) + " = " + "{0:.2f}".format(float(goodbc_goodalign_count/count)))
    print("goodbc_perfectalignment count: " + str(goodbc_perfectalign_count) + " = " + "{0:.2f}".format(float(goodbc_perfectalign_count/count)))
    print("")

    filename = file.split("/")[-1].split("-merged.fastq")[0]    
    d = pd.DataFrame(list(combined_bccounts.items()), columns = ['Barcode', 'Counts'])
    d = d.set_index('Barcode')
    d.to_csv(out_dir+"/01_bccounts/01_dna_bccounts/" + filename + "_bccounts.csv")
    d2 = pd.read_csv(out_dir+"/01_bccounts/01_dna_bccounts/" + filename + "_bccounts.csv")

    BC = ""
    BC_lst = []
    for BC_rev in d2['Barcode']:
        BC = Seq(BC_rev).reverse_complement()
        BC = ''.join(BC)
        BC_lst.append(BC)
    d2['Barcode'] = BC_lst
    d2.to_csv(out_dir+"/01_bccounts/01_dna_bccounts/" + filename + "_bccounts.csv", index=False)
    
    #write to files
    write_dict(combined_bccounts,out_dir+"/02_log/01_bccounts/" + filename + ".txt")
    write_dict(lowq,out_dir+"/02_log/02_lowq/" + filename + ".txt")
    write_dict(missingadapter,out_dir+"/02_log/03_missingadapter/" + filename + ".txt")
    write_dict(badbc,out_dir+"/02_log/04_badbc/" + filename + ".txt")
    write_dict(goodbc_badalign,out_dir+"/02_log/05_goodbc_badalign/" + filename + ".txt")
    write_dict(frag,out_dir+"/02_log/06_frag/" + filename + ".txt")
    write_dict(goodbc_perfectalign,out_dir+"/02_log/07_goodbc_perfectalign/" + filename + ".txt")
    write_dict(goodbc_goodalign,out_dir+"/02_log/08_goodbc_goodalign/" + filename + ".txt")
    write_dict(goodbc_perfectalign_bccounts,out_dir+"/02_log/09_goodbc_perfectalign_bccounts/" + filename + ".txt")
    write_dict(goodbc_goodalign_bccounts,out_dir+"/02_log/10_goodbc_goodalign_bccounts/" + filename + ".txt")

    with open(out_dir+"/02_log/11_log_files/" + filename + ".txt","w+") as f:
        f.write("total count: " + str(count) + "\n")
        f.write("lowq count: " +str(lowq_count) + " = " + "{0:.2f}".format(float(lowq_count/count)) + "\n")
        f.write("missing_adapter count: "+ str(missingadapter_count) + " = " + "{0:.2f}".format(float(missingadapter_count/count)) + "\n")
        f.write("frag count: " + str(frag_count) + " = " + "{0:.2f}".format(float(frag_count/count)) + "\n")
        f.write("badbc count: " + str(badbc_count) + " = " + "{0:.2f}".format(float(badbc_count/count)) + "\n")
        f.write("goodbc_badalignment count: " + str(goodbc_badalign_count) + " = " + "{0:.2f}".format(float(goodbc_badalign_count/count)) + "\n")
        f.write("goodbc_goodalignment count: " + str(goodbc_goodalign_count) + " = " + "{0:.2f}".format(float(goodbc_goodalign_count/count)) + "\n")
        f.write("goodbc_perfectalignment count: " + str(goodbc_perfectalign_count) + " = " + "{0:.2f}".format(float(goodbc_perfectalign_count/count)) + "\n")
       
    write_list(exp_err_dist,out_dir+"/02_log/11_log_files/exp_err_dist_" + filename + ".txt")
    write_list(read_trim_len_dist,out_dir+"/02_log/11_log_files/read_trim_len_dist_" + filename + ".txt")
    write_list(prom_len_dist,out_dir+"/02_log/11_log_files/prom_len_dist_" + filename + ".txt")
    write_list(score_dist,out_dir+"/02_log/11_log_files/score_dist_" + filename + ".txt")
        
    print("Parsing of " + file + " finished at:" + str(datetime.now()))

# number of threads for later multi-threading
num_threads = 36

DNA_dir = 'DNA_folder'
RNA_dir = 'RNA_folder'
out_dir = 'barcode_counts'

ref_csv = 'full_promoter_library_5reps_metadata.csv'

lib = pd.read_csv(ref_csv, index_col = 0) # barcode sequences in metadata is reverse complemented
lib = lib.set_index('Barcode')

def extractPromSeq(seq,f_adapter,r_adapter):
    seq_t = seq.upper()
    seq_t = seq_t.split(f_adapter,1)[1]
    seq_t = seq_t.split(r_adapter,1)[0][:-12]
    return seq_t

f_trim = 'TGCTGCGAGTC' #DNAseq_primer
r_trim = 'GTGGTCCGGTCT' #UTRseq

ref_bc = []
fwd_bc = ""
prom_lib = {}
for bc in lib.index:
#     fwd_bc = Seq(bc).reverse_complement()
    fwd_bc = Seq(bc)
    fwd_bc = ''.join(fwd_bc)
    ref_bc.append(fwd_bc)
#    prom_lib[fwd_bc] = lib['Sequence'][bc][18:-30]
    prom_lib[fwd_bc] = extractPromSeq(lib['Sequence'][bc],f_trim,r_trim)
sref_bc = set(ref_bc)

# list up DNA-seq data
DNA_file_list = []
for subdir, dirs, files in os.walk(DNA_dir):
    for file in files:
        if file.endswith("-merged.fastq"):
            DNA_file_list.append(os.path.join(subdir, file))
DNA_file_list = natsorted(DNA_file_list)

pool = ThreadPool(num_threads)
pairwise_results = pool.map(parse_dna,DNA_file_list)
pool.close()
pool.join()
