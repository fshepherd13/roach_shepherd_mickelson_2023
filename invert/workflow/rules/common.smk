import glob
import pandas as pd


def get_r1(wildcards):
    #given ids, expand with glob to retrieve the full file name with added nucleotides
    return glob.glob(config["in_dir"]+"/"+wildcards.sample+'_*R1*.fastq.gz')

def get_r2(wildcards):
    return glob.glob(config["in_dir"]+"/"+wildcards.sample+'_*R2*.fastq.gz')
    
def get_target_files(dictionary):
    list=[]
    for k,v in dictionary.items():
        for values in v:
            list.append("../results/invert_by_group/"+k+"/"+values+"_final_results.csv")
    return(list)

def get_cufflinks_labels(samples):
    return(",".join(samples))