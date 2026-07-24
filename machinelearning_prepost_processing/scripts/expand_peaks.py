import pandas as pd
def clean_bed(bed_dir,output_dir):
    '''
    This function cleans duplicated peaks and merge close peaks
    '''
    import pybedtools as bt
    bed=bt.BedTool(bed_dir)
    bed.sort().merge(d=200).saveas(output_dir)


def expand_peaks(bed_dir,output_dir,expand_width):
    '''
    This function drops duplicated peaks, calculate the mid point of the peak and expand the peaks to 500 bp peaks. Eventually, a CSV file is saved.
    '''
    import pandas as pd
    #read bed

    bed = pd.read_csv(bed_dir,sep="\t",header=None)
    #drop duplication
    bed= bed.drop_duplicates(subset=[0,1,2], keep='first', inplace=False)
    #calculate the width of peaks
    bed["length"] = bed[2]-bed[1]
    #calculate the mid point
    bed["mid_point"] = bed["length"]//2+bed[1]
    #deepcopy
    bed_copy = bed.copy(deep=True)
    bed_copy[1]=bed_copy["mid_point"]- int(expand_width//2)
    bed_copy[2]=bed_copy["mid_point"]+ int(expand_width//2)
    bed_copy["length"] = bed_copy[2] - bed_copy[1]
    #output expanded peaks

    bed_copy.drop(["length","mid_point"],axis=1).to_csv(output_dir,index=False, sep="\t",header=None)

def expand_peaks(bed_dir,output_dir,expand_width):
    '''
    This function drops duplicated peaks, calculate the mid point of the peak and expand the peaks to 500 bp peaks. Eventually, a CSV file is saved.
    '''
    import pandas as pd
    #read bed

    bed = pd.read_csv(bed_dir,sep="\t",header=None)
    #drop duplication
    bed= bed.drop_duplicates(subset=[0,1,2], keep='first', inplace=False)
    #calculate the width of peaks
    bed["length"] = bed[2]-bed[1]
    #calculate the mid point
    bed["mid_point"] = bed["length"]//2+bed[1]
    #deepcopy
    bed_copy = bed.copy(deep=True)
    bed_copy[1]=bed_copy["mid_point"]- int(expand_width//2)
    bed_copy[2]=bed_copy["mid_point"]+ int(expand_width//2)
    bed_copy["length"] = bed_copy[2] - bed_copy[1]
    #output expanded peaks

    bed_copy.drop(["length","mid_point"],axis=1).to_csv(output_dir,index=False, sep="\t",header=None)

def fasta_to_df(dic_fasta):
    '''
    This function convert a fasta file to dataframe
    '''
    #create a dictionary
    seq_pair = {}
    for i in dic_fasta:
        key=i
        seq=dic_fasta[key].seq
        seq_pair[key]=str(seq)
    #create a dataframe
    df_seq=pd.DataFrame.from_dict(seq_pair,orient="index")
    df_seq[0]=df_seq[0].str.upper()
    return df_seq

def create_XY(fasta_dir,bed_dir,XY_file):
    from Bio import SeqIO
    import copy



    #read FASTA
    fasta=list(SeqIO.parse(fasta_dir,"fasta"))

    #Create reverse complementry sequence
    fasta_rev = copy.deepcopy(fasta)
    for i in range(len(fasta_rev)):
        r_seq = fasta[i].seq.complement()
        fasta_rev[i].seq = r_seq
        fasta_rev[i].id ="rev_" + fasta[i].id
    dic_fasta_rev = SeqIO.to_dict(fasta_rev)

    #create forward df
    fasta=SeqIO.parse(fasta_dir,"fasta")
    dic_fasta = SeqIO.to_dict(fasta)
    df_seq_f=fasta_to_df(dic_fasta)

    #create reverse df
    df_seq_r=fasta_to_df(dic_fasta_rev)
    df_rf=df_seq_f.append(df_seq_r)
    bed = pd.read_csv(bed_dir,sep="\t",header=None)
    bed[[1]] = bed[[1]].astype(str)
    bed[[2]] = bed[[2]].astype(str)
    bed["coordinate"]=bed[0].str.cat(bed[1], sep=":").str.cat(bed[2],sep="-")
    bed["coordinate_rev"]="rev_"+bed["coordinate"]
    bed=bed.set_index("coordinate")

    #bed=bed.set_index("coordinate_rev")
    UCLA3_f=pd.concat([df_rf, bed[6]], axis=1, join='inner')
    bed=bed.set_index("coordinate_rev")
    UCLA3_r=pd.concat([df_rf, bed[6]], axis=1, join='inner')

    #append r to f
    UCLA3_r.append(UCLA3_f).to_csv(XY_file,header=None)
