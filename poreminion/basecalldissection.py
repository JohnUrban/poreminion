## events and base calling

## similar to fragstats
## get:
#name,
#numevents,
#hascomp
#hastemp
#num temp,
#num comp,
#num called temp,
#num called comp,
#temp start/end,
#comp start end
# 2d align length
# 2d align score
# temp drift, shift, etc
# comp drift, shift, etc
## start time, end time


def get_2d_align_len(f5):
    return f5['/Analyses/Basecall_2D_000/Summary/hairpin_align/alignment_length']

def get_2d_align_score(f5):
    return f5['/Analyses/Basecall_2D_000/Summary/hairpin_align/alignment_score']

def get_strand_start_index(f5, strand):
    if strand == "template":
        return f5['/Analyses/Basecall_2D_000/Summary/split_hairpin/start_index_temp']
    elif strand == "complement":
        return f5['/Analyses/Basecall_2D_000/Summary/split_hairpin/start_index_comp']

def get_strand_end_index(f5, strand):
    if strand == "template":
        return f5['/Analyses/Basecall_2D_000/Summary/split_hairpin/end_index_temp']
    elif strand == "complement":
        return f5['/Analyses/Basecall_2D_000/Summary/split_hairpin/end_index_comp']
