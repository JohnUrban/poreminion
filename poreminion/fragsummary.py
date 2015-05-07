from fragstats import *

## TODO 4/23/15-- break this hard-coded monster up into functions
## there is a lot of repeated pieces of code - tighten it up, shrink it

def summarize_fragstats(fragstats_df, extensive=False, timecheck=False):
    ## fragstats_df is dataframe from make_fragstats_dataframe()
    has2d = fragstats_df['has2d'] == 1
    hascomp = fragstats_df['hascomp'] == 1
    no2d = fragstats_df['has2d'] == 0
    temponly = fragstats_df['hascomp'] == 0
    
    n_molecules = len(fragstats_df['name'])
    n_temp_only = sum(temponly)
    n_comp = sum(hascomp)
    n_no_2d = sum(no2d)
    n_comp_has2d = sum(hascomp[has2d]) ## should be same as number with 2D
    n_comp_no2d = sum(hascomp[no2d])
    n_2d = sum(has2d)
    print "n_molecules\t" + str(n_molecules)
    print "n_template_only\t" + str(n_temp_only)
    print "pct_template_only\t" + str(100.0*n_temp_only/n_molecules)
    print "n_has_comp\t" + str(n_comp)
    print "pct_has_comp\t" + str(100.0*n_comp/n_molecules)
    print "n_has_2d\t" + str(n_2d)
    print "pct_has_2d\t" + str(100.0*n_2d/n_molecules)
    print "n_does_NOT_have_2d\t" + str(n_no_2d)
    print "pct_does_NOT_have_2d\t" + str(100.0*n_no_2d/n_molecules)
    print "pct_of_molecules_that_do_NOT_have_2D_because_template_only", 100.0*n_temp_only/n_no_2d
    print "pct_of_molecules_that_do_NOT_have_2D_that_DO_have_comlement", 100.0*n_comp_no2d/n_no_2d
    print "n_with_comp_that_DO_have_2D", n_comp_has2d, "(should be same as num with 2D)"
    print "n_with_comp_that_do_NOT_have_2D", n_comp_no2d
    print "pct_of_molecules_with_comp_that_DO_have_2D", 100.0*n_comp_has2d/n_comp
    print "pct_of_molecules_with_comp_that_do_NOT_have_2D", 100.0*n_comp_no2d/n_comp
    print
    

    ## MOLECULE
    sum_molecule_lengths = sum(fragstats_df['fragsize'])
    print "sum_molecule_lengths\t" + str(sum_molecule_lengths)
    x = [25,50,75]
    molecule_nx_values = NX(list(fragstats_df['fragsize']),x)
    for e in x:
        print "Molecule N%s\t%d" % (str(e), molecule_nx_values[e])
    print "mean_molecule_size", np.mean(fragstats_df['fragsize'])
    print "median_molecule_size", np.median(fragstats_df['fragsize'])
    print "max_molecule_size", max(fragstats_df['fragsize'])
    print "min_molecule_size", min(fragstats_df['fragsize'])
    n_molecules_gt_10kb = sum(fragstats_df['fragsize'] > 10e3)
    print "n_molecules_gt_10kb", n_molecules_gt_10kb
    print "pct_molecules_gt_10kb", 100.0*n_molecules_gt_10kb/n_molecules
    sum_molecules_gt_10kb = sum(fragstats_df['fragsize'][fragstats_df['fragsize'] > 10e3])
    print "sum_molecules_gt_10kb", sum_molecules_gt_10kb
    print "pct_summed_molecules_from_molecules_gt_10kb", 100.0*sum_molecules_gt_10kb/sum_molecule_lengths
    print
    

    ## HQ 2D
    q_2d_ge_9 = fragstats_df['meanscore2d'] >= 9
    sum_HQ_2d_lengths = sum(fragstats_df['seqlen2d'][has2d][q_2d_ge_9])
    print "sum_HQ_(Q>=9)_2d_lengths\t" + str(sum_HQ_2d_lengths)
    twod_HQ_nx_values = NX(list(fragstats_df['seqlen2d'][has2d][q_2d_ge_9]),x)
    for e in x:
        print "HQ_2D N%s\t%d" % (str(e), twod_HQ_nx_values[e])    
    mean_HQ_2d_length = fragstats_df['seqlen2d'][has2d][q_2d_ge_9].mean()
    print ("\t").join([str(e) for e in ["mean_HQ_2d_length", mean_HQ_2d_length]])
    median_HQ_2d_length = fragstats_df['seqlen2d'][has2d][q_2d_ge_9].median()
    print ("\t").join([str(e) for e in ["median_HQ_2d_length", median_HQ_2d_length]])
    print

    ## 2D
    sum_2d_lengths = sum(fragstats_df['seqlen2d'][has2d])
    print "sum_2d_lengths\t" + str(sum_2d_lengths)
    twod_nx_values = NX(list(fragstats_df['seqlen2d'][has2d]),x)
    for e in x:
        print "2D N%s\t%d" % (str(e), twod_nx_values[e])    
    mean_2d_length = fragstats_df['seqlen2d'][has2d].mean()
    print ("\t").join([str(e) for e in ["mean_2d_length", mean_2d_length]])
    median_2d_length = fragstats_df['seqlen2d'][has2d].median()
    print ("\t").join([str(e) for e in ["median_2d_length", median_2d_length]])
    print

    ## 1D
    sum_1d_lengths = sum(fragstats_df['seqlentemp'].append(fragstats_df['seqlencomp'][hascomp]))
    print "sum_1d_lengths\t" + str(sum_1d_lengths)   
    oned_nx_values = NX(list(fragstats_df['seqlentemp'].append(fragstats_df['seqlencomp'][hascomp])),x)
    for e in x:
        print "1D N%s\t%d" % (str(e), oned_nx_values[e])
    mean_1d_length = fragstats_df['seqlentemp'].append(fragstats_df['seqlencomp'][hascomp]).mean()
    print ("\t").join([str(e) for e in ["mean_1d_length", mean_1d_length]])
    median_1d_length = fragstats_df['seqlentemp'].append(fragstats_df['seqlencomp'][hascomp]).median()
    print ("\t").join([str(e) for e in ["median_1d_length", median_1d_length]])
    print

    ## Template
    sum_temp_lengths = sum(fragstats_df['seqlentemp'])
    print "sum_template_lengths\t" + str(sum_temp_lengths)   
    temp_nx_values = NX(list(fragstats_df['seqlentemp']),x)
    for e in x:
        print "Template N%s\t%d" % (str(e), temp_nx_values[e])
    mean_temp_length = fragstats_df['seqlentemp'].mean()
    print ("\t").join([str(e) for e in ["mean_template_length", mean_temp_length]])
    median_temp_length = fragstats_df['seqlentemp'].median()
    print ("\t").join([str(e) for e in ["median_template_length", median_temp_length]])
    print

    ## Complement
    sum_comp_lengths = sum(fragstats_df['seqlencomp'][hascomp])
    print "sum_complement_lengths\t" + str(sum_comp_lengths)   
    comp_nx_values = NX(list(fragstats_df['seqlencomp'][hascomp]),x)
    for e in x:
        print "Complement N%s\t%d" % (str(e), comp_nx_values[e])
    mean_comp_length = fragstats_df['seqlencomp'][hascomp].mean()
    print ("\t").join([str(e) for e in ["mean_complement_length", mean_comp_length]])
    median_comp_length = fragstats_df['seqlencomp'][hascomp].median()
    print ("\t").join([str(e) for e in ["median_complement_length", median_comp_length]])
    print

    ## Max, Min, and Q filtering:
    print ("\t").join(["metric", "length", "Q", "name", "read_type"])
    
    ## 2D -- TODO -- dont need to use logical indices after getting index (e.g. idxmax()) -- see 1D/T/C approach below
    min_2d_length_index = fragstats_df['fragsize'][has2d].idxmin()
    min_2d_length = fragstats_df['fragsize'][has2d][min_2d_length_index]
    min_2d_Q = fragstats_df['meanscore2d'][has2d][min_2d_length_index]
    min_2d_name = fragstats_df['name'][has2d][min_2d_length_index] 
    print ("\t").join([str(e) for e in ["min_2d_length", min_2d_length, min_2d_Q, min_2d_name, "2D"]])
    
    max_2d_length_index = fragstats_df['fragsize'][has2d].idxmax()
    max_2d_length = fragstats_df['fragsize'][has2d][max_2d_length_index]
    max_2d_Q = fragstats_df['meanscore2d'][has2d][max_2d_length_index]
    max_2d_name = fragstats_df['name'][has2d][max_2d_length_index] 
    print ("\t").join([str(e) for e in ["max_2d_length", max_2d_length, max_2d_Q, max_2d_name, "2D"]])
    
##    q_2d_ge_9 = fragstats_df['meanscore2d'] >= 9
    max_2d_Q_ge_9_index = fragstats_df['fragsize'][has2d][q_2d_ge_9].idxmax()
    max_2d_Q_ge_9_length = fragstats_df['fragsize'][has2d][q_2d_ge_9][max_2d_Q_ge_9_index]
    max_2d_Q_ge_9_Q = fragstats_df['meanscore2d'][has2d][q_2d_ge_9][max_2d_Q_ge_9_index]
    max_2d_Q_ge_9_name = fragstats_df['name'][has2d][q_2d_ge_9][max_2d_Q_ge_9_index]
    print ("\t").join([str(e) for e in ["max_2d_length_Q_ge_9", max_2d_Q_ge_9_length, max_2d_Q_ge_9_Q, max_2d_Q_ge_9_name, "2D"]])

    q_2d_ge_8_5 = fragstats_df['meanscore2d'] >= 8.5
    max_2d_Q_ge_8_5_index = fragstats_df['fragsize'][has2d][q_2d_ge_8_5].idxmax()
    max_2d_Q_ge_8_5_length = fragstats_df['fragsize'][has2d][q_2d_ge_8_5][max_2d_Q_ge_8_5_index]
    max_2d_Q_ge_8_5_Q = fragstats_df['meanscore2d'][has2d][q_2d_ge_8_5][max_2d_Q_ge_8_5_index]
    max_2d_Q_ge_8_5_name = fragstats_df['name'][has2d][q_2d_ge_8_5][max_2d_Q_ge_8_5_index]
    print ("\t").join([str(e) for e in ["max_2d_length_Q_ge_8.5", max_2d_Q_ge_8_5_length, max_2d_Q_ge_8_5_Q, max_2d_Q_ge_8_5_name, "2D"]])

    q_2d_ge_8 = fragstats_df['meanscore2d'] >= 8
    max_2d_Q_ge_8_index = fragstats_df['fragsize'][has2d][q_2d_ge_8].idxmax()
    max_2d_Q_ge_8_length = fragstats_df['fragsize'][has2d][q_2d_ge_8][max_2d_Q_ge_8_index]
    max_2d_Q_ge_8_Q = fragstats_df['meanscore2d'][has2d][q_2d_ge_8][max_2d_Q_ge_8_index]
    max_2d_Q_ge_8_name = fragstats_df['name'][has2d][q_2d_ge_8][max_2d_Q_ge_8_index]
    print ("\t").join([str(e) for e in ["max_2d_length_Q_ge_8", max_2d_Q_ge_8_length, max_2d_Q_ge_8_Q, max_2d_Q_ge_8_name, "2D"]])

    q_2d_ge_7_5 = fragstats_df['meanscore2d'] >= 7.5
    max_2d_Q_ge_7_5_index = fragstats_df['fragsize'][has2d][q_2d_ge_7_5].idxmax()
    max_2d_Q_ge_7_5_length = fragstats_df['fragsize'][has2d][q_2d_ge_7_5][max_2d_Q_ge_7_5_index]
    max_2d_Q_ge_7_5_Q = fragstats_df['meanscore2d'][has2d][q_2d_ge_7_5][max_2d_Q_ge_7_5_index]
    max_2d_Q_ge_7_5_name = fragstats_df['name'][has2d][q_2d_ge_7_5][max_2d_Q_ge_7_5_index]
    print ("\t").join([str(e) for e in ["max_2d_length_Q_ge_7.5", max_2d_Q_ge_7_5_length, max_2d_Q_ge_7_5_Q, max_2d_Q_ge_7_5_name, "2D"]])
    print


## 1D, Template, Complement
    #min
    min_template_length_index = fragstats_df['seqlentemp'].idxmin()
    min_complement_length_index = fragstats_df['seqlencomp'][hascomp].idxmin()
    
    min_template_length = fragstats_df['seqlentemp'][min_template_length_index]
    min_complement_length = fragstats_df['seqlencomp'][hascomp][min_complement_length_index]
    
    min_template_Q = fragstats_df['meanscoretemp'][min_template_length_index]
    min_complement_Q = fragstats_df['meanscorecomp'][hascomp][min_complement_length_index]
    
    min_template_name = fragstats_df['name'][min_template_length_index]
    min_complement_name = fragstats_df['name'][hascomp][min_complement_length_index]

    #max
    max_template_length_index = fragstats_df['seqlentemp'].idxmax()
    max_complement_length_index = fragstats_df['seqlencomp'][hascomp].idxmax()
    
    max_template_length = fragstats_df['seqlentemp'][max_template_length_index]
    max_complement_length = fragstats_df['seqlencomp'][hascomp][max_complement_length_index]
    
    max_template_Q = fragstats_df['meanscoretemp'][max_template_length_index]
    max_complement_Q = fragstats_df['meanscorecomp'][hascomp][max_complement_length_index]
    
    max_template_name = fragstats_df['name'][max_template_length_index]
    max_complement_name = fragstats_df['name'][hascomp][max_complement_length_index]


    #max, Q >= 4
    q_template_ge_4 = fragstats_df['meanscoretemp'] >= 4
    q_complement_ge_4 = fragstats_df['meanscorecomp'] >= 4
    
    max_template_length_Q_ge_4_index = fragstats_df['seqlentemp'][q_template_ge_4].idxmax()
    max_complement_length_Q_ge_4_index = fragstats_df['seqlencomp'][hascomp][q_complement_ge_4].idxmax()

    max_template_Q_ge_4_length = fragstats_df['seqlentemp'][max_template_length_Q_ge_4_index]
    max_complement_Q_ge_4_length = fragstats_df['seqlencomp'][hascomp][max_complement_length_Q_ge_4_index]

    max_template_Q_ge_4_Q = fragstats_df['meanscoretemp'][max_template_length_Q_ge_4_index]
    max_complement_Q_ge_4_Q = fragstats_df['meanscorecomp'][hascomp][max_complement_length_Q_ge_4_index]
    
    max_template_Q_ge_4_name = fragstats_df['name'][max_template_length_Q_ge_4_index]
    max_complement_Q_ge_4_name = fragstats_df['name'][hascomp][max_complement_length_Q_ge_4_index]

    #max, Q >= 3.5
    q_template_ge_3_5 = fragstats_df['meanscoretemp'] >= 3.5
    q_complement_ge_3_5 = fragstats_df['meanscorecomp'] >= 3.5
    
    max_template_length_Q_ge_3_5_index = fragstats_df['seqlentemp'][q_template_ge_3_5].idxmax()
    max_complement_length_Q_ge_3_5_index = fragstats_df['seqlencomp'][hascomp][q_complement_ge_3_5].idxmax()

    max_template_Q_ge_3_5_length = fragstats_df['seqlentemp'][max_template_length_Q_ge_3_5_index]
    max_complement_Q_ge_3_5_length = fragstats_df['seqlencomp'][hascomp][max_complement_length_Q_ge_3_5_index]

    max_template_Q_ge_3_5_Q = fragstats_df['meanscoretemp'][max_template_length_Q_ge_3_5_index]
    max_complement_Q_ge_3_5_Q = fragstats_df['meanscorecomp'][hascomp][max_complement_length_Q_ge_3_5_index]
    
    max_template_Q_ge_3_5_name = fragstats_df['name'][max_template_length_Q_ge_3_5_index]
    max_complement_Q_ge_3_5_name = fragstats_df['name'][hascomp][max_complement_length_Q_ge_3_5_index]


    #max, Q >= 3
    q_template_ge_3 = fragstats_df['meanscoretemp'] >= 3
    q_complement_ge_3 = fragstats_df['meanscorecomp'] >= 3
    
    max_template_length_Q_ge_3_index = fragstats_df['seqlentemp'][q_template_ge_3].idxmax()
    max_complement_length_Q_ge_3_index = fragstats_df['seqlencomp'][hascomp][q_complement_ge_3].idxmax()

    max_template_Q_ge_3_length = fragstats_df['seqlentemp'][max_template_length_Q_ge_3_index]
    max_complement_Q_ge_3_length = fragstats_df['seqlencomp'][hascomp][max_complement_length_Q_ge_3_index]

    max_template_Q_ge_3_Q = fragstats_df['meanscoretemp'][max_template_length_Q_ge_3_index]
    max_complement_Q_ge_3_Q = fragstats_df['meanscorecomp'][hascomp][max_complement_length_Q_ge_3_index]
    
    max_template_Q_ge_3_name = fragstats_df['name'][max_template_length_Q_ge_3_index]
    max_complement_Q_ge_3_name = fragstats_df['name'][hascomp][max_complement_length_Q_ge_3_index]

    #max, Q >= 2.5
    q_template_ge_2_5 = fragstats_df['meanscoretemp'] >= 2.5
    q_complement_ge_2_5 = fragstats_df['meanscorecomp'] >= 2.5
    
    max_template_length_Q_ge_2_5_index = fragstats_df['seqlentemp'][q_template_ge_2_5].idxmax()
    max_complement_length_Q_ge_2_5_index = fragstats_df['seqlencomp'][hascomp][q_complement_ge_2_5].idxmax()

    max_template_Q_ge_2_5_length = fragstats_df['seqlentemp'][max_template_length_Q_ge_2_5_index]
    max_complement_Q_ge_2_5_length = fragstats_df['seqlencomp'][hascomp][max_complement_length_Q_ge_2_5_index]

    max_template_Q_ge_2_5_Q = fragstats_df['meanscoretemp'][max_template_length_Q_ge_2_5_index]
    max_complement_Q_ge_2_5_Q = fragstats_df['meanscorecomp'][hascomp][max_complement_length_Q_ge_2_5_index]
    
    max_template_Q_ge_2_5_name = fragstats_df['name'][max_template_length_Q_ge_2_5_index]
    max_complement_Q_ge_2_5_name = fragstats_df['name'][hascomp][max_complement_length_Q_ge_2_5_index]

    # 1D
    if min_template_length < min_complement_length:
        print ("\t").join([str(e) for e in ["min_1d_length", min_template_length, min_template_Q, min_template_name, "template"]])
    else:
        print ("\t").join([str(e) for e in ["min_complement_length", min_complement_length, min_complement_Q, min_complement_name, "complement"]])
    if max_template_length > max_complement_length:
        print ("\t").join([str(e) for e in ["max_1d_length", max_template_length, max_template_Q, max_template_name, "template"]])
    else:
        print ("\t").join([str(e) for e in ["max_complement_length", max_complement_length, max_complement_Q, max_complement_name, "complement"]])
    if max_template_Q_ge_4_length > max_complement_Q_ge_4_length:
        print ("\t").join([str(e) for e in ["max_1d_Q_ge_4_length", max_template_Q_ge_4_length, max_template_Q_ge_4_Q, max_template_Q_ge_4_name, "template"]])
    else:
        print ("\t").join([str(e) for e in ["max_complement_Q_ge_4_length", max_complement_Q_ge_4_length, max_complement_Q_ge_4_Q, max_complement_Q_ge_4_name, "complement"]])
    if max_template_Q_ge_3_5_length > max_complement_Q_ge_3_5_length:
        print ("\t").join([str(e) for e in ["max_1d_Q_ge_3.5_length", max_template_Q_ge_3_5_length, max_template_Q_ge_3_5_Q, max_template_Q_ge_3_5_name, "template"]])
    else:
        print ("\t").join([str(e) for e in ["max_complement_Q_ge_3.5_length", max_complement_Q_ge_3_5_length, max_complement_Q_ge_3_5_Q, max_complement_Q_ge_3_5_name, "complement"]])
    if max_template_Q_ge_3_length > max_complement_Q_ge_3_length:
        print ("\t").join([str(e) for e in ["max_1d_Q_ge_3_length", max_template_Q_ge_3_length, max_template_Q_ge_3_Q, max_template_Q_ge_3_name, "template"]])
    else:
        print ("\t").join([str(e) for e in ["max_complement_Q_ge_3_length", max_complement_Q_ge_3_length, max_complement_Q_ge_3_Q, max_complement_Q_ge_3_name, "complement"]])
    if max_template_Q_ge_2_5_length > max_complement_Q_ge_2_5_length:
        print ("\t").join([str(e) for e in ["max_1d_Q_ge_2.5_length", max_template_Q_ge_2_5_length, max_template_Q_ge_2_5_Q, max_template_Q_ge_2_5_name, "template"]])
    else:
        print ("\t").join([str(e) for e in ["max_complement_Q_ge_2.5_length", max_complement_Q_ge_2_5_length, max_complement_Q_ge_2_5_Q, max_complement_Q_ge_2_5_name, "complement"]])
     
    print


    # Template
    print ("\t").join([str(e) for e in ["min_template_length", min_template_length, min_template_Q, min_template_name, "template"]])
    print ("\t").join([str(e) for e in ["max_template_length", max_template_length, max_template_Q, max_template_name, "template"]])
    print ("\t").join([str(e) for e in ["max_template_Q_ge_4_length", max_template_Q_ge_4_length, max_template_Q_ge_4_Q, max_template_Q_ge_4_name, "template"]])
    print ("\t").join([str(e) for e in ["max_template_Q_ge_3.5_length", max_template_Q_ge_3_5_length, max_template_Q_ge_3_5_Q, max_template_Q_ge_3_5_name, "template"]])
    print ("\t").join([str(e) for e in ["max_template_Q_ge_3_length", max_template_Q_ge_3_length, max_template_Q_ge_3_Q, max_template_Q_ge_3_name, "template"]])
    print ("\t").join([str(e) for e in ["max_template_Q_ge_2.5_length", max_template_Q_ge_2_5_length, max_template_Q_ge_2_5_Q, max_template_Q_ge_2_5_name, "template"]])
    print
    
    # Complement
    print ("\t").join([str(e) for e in ["min_complement_length", min_complement_length, min_complement_Q, min_complement_name, "complement"]])
    print ("\t").join([str(e) for e in ["max_complement_length", max_complement_length, max_complement_Q, max_complement_name, "complement"]])
    print ("\t").join([str(e) for e in ["max_complement_Q_ge_4_length", max_complement_Q_ge_4_length, max_complement_Q_ge_4_Q, max_complement_Q_ge_4_name, "complement"]])
    print ("\t").join([str(e) for e in ["max_complement_Q_ge_3.5_length", max_complement_Q_ge_3_5_length, max_complement_Q_ge_3_5_Q, max_complement_Q_ge_3_5_name, "complement"]])
    print ("\t").join([str(e) for e in ["max_complement_Q_ge_3_length", max_complement_Q_ge_3_length, max_complement_Q_ge_3_Q, max_complement_Q_ge_3_5_name, "complement"]])
    print ("\t").join([str(e) for e in ["max_complement_Q_ge_2.5_length", max_complement_Q_ge_2_5_length, max_complement_Q_ge_2_5_Q, max_complement_Q_ge_2_5_name, "complement"]])
    print

    ## longest 10
    ## 2D
    top10 = fragstats_df.sort(["seqlen2d","meanscore2d"], ascending=False)[:10][["seqlen2d","meanscore2d","name"]]
    seqlens = list(top10["seqlen2d"])
    meanscores = list(top10["meanscore2d"])
    names = list(top10["name"])
    print ("\t").join(["rank", "length", "Q", "name", "read_type", "analysis"])
    for i in range(len(names)):
        print ("\t").join([str(e) for e in [i+1, seqlens[i], meanscores[i], names[i], "2D", "Top_10_2D_reads"]])
    print

    ## Template
    top10 = fragstats_df.sort(["seqlentemp","meanscoretemp"], ascending=False)[:10][["seqlentemp","meanscoretemp","name"]]
    seqlens = list(top10["seqlentemp"])
    meanscores = list(top10["meanscoretemp"])
    names = list(top10["name"])
    print ("\t").join(["rank", "length", "Q", "name", "read_type", "analysis"])
    for i in range(len(names)):
        print ("\t").join([str(e) for e in [i+1, seqlens[i], meanscores[i], names[i], "template", "Top_10_Template_reads"]])
    print

    ## Complement
    top10 = fragstats_df.sort(["seqlencomp","meanscorecomp"], ascending=False)[:10][["seqlencomp","meanscorecomp","name"]]
    seqlens = list(top10["seqlencomp"])
    meanscores = list(top10["meanscorecomp"])
    names = list(top10["name"])
    print ("\t").join(["rank", "length", "Q", "name", "read_type", "analysis"])
    for i in range(len(names)):
        print ("\t").join([str(e) for e in [i+1, seqlens[i], meanscores[i], names[i], "complement", "Top_10_Complement_reads"]])
    print
    
    
## Q score distribution
    ## 2D
    mean_2d_Q = fragstats_df['meanscore2d'][has2d].mean()
    median_2d_Q = fragstats_df['meanscore2d'][has2d].median()
    std_2d_Q = fragstats_df['meanscore2d'][has2d].std()
    min_2d_Q_idx = fragstats_df['meanscore2d'][has2d].idxmin()
    max_2d_Q_idx = fragstats_df['meanscore2d'][has2d].idxmax()
    min_2d_Q = fragstats_df['meanscore2d'][min_2d_Q_idx]
    max_2d_Q = fragstats_df['meanscore2d'][max_2d_Q_idx]
    min_2d_Q_length = fragstats_df['seqlen2d'][min_2d_Q_idx]
    max_2d_Q_length = fragstats_df['seqlen2d'][max_2d_Q_idx]
    min_2d_Q_name = fragstats_df['name'][min_2d_Q_idx]
    max_2d_Q_name = fragstats_df['name'][max_2d_Q_idx]
    
    ## Template
    mean_temp_Q = fragstats_df['meanscoretemp'].mean()
    median_temp_Q = fragstats_df['meanscoretemp'].median()
    std_temp_Q = fragstats_df['meanscoretemp'].std()
    min_temp_Q_idx = fragstats_df['meanscoretemp'].idxmin()
    max_temp_Q_idx = fragstats_df['meanscoretemp'].idxmax()
    min_temp_Q = fragstats_df['meanscoretemp'][min_temp_Q_idx]
    max_temp_Q = fragstats_df['meanscoretemp'][max_temp_Q_idx]
    min_temp_Q_length = fragstats_df['seqlentemp'][min_temp_Q_idx]
    max_temp_Q_length = fragstats_df['seqlentemp'][max_temp_Q_idx]
    min_temp_Q_name = fragstats_df['name'][min_temp_Q_idx]
    max_temp_Q_name = fragstats_df['name'][max_temp_Q_idx]
    
    ## Complement
    mean_comp_Q = fragstats_df['meanscorecomp'][hascomp].mean()
    median_comp_Q = fragstats_df['meanscorecomp'][hascomp].median()
    std_comp_Q = fragstats_df['meanscorecomp'][hascomp].std()
    min_comp_Q_idx = fragstats_df['meanscorecomp'][hascomp].idxmin()
    max_comp_Q_idx = fragstats_df['meanscorecomp'][hascomp].idxmax()
    min_comp_Q = fragstats_df['meanscorecomp'][min_comp_Q_idx]
    max_comp_Q = fragstats_df['meanscorecomp'][max_comp_Q_idx]
    min_comp_Q_length = fragstats_df['seqlencomp'][min_comp_Q_idx]
    max_comp_Q_length = fragstats_df['seqlencomp'][max_comp_Q_idx]
    min_comp_Q_name = fragstats_df['name'][min_comp_Q_idx]
    max_comp_Q_name = fragstats_df['name'][max_comp_Q_idx]
    
    ## 1D
    mean_1d_Q = fragstats_df['meanscoretemp'].append(fragstats_df['meanscorecomp'][hascomp]).mean()
    median_1d_Q = fragstats_df['meanscoretemp'].append(fragstats_df['meanscorecomp'][hascomp]).median()
    std_1d_Q = fragstats_df['meanscoretemp'].append(fragstats_df['meanscorecomp'][hascomp]).std()
    if min_temp_Q < min_comp_Q:
        min_1d_Q_idx = min_temp_Q_idx
        min_1d_Q = min_temp_Q
        min_1d_Q_length = min_temp_Q_length
        min_1d_Q_name = min_temp_Q_name
        min_1d_read_type = "template"
    else:
        min_1d_Q_idx = min_comp_Q_idx
        min_1d_Q = min_comp_Q
        min_1d_Q_length = min_comp_Q_length
        min_1d_Q_name = min_comp_Q_name
        min_1d_read_type = "complement"
    if max_temp_Q > max_comp_Q:
        max_1d_Q_idx = max_temp_Q_idx
        max_1d_Q = max_temp_Q
        max_1d_Q_length = max_temp_Q_length
        max_1d_Q_name = max_temp_Q_name
        max_1d_read_type = "template"
    else:
        max_1d_Q_idx = max_comp_Q_idx
        max_1d_Q = max_comp_Q
        max_1d_Q_length = max_comp_Q_length
        max_1d_Q_name = max_comp_Q_name
        max_1d_read_type = "complement"
##    print ("\t").join(["read_type", "mean_Q", "median_Q", "std_dev_Q", "min_Q", "max_Q", "read_length_min_Q", "read_length_max_Q", "read_name_min_Q", "read_name_max_Q"])
##    print ("\t").join([str(e) for e in ["2D", mean_2d_Q, median_2d_Q, std_2d_Q, min_2d_Q, max_2d_Q]])
##    print ("\t").join([str(e) for e in ["1D", mean_1d_Q, median_1d_Q, std_1d_Q, min_1d_Q, max_1d_Q]])
##    print ("\t").join([str(e) for e in ["Template", mean_temp_Q, median_temp_Q, std_temp_Q, min_temp_Q, max_temp_Q]])
##    print ("\t").join([str(e) for e in ["Complement", mean_comp_Q, median_comp_Q, std_comp_Q, min_comp_Q, max_comp_Q]])
##    print
    print ("\t").join(["metric", "Q", "length", "read_type", "name"])
    print ("\t").join([str(e) for e in ["median_Q_2d", median_2d_Q, "-", "2D", "-"]])
    print ("\t").join([str(e) for e in ["mean_Q_2d", mean_2d_Q, "-", "2D", "-"]])
    print ("\t").join([str(e) for e in ["std_dev_Q_2d", std_2d_Q, "-", "2D", "-"]])
    print ("\t").join([str(e) for e in ["min_Q_2d", min_2d_Q, min_2d_Q_length, "2D", min_2d_Q_name]])
    print ("\t").join([str(e) for e in ["max_Q_2d", max_2d_Q, max_2d_Q_length, "2D", max_2d_Q_name]])
    print
    print ("\t").join(["metric", "Q", "length", "read_type", "name"])    
    print ("\t").join([str(e) for e in ["median_Q_1d", median_1d_Q, "-", "1D", "-"]])
    print ("\t").join([str(e) for e in ["mean_Q_1d", mean_1d_Q, "-", "1D", "-"]])
    print ("\t").join([str(e) for e in ["std_dev_Q_1d", std_1d_Q, "-", "1D", "-"]])
    print ("\t").join([str(e) for e in ["min_Q_1d", min_1d_Q, min_1d_Q_length, min_1d_read_type, min_1d_Q_name]])
    print ("\t").join([str(e) for e in ["max_Q_1d", max_1d_Q, max_1d_Q_length, max_1d_read_type, max_1d_Q_name]])
    print
    print ("\t").join(["metric", "Q", "length", "read_type", "name"])
    print ("\t").join([str(e) for e in ["median_Q_template", median_temp_Q, "-", "template", "-"]])
    print ("\t").join([str(e) for e in ["mean_Q_template", mean_temp_Q, "-", "template", "-"]])
    print ("\t").join([str(e) for e in ["std_dev_Q_template", std_temp_Q, "-", "template", "-"]])
    print ("\t").join([str(e) for e in ["min_Q_template", min_temp_Q, min_temp_Q_length, "template", min_temp_Q_name]])
    print ("\t").join([str(e) for e in ["max_Q_template", max_temp_Q, max_temp_Q_length, "template", max_temp_Q_name]])
    print
    print ("\t").join(["metric", "Q", "length", "read_type", "name"])
    print ("\t").join([str(e) for e in ["median_Q_complement", median_comp_Q, "-", "complement", "-"]])
    print ("\t").join([str(e) for e in ["mean_Q_complement", mean_comp_Q, "-", "complement", "-"]])
    print ("\t").join([str(e) for e in ["std_dev_Q_complement", std_comp_Q, "-", "complement", "-"]])
    print ("\t").join([str(e) for e in ["min_Q_complement", min_comp_Q, min_comp_Q_length, "complement", min_comp_Q_name]])
    print ("\t").join([str(e) for e in ["max_Q_complement", max_comp_Q, max_comp_Q_length, "complement", max_comp_Q_name]])
    print

    
    ## RATIO
    print "Template:Complement events ratio stats:"
    mean_tc_ratio_hascomp = fragstats_df['log2_tc_ratio'][fragstats_df['hascomp'] == 1].mean()
    min_tc_ratio_hascomp = fragstats_df['log2_tc_ratio'][fragstats_df['hascomp'] == 1].min()
    max_tc_ratio_hascomp = fragstats_df['log2_tc_ratio'][fragstats_df['hascomp'] == 1].max()
    mean_tc_ratio_has2d = fragstats_df['log2_tc_ratio'][fragstats_df['has2d'] == 1].mean()
    min_tc_ratio_has2d = fragstats_df['log2_tc_ratio'][fragstats_df['has2d'] == 1].min()
    max_tc_ratio_has2d = fragstats_df['log2_tc_ratio'][fragstats_df['has2d'] == 1].max()
    mean_tc_ratio_hascomp_no2d = fragstats_df['log2_tc_ratio'][fragstats_df['hascomp'] == 1][fragstats_df['has2d'] == 0].mean()
    min_tc_ratio_hascomp_no2d = fragstats_df['log2_tc_ratio'][fragstats_df['hascomp'] == 1][fragstats_df['has2d'] == 0].min()
    max_tc_ratio_hascomp_no2d = fragstats_df['log2_tc_ratio'][fragstats_df['hascomp'] == 1][fragstats_df['has2d'] == 0].max()
    print "mean_log2_tc_ratio_hascomp\t" + str(mean_tc_ratio_hascomp)
    print "min_log2_tc_ratio_hascomp\t" + str(min_tc_ratio_hascomp)
    print "max_log2_tc_ratio_hascomp\t" + str(max_tc_ratio_hascomp)
    print "mean_log2_tc_ratio_has2d\t" + str(mean_tc_ratio_has2d)
    print "min_log2_tc_ratio_has2d\t" + str(min_tc_ratio_has2d)
    print "max_log2_tc_ratio_has2d\t" + str(max_tc_ratio_has2d)
    print "mean_log2_tc_ratio_hascomp_no2d\t" + str(mean_tc_ratio_hascomp_no2d)
    print "min_log2_tc_ratio_hascomp_no2d\t" + str(min_tc_ratio_hascomp_no2d)
    print "max_log2_tc_ratio_hascomp_no2d\t" + str(max_tc_ratio_hascomp_no2d)
    print

    
    if extensive:
        #slope
        mean_slope = fragstats_df['slope'].mean()
        median_slope = fragstats_df['slope'].median()
        sd_slope = fragstats_df['slope'].std()
        min_slope_idx = fragstats_df['slope'].idxmin()
        max_slope_idx = fragstats_df['slope'].idxmax()
        print "median_slope", median_slope
        print "mean_slope", mean_slope
        print "sd_slope", sd_slope
        print
        print ("\t").join(["#metric", "slope_value", "numevents_from_molecule", "seq_len_2d", "seq_len_template", "seq_len_complement", "Q_2d", "Q_template", "Q_complement", "name"])
        print "min_slope", fragstats_df['slope'][min_slope_idx], fragstats_df['numevents'][min_slope_idx], nonetodash(fragstats_df['seqlen2d'][min_slope_idx]), nonetodash(fragstats_df['seqlentemp'][min_slope_idx]), nonetodash(fragstats_df['seqlencomp'][min_slope_idx]),  nonetodash(fragstats_df['meanscore2d'][min_slope_idx]), nonetodash(fragstats_df['meanscoretemp'][min_slope_idx]), nonetodash(fragstats_df['meanscorecomp'][min_slope_idx]), fragstats_df['name'][min_slope_idx]
        print "max_slope", fragstats_df['slope'][max_slope_idx], fragstats_df['numevents'][max_slope_idx], nonetodash(fragstats_df['seqlen2d'][max_slope_idx]), nonetodash(fragstats_df['seqlentemp'][max_slope_idx]), nonetodash(fragstats_df['seqlencomp'][max_slope_idx]), nonetodash(fragstats_df['meanscore2d'][max_slope_idx]), nonetodash(fragstats_df['meanscoretemp'][max_slope_idx]), nonetodash(fragstats_df['meanscorecomp'][max_slope_idx]), fragstats_df['name'][max_slope_idx]
        print
        
        # t moves -- NOTE: sum(all moves) = num_called_events
        # thus, 100*moves_x/num_called_events is percent of given move x
        med_pct_tevents_move_0 = 100.0*(fragstats_df['tmove_0']/fragstats_df['numcalledeventstemp']).median()
        mean_pct_tevents_move_0 = 100.0*(fragstats_df['tmove_0']/fragstats_df['numcalledeventstemp']).mean()
        std_pct_tevents_move_0 = 100.0*(fragstats_df['tmove_0']/fragstats_df['numcalledeventstemp']).std()
        minidx_pct_tevents_move_0 = (fragstats_df['tmove_0']/fragstats_df['numcalledeventstemp']).idxmin()
        maxidx_pct_tevents_move_0 = (fragstats_df['tmove_0']/fragstats_df['numcalledeventstemp']).idxmax()
        min_pct_tevents_move_0 = 100.0*(fragstats_df['tmove_0']/fragstats_df['numcalledeventstemp']).min()
        max_pct_tevents_move_0 = 100.0*(fragstats_df['tmove_0']/fragstats_df['numcalledeventstemp']).max()

        med_pct_tevents_move_1 = 100.0*(fragstats_df['tmove_1']/fragstats_df['numcalledeventstemp']).median()
        mean_pct_tevents_move_1 = 100.0*(fragstats_df['tmove_1']/fragstats_df['numcalledeventstemp']).mean()
        std_pct_tevents_move_1 = 100.0*(fragstats_df['tmove_1']/fragstats_df['numcalledeventstemp']).std()
        minidx_pct_tevents_move_1 = (fragstats_df['tmove_1']/fragstats_df['numcalledeventstemp']).idxmin()
        maxidx_pct_tevents_move_1 = (fragstats_df['tmove_1']/fragstats_df['numcalledeventstemp']).idxmax()
        min_pct_tevents_move_1 = 100.0*(fragstats_df['tmove_1']/fragstats_df['numcalledeventstemp']).min()
        max_pct_tevents_move_1 = 100.0*(fragstats_df['tmove_1']/fragstats_df['numcalledeventstemp']).max()
        
        med_pct_tevents_move_2 = 100.0*(fragstats_df['tmove_2']/fragstats_df['numcalledeventstemp']).median()
        mean_pct_tevents_move_2 = 100.0*(fragstats_df['tmove_2']/fragstats_df['numcalledeventstemp']).mean()
        std_pct_tevents_move_2 = 100.0*(fragstats_df['tmove_2']/fragstats_df['numcalledeventstemp']).std()
        minidx_pct_tevents_move_2 = (fragstats_df['tmove_2']/fragstats_df['numcalledeventstemp']).idxmin()
        maxidx_pct_tevents_move_2 = (fragstats_df['tmove_2']/fragstats_df['numcalledeventstemp']).idxmax()
        min_pct_tevents_move_2 = 100.0*(fragstats_df['tmove_2']/fragstats_df['numcalledeventstemp']).min()
        max_pct_tevents_move_2 = 100.0*(fragstats_df['tmove_2']/fragstats_df['numcalledeventstemp']).max()
        
        med_pct_tevents_move_3 = 100.0*(fragstats_df['tmove_3']/fragstats_df['numcalledeventstemp']).median()
        mean_pct_tevents_move_3 = 100.0*(fragstats_df['tmove_3']/fragstats_df['numcalledeventstemp']).mean()
        std_pct_tevents_move_3 = 100.0*(fragstats_df['tmove_3']/fragstats_df['numcalledeventstemp']).std()
        minidx_pct_tevents_move_3 = (fragstats_df['tmove_3']/fragstats_df['numcalledeventstemp']).idxmin()
        maxidx_pct_tevents_move_3 = (fragstats_df['tmove_3']/fragstats_df['numcalledeventstemp']).idxmax()
        min_pct_tevents_move_3 = 100.0*(fragstats_df['tmove_3']/fragstats_df['numcalledeventstemp']).min()
        max_pct_tevents_move_3 = 100.0*(fragstats_df['tmove_3']/fragstats_df['numcalledeventstemp']).max()
        
        med_pct_tevents_move_4 = 100.0*(fragstats_df['tmove_4']/fragstats_df['numcalledeventstemp']).median()
        mean_pct_tevents_move_4 = 100.0*(fragstats_df['tmove_4']/fragstats_df['numcalledeventstemp']).mean()
        std_pct_tevents_move_4 = 100.0*(fragstats_df['tmove_4']/fragstats_df['numcalledeventstemp']).std()
        minidx_pct_tevents_move_4 = (fragstats_df['tmove_4']/fragstats_df['numcalledeventstemp']).idxmin()
        maxidx_pct_tevents_move_4 = (fragstats_df['tmove_4']/fragstats_df['numcalledeventstemp']).idxmax()
        min_pct_tevents_move_4 = 100.0*(fragstats_df['tmove_4']/fragstats_df['numcalledeventstemp']).min()
        max_pct_tevents_move_4 = 100.0*(fragstats_df['tmove_4']/fragstats_df['numcalledeventstemp']).max()
        
        med_pct_tevents_move_5 = 100.0*(fragstats_df['tmove_5']/fragstats_df['numcalledeventstemp']).median()
        mean_pct_tevents_move_5 = 100.0*(fragstats_df['tmove_5']/fragstats_df['numcalledeventstemp']).mean()
        std_pct_tevents_move_5 = 100.0*(fragstats_df['tmove_5']/fragstats_df['numcalledeventstemp']).std()
        minidx_pct_tevents_move_5 = (fragstats_df['tmove_5']/fragstats_df['numcalledeventstemp']).idxmin()
        maxidx_pct_tevents_move_5 = (fragstats_df['tmove_5']/fragstats_df['numcalledeventstemp']).idxmax()
        min_pct_tevents_move_5 = 100.0*(fragstats_df['tmove_5']/fragstats_df['numcalledeventstemp']).min()
        max_pct_tevents_move_5 = 100.0*(fragstats_df['tmove_5']/fragstats_df['numcalledeventstemp']).max()
        
        
        # c moves -- add hascomp?
        med_pct_cevents_move_0 = 100.0*(fragstats_df['cmove_0']/fragstats_df['numcalledeventscomp']).median()
        mean_pct_cevents_move_0 = 100.0*(fragstats_df['cmove_0']/fragstats_df['numcalledeventscomp']).mean()
        std_pct_cevents_move_0 = 100.0*(fragstats_df['cmove_0']/fragstats_df['numcalledeventscomp']).std()
        minidx_pct_cevents_move_0 = (fragstats_df['cmove_0']/fragstats_df['numcalledeventscomp']).idxmin()
        maxidx_pct_cevents_move_0 = (fragstats_df['cmove_0']/fragstats_df['numcalledeventscomp']).idxmax()
        min_pct_cevents_move_0 = 100.0*(fragstats_df['cmove_0']/fragstats_df['numcalledeventscomp']).min()
        max_pct_cevents_move_0 = 100.0*(fragstats_df['cmove_0']/fragstats_df['numcalledeventscomp']).max()

        med_pct_cevents_move_1 = 100.0*(fragstats_df['cmove_1']/fragstats_df['numcalledeventscomp']).median()
        mean_pct_cevents_move_1 = 100.0*(fragstats_df['cmove_1']/fragstats_df['numcalledeventscomp']).mean()
        std_pct_cevents_move_1 = 100.0*(fragstats_df['cmove_1']/fragstats_df['numcalledeventscomp']).std()
        minidx_pct_cevents_move_1 = (fragstats_df['cmove_1']/fragstats_df['numcalledeventscomp']).idxmin()
        maxidx_pct_cevents_move_1 = (fragstats_df['cmove_1']/fragstats_df['numcalledeventscomp']).idxmax()
        min_pct_cevents_move_1 = 100.0*(fragstats_df['cmove_1']/fragstats_df['numcalledeventscomp']).min()
        max_pct_cevents_move_1 = 100.0*(fragstats_df['cmove_1']/fragstats_df['numcalledeventscomp']).max()
        
        med_pct_cevents_move_2 = 100.0*(fragstats_df['cmove_2']/fragstats_df['numcalledeventscomp']).median()
        mean_pct_cevents_move_2 = 100.0*(fragstats_df['cmove_2']/fragstats_df['numcalledeventscomp']).mean()
        std_pct_cevents_move_2 = 100.0*(fragstats_df['cmove_2']/fragstats_df['numcalledeventscomp']).std()
        minidx_pct_cevents_move_2 = (fragstats_df['cmove_2']/fragstats_df['numcalledeventscomp']).idxmin()
        maxidx_pct_cevents_move_2 = (fragstats_df['cmove_2']/fragstats_df['numcalledeventscomp']).idxmax()
        min_pct_cevents_move_2 = 100.0*(fragstats_df['cmove_2']/fragstats_df['numcalledeventscomp']).min()
        max_pct_cevents_move_2 = 100.0*(fragstats_df['cmove_2']/fragstats_df['numcalledeventscomp']).max()
        
        med_pct_cevents_move_3 = 100.0*(fragstats_df['cmove_3']/fragstats_df['numcalledeventscomp']).median()
        mean_pct_cevents_move_3 = 100.0*(fragstats_df['cmove_3']/fragstats_df['numcalledeventscomp']).mean()
        std_pct_cevents_move_3 = 100.0*(fragstats_df['cmove_3']/fragstats_df['numcalledeventscomp']).std()
        minidx_pct_cevents_move_3 = (fragstats_df['cmove_3']/fragstats_df['numcalledeventscomp']).idxmin()
        maxidx_pct_cevents_move_3 = (fragstats_df['cmove_3']/fragstats_df['numcalledeventscomp']).idxmax()
        min_pct_cevents_move_3 = 100.0*(fragstats_df['cmove_3']/fragstats_df['numcalledeventscomp']).min()
        max_pct_cevents_move_3 = 100.0*(fragstats_df['cmove_3']/fragstats_df['numcalledeventscomp']).max()
        
        med_pct_cevents_move_4 = 100.0*(fragstats_df['cmove_4']/fragstats_df['numcalledeventscomp']).median()
        mean_pct_cevents_move_4 = 100.0*(fragstats_df['cmove_4']/fragstats_df['numcalledeventscomp']).mean()
        std_pct_cevents_move_4 = 100.0*(fragstats_df['cmove_4']/fragstats_df['numcalledeventscomp']).std()
        minidx_pct_cevents_move_4 = (fragstats_df['cmove_4']/fragstats_df['numcalledeventscomp']).idxmin()
        maxidx_pct_cevents_move_4 = (fragstats_df['cmove_4']/fragstats_df['numcalledeventscomp']).idxmax()
        min_pct_cevents_move_4 = 100.0*(fragstats_df['cmove_4']/fragstats_df['numcalledeventscomp']).min()
        max_pct_cevents_move_4 = 100.0*(fragstats_df['cmove_4']/fragstats_df['numcalledeventscomp']).max()
        
        med_pct_cevents_move_5 = 100.0*(fragstats_df['cmove_5']/fragstats_df['numcalledeventscomp']).median()
        mean_pct_cevents_move_5 = 100.0*(fragstats_df['cmove_5']/fragstats_df['numcalledeventscomp']).mean()
        std_pct_cevents_move_5 = 100.0*(fragstats_df['cmove_5']/fragstats_df['numcalledeventscomp']).std()
        minidx_pct_cevents_move_5 = (fragstats_df['cmove_5']/fragstats_df['numcalledeventscomp']).idxmin()
        maxidx_pct_cevents_move_5 = (fragstats_df['cmove_5']/fragstats_df['numcalledeventscomp']).idxmax()
        min_pct_cevents_move_5 = 100.0*(fragstats_df['cmove_5']/fragstats_df['numcalledeventscomp']).min()
        max_pct_cevents_move_5 = 100.0*(fragstats_df['cmove_5']/fragstats_df['numcalledeventscomp']).max()

        print ("\t").join(["metric", "median", "mean", "std_dev", "min", "max"])
        print ("\t").join([str(e) for e in ["template_0_moves", med_pct_tevents_move_0, mean_pct_tevents_move_0, std_pct_tevents_move_0, min_pct_tevents_move_0, max_pct_tevents_move_0]])
        print ("\t").join([str(e) for e in ["template_1_moves", med_pct_tevents_move_1, mean_pct_tevents_move_1, std_pct_tevents_move_1, min_pct_tevents_move_1, max_pct_tevents_move_1]])
        print ("\t").join([str(e) for e in ["template_2_moves", med_pct_tevents_move_2, mean_pct_tevents_move_2, std_pct_tevents_move_2, min_pct_tevents_move_2, max_pct_tevents_move_2]])
        print ("\t").join([str(e) for e in ["template_3_moves", med_pct_tevents_move_3, mean_pct_tevents_move_3, std_pct_tevents_move_3, min_pct_tevents_move_3, max_pct_tevents_move_3]])
        print ("\t").join([str(e) for e in ["template_4_moves", med_pct_tevents_move_4, mean_pct_tevents_move_4, std_pct_tevents_move_4, min_pct_tevents_move_4, max_pct_tevents_move_4]])
        print ("\t").join([str(e) for e in ["template_5_moves", med_pct_tevents_move_5, mean_pct_tevents_move_5, std_pct_tevents_move_5, min_pct_tevents_move_5, max_pct_tevents_move_5]])
        print ("\t").join([str(e) for e in ["complement_0_moves", med_pct_cevents_move_0, mean_pct_cevents_move_0, std_pct_cevents_move_0, min_pct_cevents_move_0, max_pct_cevents_move_0]])
        print ("\t").join([str(e) for e in ["complement_1_moves", med_pct_cevents_move_1, mean_pct_cevents_move_1, std_pct_cevents_move_1, min_pct_cevents_move_1, max_pct_cevents_move_1]])
        print ("\t").join([str(e) for e in ["complement_2_moves", med_pct_cevents_move_2, mean_pct_cevents_move_2, std_pct_cevents_move_2, min_pct_cevents_move_2, max_pct_cevents_move_2]])
        print ("\t").join([str(e) for e in ["complement_3_moves", med_pct_cevents_move_3, mean_pct_cevents_move_3, std_pct_cevents_move_3, min_pct_cevents_move_3, max_pct_cevents_move_3]])
        print ("\t").join([str(e) for e in ["complement_4_moves", med_pct_cevents_move_4, mean_pct_cevents_move_4, std_pct_cevents_move_4, min_pct_cevents_move_4, max_pct_cevents_move_4]])
        print ("\t").join([str(e) for e in ["complement_5_moves", med_pct_cevents_move_5, mean_pct_cevents_move_5, std_pct_cevents_move_5, min_pct_cevents_move_5, max_pct_cevents_move_5]])
        print
        print minidx_pct_tevents_move_0, maxidx_pct_tevents_move_0
        minidx_pct_tevents_move_0, maxidx_pct_tevents_move_0 = int(minidx_pct_tevents_move_0), int(maxidx_pct_tevents_move_0)

        print ("\t").join(["#metric", "strand", "value", "numevents_from_molecule", "seq_len_2d", "seq_len_template", "seq_len_complement", "Q_2d", "Q_template", "Q_complement", "name"])
        print "min_pct_0_moves", "template", min_pct_tevents_move_0, fragstats_df['numevents'][minidx_pct_tevents_move_0], nonetodash(fragstats_df['seqlen2d'][minidx_pct_tevents_move_0]), nonetodash(fragstats_df['seqlentemp'][minidx_pct_tevents_move_0]), nonetodash(fragstats_df['seqlencomp'][minidx_pct_tevents_move_0]), nonetodash(fragstats_df['meanscore2d'][minidx_pct_tevents_move_0]), nonetodash(fragstats_df['meanscoretemp'][minidx_pct_tevents_move_0]), nonetodash(fragstats_df['meanscorecomp'][minidx_pct_tevents_move_0]), fragstats_df['name'][minidx_pct_tevents_move_0]
        print "max_pct_0_moves", "template", max_pct_tevents_move_0, fragstats_df['numevents'][maxidx_pct_tevents_move_0], nonetodash(fragstats_df['seqlen2d'][maxidx_pct_tevents_move_0]), nonetodash(fragstats_df['seqlentemp'][maxidx_pct_tevents_move_0]), nonetodash(fragstats_df['seqlencomp'][maxidx_pct_tevents_move_0]), nonetodash(fragstats_df['meanscore2d'][maxidx_pct_tevents_move_0]), nonetodash(fragstats_df['meanscoretemp'][maxidx_pct_tevents_move_0]), nonetodash(fragstats_df['meanscorecomp'][maxidx_pct_tevents_move_0]), fragstats_df['name'][maxidx_pct_tevents_move_0]
        print "min_pct_1_moves", "template", min_pct_tevents_move_1, fragstats_df['numevents'][minidx_pct_tevents_move_1], nonetodash(fragstats_df['seqlen2d'][minidx_pct_tevents_move_1]), nonetodash(fragstats_df['seqlentemp'][minidx_pct_tevents_move_1]), nonetodash(fragstats_df['seqlencomp'][minidx_pct_tevents_move_1]), nonetodash(fragstats_df['meanscore2d'][minidx_pct_tevents_move_1]), nonetodash(fragstats_df['meanscoretemp'][minidx_pct_tevents_move_1]), nonetodash(fragstats_df['meanscorecomp'][minidx_pct_tevents_move_1]), fragstats_df['name'][minidx_pct_tevents_move_1]
        print "max_pct_1_moves", "template", max_pct_tevents_move_1, fragstats_df['numevents'][maxidx_pct_tevents_move_1], nonetodash(fragstats_df['seqlen2d'][maxidx_pct_tevents_move_1]), nonetodash(fragstats_df['seqlentemp'][maxidx_pct_tevents_move_1]), nonetodash(fragstats_df['seqlencomp'][maxidx_pct_tevents_move_1]), nonetodash(fragstats_df['meanscore2d'][maxidx_pct_tevents_move_1]), nonetodash(fragstats_df['meanscoretemp'][maxidx_pct_tevents_move_1]), nonetodash(fragstats_df['meanscorecomp'][maxidx_pct_tevents_move_1]), fragstats_df['name'][maxidx_pct_tevents_move_1]
        print "min_pct_2_moves", "template", min_pct_tevents_move_2, fragstats_df['numevents'][minidx_pct_tevents_move_2], nonetodash(fragstats_df['seqlen2d'][minidx_pct_tevents_move_2]), nonetodash(fragstats_df['seqlentemp'][minidx_pct_tevents_move_2]), nonetodash(fragstats_df['seqlencomp'][minidx_pct_tevents_move_2]), nonetodash(fragstats_df['meanscore2d'][minidx_pct_tevents_move_2]), nonetodash(fragstats_df['meanscoretemp'][minidx_pct_tevents_move_2]), nonetodash(fragstats_df['meanscorecomp'][minidx_pct_tevents_move_2]), fragstats_df['name'][minidx_pct_tevents_move_2]
        print "max_pct_2_moves", "template", max_pct_tevents_move_2, fragstats_df['numevents'][maxidx_pct_tevents_move_2], nonetodash(fragstats_df['seqlen2d'][maxidx_pct_tevents_move_2]), nonetodash(fragstats_df['seqlentemp'][maxidx_pct_tevents_move_2]), nonetodash(fragstats_df['seqlencomp'][maxidx_pct_tevents_move_2]), nonetodash(fragstats_df['meanscore2d'][maxidx_pct_tevents_move_2]), nonetodash(fragstats_df['meanscoretemp'][maxidx_pct_tevents_move_2]), nonetodash(fragstats_df['meanscorecomp'][maxidx_pct_tevents_move_2]), fragstats_df['name'][maxidx_pct_tevents_move_2]
        print "min_pct_3_moves", "template", min_pct_tevents_move_3, fragstats_df['numevents'][minidx_pct_tevents_move_3], nonetodash(fragstats_df['seqlen2d'][minidx_pct_tevents_move_3]), nonetodash(fragstats_df['seqlentemp'][minidx_pct_tevents_move_3]), nonetodash(fragstats_df['seqlencomp'][minidx_pct_tevents_move_3]), nonetodash(fragstats_df['meanscore2d'][minidx_pct_tevents_move_3]), nonetodash(fragstats_df['meanscoretemp'][minidx_pct_tevents_move_3]), nonetodash(fragstats_df['meanscorecomp'][minidx_pct_tevents_move_3]), fragstats_df['name'][minidx_pct_tevents_move_3]
        print "max_pct_3_moves", "template", max_pct_tevents_move_3, fragstats_df['numevents'][maxidx_pct_tevents_move_3], nonetodash(fragstats_df['seqlen2d'][maxidx_pct_tevents_move_3]), nonetodash(fragstats_df['seqlentemp'][maxidx_pct_tevents_move_3]), nonetodash(fragstats_df['seqlencomp'][maxidx_pct_tevents_move_3]), nonetodash(fragstats_df['meanscore2d'][maxidx_pct_tevents_move_3]), nonetodash(fragstats_df['meanscoretemp'][maxidx_pct_tevents_move_3]), nonetodash(fragstats_df['meanscorecomp'][maxidx_pct_tevents_move_3]), fragstats_df['name'][maxidx_pct_tevents_move_3]
        print "min_pct_4_moves", "template", min_pct_tevents_move_4, fragstats_df['numevents'][minidx_pct_tevents_move_4], nonetodash(fragstats_df['seqlen2d'][minidx_pct_tevents_move_4]), nonetodash(fragstats_df['seqlentemp'][minidx_pct_tevents_move_4]), nonetodash(fragstats_df['seqlencomp'][minidx_pct_tevents_move_4]), nonetodash(fragstats_df['meanscore2d'][minidx_pct_tevents_move_4]), nonetodash(fragstats_df['meanscoretemp'][minidx_pct_tevents_move_4]), nonetodash(fragstats_df['meanscorecomp'][minidx_pct_tevents_move_4]), fragstats_df['name'][minidx_pct_tevents_move_4]
        print "max_pct_4_moves", "template", max_pct_tevents_move_4, fragstats_df['numevents'][maxidx_pct_tevents_move_4], nonetodash(fragstats_df['seqlen2d'][maxidx_pct_tevents_move_4]), nonetodash(fragstats_df['seqlentemp'][maxidx_pct_tevents_move_4]), nonetodash(fragstats_df['seqlencomp'][maxidx_pct_tevents_move_4]), nonetodash(fragstats_df['meanscore2d'][maxidx_pct_tevents_move_4]), nonetodash(fragstats_df['meanscoretemp'][maxidx_pct_tevents_move_4]), nonetodash(fragstats_df['meanscorecomp'][maxidx_pct_tevents_move_4]), fragstats_df['name'][maxidx_pct_tevents_move_4]
        print "min_pct_5_moves", "template", min_pct_tevents_move_5, fragstats_df['numevents'][minidx_pct_tevents_move_5], nonetodash(fragstats_df['seqlen2d'][minidx_pct_tevents_move_5]), nonetodash(fragstats_df['seqlentemp'][minidx_pct_tevents_move_5]), nonetodash(fragstats_df['seqlencomp'][minidx_pct_tevents_move_5]), nonetodash(fragstats_df['meanscore2d'][minidx_pct_tevents_move_5]), nonetodash(fragstats_df['meanscoretemp'][minidx_pct_tevents_move_5]), nonetodash(fragstats_df['meanscorecomp'][minidx_pct_tevents_move_5]), fragstats_df['name'][minidx_pct_tevents_move_5]
        print "max_pct_5_moves", "template", max_pct_tevents_move_5, fragstats_df['numevents'][maxidx_pct_tevents_move_5], nonetodash(fragstats_df['seqlen2d'][maxidx_pct_tevents_move_5]), nonetodash(fragstats_df['seqlentemp'][maxidx_pct_tevents_move_5]), nonetodash(fragstats_df['seqlencomp'][maxidx_pct_tevents_move_5]), nonetodash(fragstats_df['meanscore2d'][maxidx_pct_tevents_move_5]), nonetodash(fragstats_df['meanscoretemp'][maxidx_pct_tevents_move_5]), nonetodash(fragstats_df['meanscorecomp'][maxidx_pct_tevents_move_5]), fragstats_df['name'][maxidx_pct_tevents_move_5]

        print "min_pct_0_moves", "complement", min_pct_cevents_move_0, fragstats_df['numevents'][minidx_pct_cevents_move_0], nonetodash(fragstats_df['seqlen2d'][minidx_pct_cevents_move_0]), nonetodash(fragstats_df['seqlentemp'][minidx_pct_cevents_move_0]), nonetodash(fragstats_df['seqlencomp'][minidx_pct_cevents_move_0]), nonetodash(fragstats_df['meanscore2d'][minidx_pct_cevents_move_0]), nonetodash(fragstats_df['meanscoretemp'][minidx_pct_cevents_move_0]), nonetodash(fragstats_df['meanscorecomp'][minidx_pct_cevents_move_0]), fragstats_df['name'][minidx_pct_cevents_move_0]
        print "max_pct_0_moves", "complement", max_pct_cevents_move_0, fragstats_df['numevents'][maxidx_pct_cevents_move_0], nonetodash(fragstats_df['seqlen2d'][maxidx_pct_cevents_move_0]), nonetodash(fragstats_df['seqlentemp'][maxidx_pct_cevents_move_0]), nonetodash(fragstats_df['seqlencomp'][maxidx_pct_cevents_move_0]), nonetodash(fragstats_df['meanscore2d'][maxidx_pct_cevents_move_0]), nonetodash(fragstats_df['meanscoretemp'][maxidx_pct_cevents_move_0]), nonetodash(fragstats_df['meanscorecomp'][maxidx_pct_cevents_move_0]), fragstats_df['name'][maxidx_pct_cevents_move_0]
        print "min_pct_1_moves", "complement", min_pct_cevents_move_1, fragstats_df['numevents'][minidx_pct_cevents_move_1], nonetodash(fragstats_df['seqlen2d'][minidx_pct_cevents_move_1]), nonetodash(fragstats_df['seqlentemp'][minidx_pct_cevents_move_1]), nonetodash(fragstats_df['seqlencomp'][minidx_pct_cevents_move_1]), nonetodash(fragstats_df['meanscore2d'][minidx_pct_cevents_move_1]), nonetodash(fragstats_df['meanscoretemp'][minidx_pct_cevents_move_1]), nonetodash(fragstats_df['meanscorecomp'][minidx_pct_cevents_move_1]), fragstats_df['name'][minidx_pct_cevents_move_1]
        print "max_pct_1_moves", "complement", max_pct_cevents_move_1, fragstats_df['numevents'][maxidx_pct_cevents_move_1], nonetodash(fragstats_df['seqlen2d'][maxidx_pct_cevents_move_1]), nonetodash(fragstats_df['seqlentemp'][maxidx_pct_cevents_move_1]), nonetodash(fragstats_df['seqlencomp'][maxidx_pct_cevents_move_1]), nonetodash(fragstats_df['meanscore2d'][maxidx_pct_cevents_move_1]), nonetodash(fragstats_df['meanscoretemp'][maxidx_pct_cevents_move_1]), nonetodash(fragstats_df['meanscorecomp'][maxidx_pct_cevents_move_1]), fragstats_df['name'][maxidx_pct_cevents_move_1]
        print "min_pct_2_moves", "complement", min_pct_cevents_move_2, fragstats_df['numevents'][minidx_pct_cevents_move_2], nonetodash(fragstats_df['seqlen2d'][minidx_pct_cevents_move_2]), nonetodash(fragstats_df['seqlentemp'][minidx_pct_cevents_move_2]), nonetodash(fragstats_df['seqlencomp'][minidx_pct_cevents_move_2]), nonetodash(fragstats_df['meanscore2d'][minidx_pct_cevents_move_2]), nonetodash(fragstats_df['meanscoretemp'][minidx_pct_cevents_move_2]), nonetodash(fragstats_df['meanscorecomp'][minidx_pct_cevents_move_2]), fragstats_df['name'][minidx_pct_cevents_move_2]
        print "max_pct_2_moves", "complement", max_pct_cevents_move_2, fragstats_df['numevents'][maxidx_pct_cevents_move_2], nonetodash(fragstats_df['seqlen2d'][maxidx_pct_cevents_move_2]), nonetodash(fragstats_df['seqlentemp'][maxidx_pct_cevents_move_2]), nonetodash(fragstats_df['seqlencomp'][maxidx_pct_cevents_move_2]), nonetodash(fragstats_df['meanscore2d'][maxidx_pct_cevents_move_2]), nonetodash(fragstats_df['meanscoretemp'][maxidx_pct_cevents_move_2]), nonetodash(fragstats_df['meanscorecomp'][maxidx_pct_cevents_move_2]), fragstats_df['name'][maxidx_pct_cevents_move_2]
        print "min_pct_3_moves", "complement", min_pct_cevents_move_3, fragstats_df['numevents'][minidx_pct_cevents_move_3], nonetodash(fragstats_df['seqlen2d'][minidx_pct_cevents_move_3]), nonetodash(fragstats_df['seqlentemp'][minidx_pct_cevents_move_3]), nonetodash(fragstats_df['seqlencomp'][minidx_pct_cevents_move_3]), nonetodash(fragstats_df['meanscore2d'][minidx_pct_cevents_move_3]), nonetodash(fragstats_df['meanscoretemp'][minidx_pct_cevents_move_3]), nonetodash(fragstats_df['meanscorecomp'][minidx_pct_cevents_move_3]), fragstats_df['name'][minidx_pct_cevents_move_3]
        print "max_pct_3_moves", "complement", max_pct_cevents_move_3, fragstats_df['numevents'][maxidx_pct_cevents_move_3], nonetodash(fragstats_df['seqlen2d'][maxidx_pct_cevents_move_3]), nonetodash(fragstats_df['seqlentemp'][maxidx_pct_cevents_move_3]), nonetodash(fragstats_df['seqlencomp'][maxidx_pct_cevents_move_3]), nonetodash(fragstats_df['meanscore2d'][maxidx_pct_cevents_move_3]), nonetodash(fragstats_df['meanscoretemp'][maxidx_pct_cevents_move_3]), nonetodash(fragstats_df['meanscorecomp'][maxidx_pct_cevents_move_3]), fragstats_df['name'][maxidx_pct_cevents_move_3]
        print "min_pct_4_moves", "complement", min_pct_cevents_move_4, fragstats_df['numevents'][minidx_pct_cevents_move_4], nonetodash(fragstats_df['seqlen2d'][minidx_pct_cevents_move_4]), nonetodash(fragstats_df['seqlentemp'][minidx_pct_cevents_move_4]), nonetodash(fragstats_df['seqlencomp'][minidx_pct_cevents_move_4]), nonetodash(fragstats_df['meanscore2d'][minidx_pct_cevents_move_4]), nonetodash(fragstats_df['meanscoretemp'][minidx_pct_cevents_move_4]), nonetodash(fragstats_df['meanscorecomp'][minidx_pct_cevents_move_4]), fragstats_df['name'][minidx_pct_cevents_move_4]
        print "max_pct_4_moves", "complement", max_pct_cevents_move_4, fragstats_df['numevents'][maxidx_pct_cevents_move_4], nonetodash(fragstats_df['seqlen2d'][maxidx_pct_cevents_move_4]), nonetodash(fragstats_df['seqlentemp'][maxidx_pct_cevents_move_4]), nonetodash(fragstats_df['seqlencomp'][maxidx_pct_cevents_move_4]), nonetodash(fragstats_df['meanscore2d'][maxidx_pct_cevents_move_4]), nonetodash(fragstats_df['meanscoretemp'][maxidx_pct_cevents_move_4]), nonetodash(fragstats_df['meanscorecomp'][maxidx_pct_cevents_move_4]), fragstats_df['name'][maxidx_pct_cevents_move_4]
        print "min_pct_5_moves", "complement", min_pct_cevents_move_5, fragstats_df['numevents'][minidx_pct_cevents_move_5], nonetodash(fragstats_df['seqlen2d'][minidx_pct_cevents_move_5]), nonetodash(fragstats_df['seqlentemp'][minidx_pct_cevents_move_5]), nonetodash(fragstats_df['seqlencomp'][minidx_pct_cevents_move_5]), nonetodash(fragstats_df['meanscore2d'][minidx_pct_cevents_move_5]), nonetodash(fragstats_df['meanscoretemp'][minidx_pct_cevents_move_5]), nonetodash(fragstats_df['meanscorecomp'][minidx_pct_cevents_move_5]), fragstats_df['name'][minidx_pct_cevents_move_5]
        print "max_pct_5_moves", "complement", max_pct_cevents_move_5, fragstats_df['numevents'][maxidx_pct_cevents_move_5], nonetodash(fragstats_df['seqlen2d'][maxidx_pct_cevents_move_5]), nonetodash(fragstats_df['seqlentemp'][maxidx_pct_cevents_move_5]), nonetodash(fragstats_df['seqlencomp'][maxidx_pct_cevents_move_5]), nonetodash(fragstats_df['meanscore2d'][maxidx_pct_cevents_move_5]), nonetodash(fragstats_df['meanscoretemp'][maxidx_pct_cevents_move_5]), nonetodash(fragstats_df['meanscorecomp'][maxidx_pct_cevents_move_5]), fragstats_df['name'][maxidx_pct_cevents_move_5]
        
##    move tc ratios -- is there a pattern?
    if timecheck:
        n_time_errors = sum(fragstats_df['timeerror'])
        pct_time_errors = 100.0*n_time_errors/n_molecules




def nonetodash(x):
    if not x:
        return "-"
    elif np.isnan(x):
        return "-"
    else:
        return x

def run(parser, args):
    fragstats_df = make_fragstats_dataframe(args.fragfile, extensive=args.extensive, timecheck=args.checktime)
    summarize_fragstats(fragstats_df, extensive=args.extensive, timecheck=args.checktime)



