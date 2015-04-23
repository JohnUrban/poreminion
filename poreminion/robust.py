# robust.py
## fragrobust
from fragstats import *

def is_robust(fragsize, seqlen2d, seqlentemp, seqlencomp, has2d, hascomp):
    ## is fragsize estimate robust to all seq lengths?
    upper = fragsize+fragsize*0.2
    lower = fragsize-fragsize*0.2
    score = 0
    if seqlentemp <= upper and seqlentemp >= lower:
        score += 1
    if hascomp:
        if seqlencomp <= upper and seqlencomp >= lower:
            score += 1
    else: ## if no comp, then only temp, and fragsize robust by default
        score += 1 ## no need to check seqlen2d since if present it is fragsize and if not, no penalty
    if score == 2:
        robust = 1
    else:
        robust = 0
    return robust

def get_robust_scores(fragstats_df):
    ## takes in fragstats dataframe, returns dataframe with "robust" in it as well.
    robust_scores = []
    for i in range(len(fragstats_df['name'])):
        fragsize = fragstats_df['fragsize'][i]
        seqlen2d = fragstats_df['seqlen2d'][i]
        seqlentemp = fragstats_df['seqlentemp'][i]
        seqlencomp = fragstats_df['seqlencomp'][i]
        has2d = fragstats_df['has2d'][i]
        hascomp = fragstats_df['hascomp'][i]
        robust_scores.append(is_robust(fragsize, seqlen2d, seqlentemp, seqlencomp, has2d, hascomp))
    fragstats_df['robust'] = robust_scores
    return fragstats_df



def robust_message():
    print
    print "Robustness of molecule size estimation to other read sizes from same molecule is"
    print "not to be confused with goodness of the estimate since the best estimate is chosen in each situation."
    print "Molecule size is estimated as 2d read length if 2d read present, template read length if template only, longer of template and complement if both are present without a 2D."
    print "The molecule size estimate is considered robust when all read lengths present for a molecule are within +/- 20% the molecule size estimate."
    print "This means that Template-only reads get vacuously 100% robust estimates."
    print "Moreover, when the template:complement ratio is too large or too small, no 2d base-calling is performed leaving only 1D template and complement reads."
    print "It is typically these molecules that are not considered robust since the template and complement are such different sizes."
    print "Assuming hairpin parsing occurred correctly, the longer read in this situation is the best estimate of the molecule size."
    print "Larger template reads occur, for example, when there is a nick in the complement."
    print "Larger complements occur, for example, when the motor protein on the lead adapter falls off and the template zooms through the pore until the motor protein on the hairpin catches it." 
    print "Expecting the read lengths to be similar in this situation is not appropriate."
    print "Therefore, the estimates in these situations are vacuously not robust, despite being a good bet."
    print "Finally, molecules with 2D reads occur when the ratio of template to complement is relatively close to 1."
    print "Therefore, the robustness of those molecule size estimates can be thought of as vacuously robust since all read types are guaranteed to be of similar lengths."
    print "Given these caveats, robustness of molecule size estimate to other reads from the same molcule is somewhat dependent on the number of template only and molecules with 2D."
    print "Low robustness % of a dataset pretty much just means that there was a lot of molecules with template and complement, but no 2D"
    print "High robustness % of a dataset pretty much just means that there were few molecules with template and complement, but no 2D - although it cannot discern if there was a lot of template-only, a lot of 2D, or both."
    print "Looking at % of molecules with template only + % of molecules with 2D gives almost same number as % robust."
    print "In other words, this is a somewhat useless score."

    
def robust_summary(fragstats_df):
    print
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
    
    n_molecules_robust_to_molecule_size_estimate = sum(fragstats_df['robust'])
    pct_molecules_robust_to_molecule_size_estimate = 100.0*n_molecules_robust_to_molecule_size_estimate/n_molecules
    print "n_molecules_robust_to_molecule_size_estimate\t" + str(n_molecules_robust_to_molecule_size_estimate)
    print "pct_molecules_robust_to_molecule_size_estimate\t" + str(pct_molecules_robust_to_molecule_size_estimate)
    n_molecules_with_2d_robust_to_molecule_size_estimate = sum(fragstats_df['robust'][has2d])
    pct_molecules_with_2d_robust_to_molecule_size_estimate = 100.0*n_molecules_with_2d_robust_to_molecule_size_estimate/n_2d
    print "n_molecules_with_2d_robust_to_molecule_size_estimate\t" + str(n_molecules_with_2d_robust_to_molecule_size_estimate)
    print "pct_molecules_with_2d_robust_to_molecule_size_estimate\t" + str(pct_molecules_with_2d_robust_to_molecule_size_estimate)

    n_molecules_with_NO_2d_robust_to_molecule_size_estimate = sum(fragstats_df['robust'][no2d])
    pct_molecules_with_NO_2d_robust_to_molecule_size_estimate = 100.0*n_molecules_with_NO_2d_robust_to_molecule_size_estimate/n_no_2d
    print "n_molecules_with_NO_2d_robust_to_molecule_size_estimate\t" + str(n_molecules_with_NO_2d_robust_to_molecule_size_estimate)
    print "pct_molecules_with_NO_2d_robust_to_molecule_size_estimate\t" + str(pct_molecules_with_NO_2d_robust_to_molecule_size_estimate)

    n_molecules_with_template_only_robust_to_molecule_size_estimate = sum(fragstats_df['robust'][temponly])
    pct_molecules_with_template_only_robust_to_molecule_size_estimate = 100.0*n_molecules_with_template_only_robust_to_molecule_size_estimate/n_temp_only
    print "n_molecules_with_template_only_robust_to_molecule_size_estimate\t" + str(n_molecules_with_template_only_robust_to_molecule_size_estimate)
    print "pct_molecules_with_template_only_robust_to_molecule_size_estimate\t" + str(pct_molecules_with_template_only_robust_to_molecule_size_estimate)

    n_molecules_with_complement_robust_to_molecule_size_estimate = sum(fragstats_df['robust'][hascomp])
    pct_molecules_with_complement_robust_to_molecule_size_estimate = 100.0*n_molecules_with_complement_robust_to_molecule_size_estimate/n_comp
    print "n_molecules_with_complement_robust_to_molecule_size_estimate\t" + str(n_molecules_with_complement_robust_to_molecule_size_estimate)
    print "pct_molecules_with_complement_robust_to_molecule_size_estimate\t" + str(pct_molecules_with_complement_robust_to_molecule_size_estimate)

    n_molecules_with_complement_no2d_robust_to_molecule_size_estimate = sum(fragstats_df['robust'][hascomp][no2d])
    pct_molecules_with_complement_no2d_robust_to_molecule_size_estimate = 100.0*n_molecules_with_complement_no2d_robust_to_molecule_size_estimate/n_comp_no2d
    print "n_molecules_with_complement_no2d_robust_to_molecule_size_estimate\t" + str(n_molecules_with_complement_no2d_robust_to_molecule_size_estimate)
    print "pct_molecules_with_complement_no2d_robust_to_molecule_size_estimate\t" + str(pct_molecules_with_complement_no2d_robust_to_molecule_size_estimate)

    n_molecules_with_complement_has2d_robust_to_molecule_size_estimate = sum(fragstats_df['robust'][hascomp][has2d])
    pct_molecules_with_complement_has2d_robust_to_molecule_size_estimate = 100.0*n_molecules_with_complement_has2d_robust_to_molecule_size_estimate/n_comp_has2d
    print "n_molecules_with_complement_has2d_robust_to_molecule_size_estimate\t" + str(n_molecules_with_complement_has2d_robust_to_molecule_size_estimate) + "\t (Should be same as 'has 2d')"
    print "pct_molecules_with_complement_has2d_robust_to_molecule_size_estimate\t" + str(pct_molecules_with_complement_has2d_robust_to_molecule_size_estimate) + "\t (Should be same as 'has 2d')"
    print

def run(parser, args):
    if args.message:
        robust_message()
    if args.fragfile:
        fragstats_df = make_fragstats_dataframe(args.fragfile, extensive=False, timecheck=False)
        fragstats_df = get_robust_scores(fragstats_df)
        robust_summary(fragstats_df)

