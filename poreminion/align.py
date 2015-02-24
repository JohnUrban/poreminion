
## Align:
## how many of each type of read aligns? %?
## how well does each type align?

def run(parser, args):
    if args.align_command == 'blasr':
        import blasr_poreminion as align_module
        print "BLASR"
    elif args.align_command == "fitting":
        import fitting_aln as align_module
        print "FITTING"
    elif args.align_command == "blastn":
        import blastn_poreminion as align_module
        print"BLASTN"
    elif args.align_command == "last":
        import last_poreminion as align_module
        print "LAST"

    # run the chosen sub-sub-module.
    align_module.run(parser, args)




