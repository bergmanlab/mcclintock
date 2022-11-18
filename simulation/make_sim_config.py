import sys
import os
import json
import argparse

def main():

    args = parse_args()
    te_count = len(args.family)

    config = {'families': {}, 'mcclintock': {}}
    # tsd length and targets for each family
    for i in range(te_count):
        config['families'][args.family[i]] = {
            "TSD": args.tsd,
            "targets": args.bed[i]
        }
    
    # Mcc path and methods
    config['mcclintock'] = {
        "path": args.mcc,
        "methods": "ngs_te_mapper,ngs_te_mapper2,relocate,relocate2,temp,temp2,retroseq,popoolationte,popoolationte2,te-locate,teflon,tebreak"
    }

    
    with open(args.out, "w") as conf:
        json.dump(config, conf, indent=4)



def parse_args():
    parser = argparse.ArgumentParser(prog='make_sim_config.py', description="Create config file in json format for McClintock evaluation.")

    ## required ##
    parser.add_argument("--family", type=str, nargs='+', help="List of TE families. Required.", required=True)
    parser.add_argument("--bed", type=str, nargs='+', help="List of candidate insertion region files in bed format for each TE family, respectively. Or one bed file for all TE families. Required.", required=True)
    parser.add_argument("--mcc", type=str, help="Path to local McClintock repo. Required.", required=True)
    parser.add_argument("--out", type=str, help="File name of the output json. Required.", required=True)

    ## optional ##
    parser.add_argument("--tsd", type=int, help="Integer for length of target site duplication in bp. Default is 5 bp.", default=5, required=False)
    
    args = parser.parse_args()

    # parse te families and bed
    if not len(args.family) == len(args.bed):
        if not len(args.bed) == 1:
            sys.exit("ERROR: One candidate insertion region file must be specified for each TE family. Or one bed file for all TE families. \n")
        else:
            bed = os.path.abspath(args.bed[0])
            args.bed.extend([bed] * (len(args.family) - 1))
    else:
        for i in range(len(args.family)):
            args.bed[i] = os.path.abspath(args.bed[i])

    # parse mcc install path
    args.mcc = os.path.abspath(args.mcc)
    
    # parse output path
    args.out = os.path.abspath(args.out)

    # parse tsd length
    if args.tsd is None:
        args.tsd = 5

    return args



if __name__ == "__main__":
    main()

