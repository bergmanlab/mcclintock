import sys
import os
import json


def main():
    # usage: python make_new_sim_config.py $mccdir
    mcc2_path = os.path.abspath(sys.argv[1])
    cwd = os.getcwd()
    if not os.path.exists(cwd+"/config"):
        os.mkdir(cwd+"/config")

    # config file for tRNA
    config = {}
    config['families'] = {
        "TY1": {
            "TSD": 5,
            "targets": mcc2_path+"/auxiliary/simulation/yeast/200-195_tRNA_targets.bed"
        },
        "TY2": {
            "TSD": 5,
            "targets": mcc2_path+"/auxiliary/simulation/yeast/200-195_tRNA_targets.bed"
        },
        "TY3": {
            "TSD": 5,
            "targets": mcc2_path+"/auxiliary/simulation/yeast/17-12_tRNA_targets.bed"
        },
        "TY4": {
            "TSD": 5,
            "targets": mcc2_path+"/auxiliary/simulation/yeast/200-195_tRNA_targets.bed"
        }
    }

    config['mcclintock'] = {
        "path": mcc2_path,
        "methods": "ngs_te_mapper,ngs_te_mapper2,relocate,relocate2,temp,temp2,retroseq,popoolationte,popoolationte2,te-locate,teflon,tebreak,coverage"
    }

    with open(cwd+"/config/tRNA.json", "w") as conf:
        json.dump(config, conf, indent=4)

    # config file for non-TE
    config = {}
    config['families'] = {
        "TY1": {
            "TSD": 5,
            "targets": mcc2_path+"/auxiliary/simulation/yeast/non_TE_targets.bed"
        },
        "TY2": {
            "TSD": 5,
            "targets": mcc2_path+"/auxiliary/simulation/yeast/non_TE_targets.bed"
        },
        "TY3": {
            "TSD": 5,
            "targets": mcc2_path+"/auxiliary/simulation/yeast/non_TE_targets.bed"
        },
        "TY4": {
            "TSD": 5,
            "targets": mcc2_path+"/auxiliary/simulation/yeast/non_TE_targets.bed"
        }
    }

    config['mcclintock'] = {
        "path": mcc2_path,
        "methods": "ngs_te_mapper,ngs_te_mapper2,relocate,relocate2,temp,temp2,retroseq,popoolationte,popoolationte2,te-locate,teflon,tebreak,coverage"
    }

    with open(cwd+"/config/non_TE.json", "w") as conf:
        json.dump(config, conf, indent=4)


if __name__ == "__main__":
    main()
