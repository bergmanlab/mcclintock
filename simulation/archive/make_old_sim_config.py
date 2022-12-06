import sys
import os
import json

class Feature:
    def __init__(self):
        self.chrom = "?"
        self.start = -1
        self.end = -1
        self.strand = "."


def main():
    gff = sys.argv[1]
    mcc1_path = os.path.abspath(sys.argv[2])
    mcc2_path = os.path.abspath(sys.argv[3])
    cwd = os.getcwd()
    if not os.path.exists(cwd+"/config"):
        os.mkdir(cwd+"/config")
    
    if not os.path.exists(cwd+"/beds"):
        os.mkdir(cwd+"/beds")

    features = []
    with open(gff,"r") as g:
        for line in g:
            feature = Feature()
            split_line = line.split("\t")
            feature.chrom = split_line[0]
            feature.start = int(split_line[3])-1
            feature.end = int(split_line[4])
            feature.strand = split_line[6]

            features.append(feature)

    families = ["TY1", "TY2","TY3", "TY4"]
    fam_idx = 0
    feat_idx = 0
    for x in range(0,299):
        if feat_idx >= len(features):
            feat_idx = 0

        if fam_idx > 3:
            fam_idx = 0
        
        family = families[fam_idx]
        if fam_idx == 2:
            upstream_start = 17
            upstream_end = 12
        else:
            upstream_start = 200
            upstream_end = 195
        
        feature = features[feat_idx]

        if feature.strand == "+":
            new_start = feature.start - upstream_start
            new_end = new_start + 1

        else:
            new_start = feature.start + upstream_end + 1
            new_end = new_start + 1

        fam_idx += 1
        feat_idx += 1
        with open(cwd+"/beds/run_"+str(x)+".bed","w") as out:
            out_line = "\t".join([feature.chrom, str(new_start), str(new_end)])
            out.write(out_line+"\n")

        config = {}
        config['families'] = {
            family : {
                "TSD" : 5,
                "targets" : cwd+"/beds/run_"+str(x)+".bed"
            }
        }

        config['mcclintock'] = {
            "path" : mcc2_path,
            "v1_path": mcc1_path,
            "methods" : "ngs_te_mapper,relocate,temp,retroseq,popoolationte,te-locate"
        }

        with open(cwd+"/config/run_"+str(x)+".json","w") as conf:
            json.dump(config, conf, indent=4)

if __name__ == "__main__":
    main()


