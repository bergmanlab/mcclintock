import sys
import os
import subprocess
sys.path.append(snakemake.config['paths']['mcc_path'])
import scripts.mccutils as mccutils

def main():
    install_path = snakemake.config['paths']['install']+"/tools/"
    mccutils.remove(snakemake.params.zipfile)
    download_success = mccutils.download(snakemake.params.url, snakemake.params.zipfile, md5=snakemake.params.md5, max_attempts=3)

    if not download_success:
        print("relocaTE download failed... exiting...")
        sys.exit(1)

    mccutils.remove(snakemake.config['paths']['install']+"RelocaTE-ce3a2066e15f5c14e2887fdf8dce0485e1750e5b")
    command = ["unzip", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

    mccutils.remove(install_path+"RelocaTE-ce3a2066e15f5c14e2887fdf8dce0485e1750e5b")
    command = ["mv", snakemake.config['paths']['install']+"RelocaTE-ce3a2066e15f5c14e2887fdf8dce0485e1750e5b", install_path]
    mccutils.run_command(command, log=snakemake.params.log)

    mccutils.remove(install_path+"relocate")
    mccutils.mkdir(install_path+"relocate/")
    for f in os.listdir(install_path+"RelocaTE-ce3a2066e15f5c14e2887fdf8dce0485e1750e5b"):
        command = ["mv", install_path+"RelocaTE-ce3a2066e15f5c14e2887fdf8dce0485e1750e5b/"+f, install_path+"relocate/"]
        mccutils.run_command(command, log=snakemake.params.log)


    command = ["patch", "-i", snakemake.params.patch, install_path+"relocate/scripts/relocaTE_insertionFinder.pl"]
    mccutils.run_command(command, log=snakemake.params.log)

    mccutils.remove(install_path+"RelocaTE-ce3a2066e15f5c14e2887fdf8dce0485e1750e5b")
    mccutils.remove(snakemake.params.zipfile)

    output = subprocess.Popen(["which", "perl"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    perl_path = output.stdout.read()
    perl_path = perl_path.decode()

    for f in os.listdir(install_path+"relocate/scripts/"):
        if "pl" == f.split(".")[-1]:
            with open(install_path+"tmp","w") as tmp:
                with open(install_path+"relocate/scripts/"+f, "r") as script:
                    for line in script:
                        if "#!/usr/bin/perl" in line:
                            # line = "#!"+perl_path
                            line = "#!/usr/bin/env perl\n"
                        elif "defined @" in line:
                            line = line.replace("defined @", "@")
                        
                        elif "$scripts/" in line and "perl" not in line and "relocaTE.pl" in f:
                            line = line.replace("$scripts/", "perl $scripts/")
                        
                        tmp.write(line)
            
            mccutils.run_command(["mv", install_path+"tmp", install_path+"relocate/scripts/"+f])

    # write version to file
    with open(snakemake.config['paths']['install']+"/tools/relocate/version.log","w") as version:
        version.write(snakemake.params.md5)

if __name__ == "__main__":                
    main()