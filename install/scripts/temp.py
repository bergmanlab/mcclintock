import sys
import os
sys.path.append(snakemake.config['paths']['mcc_path'])
import scripts.mccutils as mccutils
import subprocess

def main():
    print("installing TEMP...")

    install_path = snakemake.config['paths']['install']+"/tools/"

    download_success = mccutils.download(snakemake.params.url, snakemake.params.zipfile, md5=snakemake.params.md5)

    if not download_success:
        print("temp download failed... retrying...")
        download_success = mccutils.download(snakemake.params.url, snakemake.params.zipfile, md5=snakemake.params.md5)
        if not download_success:
            print("temp second download attempt failed... exiting...")
            print("try running --install with --clean for clean installation")
            sys.exit(1)

    command = ["unzip", snakemake.params.zipfile]
    mccutils.run_command(command, log=snakemake.params.log)

    command = ["mv", snakemake.config['paths']['install']+"TEMP-d2500b904e2020d6a1075347b398525ede5feae1", install_path]
    mccutils.run_command(command, log=snakemake.params.log)

    for f in os.listdir(install_path+"TEMP-d2500b904e2020d6a1075347b398525ede5feae1"):
        command = ["mv", install_path+"TEMP-d2500b904e2020d6a1075347b398525ede5feae1/"+f, install_path+"temp"]
        mccutils.run_command(command, log=snakemake.params.log)


    command = ["patch", "-i", snakemake.params.patch, install_path+"temp/scripts/TEMP_Absence.sh"]
    mccutils.run_command(command,log=snakemake.params.log) 


    output = subprocess.Popen(["which", "perl"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    perl_path = output.stdout.read()
    perl_path = perl_path.decode()
    perl_path = perl_path.replace("\n","")

    with open(install_path+"temp/scripts/TEMP_Absence.sh", "r") as script:
        with open(install_path+"tmp", "w") as tmp:
            for line in script:
                if "perl" in line and ".pl" in line:
                    line = line.replace("perl", perl_path)
                
                tmp.write(line)
    
    mccutils.run_command(["mv",install_path+"tmp", install_path+"temp/scripts/TEMP_Absence.sh"])

    with open(install_path+"temp/scripts/TEMP_Insertion.sh", "r") as script:
        with open(install_path+"tmp", "w") as tmp:
            for line in script:
                if "perl" in line and ".pl" in line:
                    line = line.replace("perl", perl_path)
                
                tmp.write(line)

    mccutils.run_command(["mv",install_path+"tmp", install_path+"temp/scripts/TEMP_Insertion.sh"])

    mccutils.remove(install_path+"TEMP-d2500b904e2020d6a1075347b398525ede5feae1")
    mccutils.remove(snakemake.params.zipfile)

    print("TEMP installation complete")

if __name__ == "__main__":                
    main()