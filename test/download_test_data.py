import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/../")
import scripts.mccutils as mccutils

md5 = {
    "fq1": "852a68bc04f64ed21668f0d3942faf9b",
    "fq2": "77a455b155bbf55388267304616af6ff",
    "reference" : "8c5bc40e2c77ed5e3988a9d1c36534a3",
}

url = {
    "fq1": "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR800/SRR800842/SRR800842_1.fastq.gz",
    "fq2": "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR800/SRR800842/SRR800842_2.fastq.gz",
    "reference": "http://hgdownload.soe.ucsc.edu/goldenPath/sacCer2/bigZips/chromFa.tar.gz"
}

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    # download fastq1
    fq1 = script_dir+"/SRR800842_1.fastq.gz"
    if os.path.exists(fq1):
        print(fq1+" already exists...skipping...")
    else:
        download_success = mccutils.download(url['fq1'], fq1, md5=md5['fq1'], timeout=22000)
        if not download_success:
            print("Download of SRR800842_1.fastq.gz failed...please rerun script to try again")
            mccutils.remove(fq1)
            sys.exit(1)

    # download fastq2
    fq2 = script_dir+"/SRR800842_2.fastq.gz"
    if os.path.exists(fq2):
        print(fq2+" already exists...skipping...")
    else:
        download_success = mccutils.download(url['fq2'], fq2, md5=md5['fq2'], timeout=22000)
        if not download_success:
            print("Download of SRR800842_2.fastq.gz failed...please rerun script to try again")
            mccutils.remove(fq2)
            sys.exit(1)

    

        





if __name__ == "__main__":                
    main()