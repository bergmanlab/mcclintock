import os
import subprocess
import sys
import urllib.request
import hashlib
import socket
import shutil
import errno

class Ngs_te_mapper:
    def __init__(self):
        self.support = -1

class Temp:
    def __init__(self):
        self.classification = "None"
        self.support = -1
        self.frequency = -1
        self.junction1 = -1
        self.junction1Support = -1
        self.junction2 = -1
        self.junction2Support = -1
        self.fivePrimeSupport = -1
        self.threePrimeSupport = -1

class Telocate:
    def __init__(self):
        self.read_pair_support = -1

class Retroseq:
    def __init__(self):
        self.read_support = -1
        self.breakpoint_confidence = -1

class Relocate2:
    def __init__(self):
        self.left_support = -1
        self.right_support = -1
        self.left_junction = -1
        self.right_junction = -1

class Relocate:
    def __init__(self):
        self.left_support = -1
        self.right_support = -1

class Popoolationte2:
    def __init__(self):
        self.support_type = ""
        self.frequency = -1
        self.added = False

class Popoolationte:
    def __init__(self):
        self.support_type = ""
        self.f_read_support = 0
        self.r_read_support = 0
        self.f_read_support_percent = 0
        self.r_read_support_percent = 0

class Insertion:
    def __init__(self):
        self.chromosome = "None"
        self.start = -1
        self.end = -1
        self.name = "None"
        self.type = "None"
        self.strand = "."
        self.family = ""
        self.popoolationte = Popoolationte()
        self.popoolationte2 = Popoolationte2()
        self.relocate = Relocate()
        self.relocate2 = Relocate2()
        self.retroseq = Retroseq()
        self.telocate = Telocate()
        self.temp = Temp()
        self.ngs_te_mapper = Ngs_te_mapper()


def mkdir(indir, log=None):
    if os.path.isdir(indir) == False:
        os.mkdir(indir)
    # else:
    #     msg = "cannot make dir:"+indir+" dir exists...skipping....\n"
    #     sys.stderr.write(msg)
    #     writelog(log, msg)

def get_abs_path(in_file, log=None):
    if os.path.isfile(in_file):
        return os.path.abspath(in_file)
    else:
        msg = " ".join(["Cannot find file:",in_file,"exiting....\n"])
        sys.stderr.write(msg)
        writelog(log, msg)
        sys.exit(1)

def get_base_name(path, fastq=False):
    no_path = os.path.basename(path)
    no_ext = no_path.split(".")[0]
    
    if fastq == True:
        no_ext = no_ext.replace("_1","")
        no_ext = no_ext.replace("_2","")

    return no_ext


def run_command_stdout(cmd_list, out_file, log=None):
    msg = ""
    if log is None:
        try:
            # print(" ".join(cmd_list)+" > "+out_file)
            out = open(out_file,"w")
            subprocess.check_call(cmd_list, stdout=out)
            out.close()
        except subprocess.CalledProcessError as e:
            if e.output is not None:
                msg = str(e.output)+"\n"
            if e.stderr is not None:
                msg += str(e.stderr)+"\n"
            cmd_string = " ".join(cmd_list)
            msg += msg + cmd_string + "\n"
            sys.stderr.write(msg)
            sys.exit(1)
    
    else:
        try:
            out_log = open(log,"a")
            out_log.write(" ".join(cmd_list)+" > "+out_file+"\n")
            out = open(out_file,"w")
            subprocess.check_call(cmd_list, stdout=out, stderr=out_log)
            out.close()
            out_log.close()

        except subprocess.CalledProcessError as e:
            if e.output is not None:
                msg = str(e.output)+"\n"
            if e.stderr is not None:
                msg += str(e.stderr)+"\n"
            cmd_string = " ".join(cmd_list)
            msg += msg + cmd_string + "\n"
            writelog(log, msg)
            sys.stderr.write(msg)
            sys.exit(1)


def run_command(cmd_list, log=None):
    msg = ""
    if log is None:
        try:
            # print(" ".join(cmd_list))
            subprocess.check_call(cmd_list)
        except subprocess.CalledProcessError as e:
            if e.output is not None:
                msg = str(e.output)+"\n"
            if e.stderr is not None:
                msg += str(e.stderr)+"\n"
            cmd_string = " ".join(cmd_list)
            msg += msg + cmd_string + "\n"
            sys.stderr.write(msg)
            sys.exit(1)
    
    else:
        try:
            out = open(log,"a")
            out.write(" ".join(cmd_list)+"\n")
            subprocess.check_call(cmd_list, stdout=out, stderr=out)
            out.close()

        except subprocess.CalledProcessError as e:
            if e.output is not None:
                msg = str(e.output)+"\n"
            if e.stderr is not None:
                msg += str(e.stderr)+"\n"
            cmd_string = " ".join(cmd_list)
            msg += msg + cmd_string + "\n"
            writelog(log, msg)
            sys.stderr.write(msg)
            sys.exit(1)
    

def writelog(log, msg):
    if log is not None:
        with open(log, "a") as out:
            out.write(msg)


def download(url, out_file, md5=None, timeout=15, _attempt=1, max_attempts=1):
    if _attempt > max_attempts:
        return False
    socket.setdefaulttimeout(timeout)
    print("downloading ", url, "to", out_file)
    try:
        urllib.request.urlretrieve(url, out_file)
    
    except KeyboardInterrupt:
        sys.exit(1)

    except:
        print(sys.exc_info())
        print("download failed...")
        return download(url, out_file, md5=md5, timeout=timeout, _attempt=_attempt+1, max_attempts=max_attempts)

    print("download complete")

    if md5 == None:
        return True
    else:
        print("checking md5 of ", out_file)
        with open(out_file,"rb") as out:
            data = out.read()

            this_md5 = hashlib.md5(data).hexdigest()

            if this_md5 != md5:
                print("MD5 of", out_file, " : ", this_md5, "does not match currect MD5 hash: ", md5)
                return download(url, out_file, md5=md5, timeout=timeout, _attempt=_attempt+1, max_attempts=max_attempts)
            else:
                print("MD5 hash of ", out_file, "matches expected")
                return True


def remove(infile):
    if os.path.exists(infile):
        try:
            if os.path.isfile(infile):
                os.remove(infile)
            elif os.path.isdir(infile):
                shutil.rmtree(infile)
        except OSError as e:
            print("Error: %s : %s" % (infile, e.strerror))


def get_median_insert_size(infile):
    median_insert_size = 0
    with open(infile,"r") as inf:
        for line in inf:
            insert = line.split("=")[1]
            insert = insert.replace("\n","")
            median_insert_size = int(float(insert))
    
    return median_insert_size



def check_file_exists(infile):
    if os.path.exists(infile):
        return True
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), infile)


def replace_special_chars(string, encode="_"):
    special_chars = [";","&","(",")","|","*","?","[","]","~","{","}","<","!","^",'"',"'","\\","$","/","+","-"]
    for char in special_chars:
        string = string.replace(char,encode)
    
    return string


def replace_special_chars_fasta(fasta, new_fasta):
    with open(fasta,"r") as fa:
        with open(new_fasta, "w") as out:
            for line in fa:
                if ">" in line:
                    line = replace_special_chars(line)
                
                out.write(line)
    
    return new_fasta

def replace_special_chars_taxonomy(infile, outfile):
    with open(infile,"r") as i:
        with open(outfile, "w") as o:
            for line in i:
                line = replace_special_chars(line)
                o.write(line)
    
    return outfile


def estimate_read_length(fq, reads=10000):
    lengths = []
    with open(fq,"r") as f:
        for x, line in enumerate(f):
            if x%4 == 1:
                lengths.append(len(line.replace('\n',"")))
            
            if x >= reads:
                break
    
    length = sum(lengths)//len(lengths)

    return length


def log(step, msg, log=None):
    step = step.upper()
    max_step = 15
    if len(step) > max_step:
        step = step[:max_step]
    
    step_buffer = (max_step - len(step)) + 2
    step += " "*step_buffer

    max_msg = 103

    lines = msg
    lines = step + lines

    if log is not None:
        lines += " &> " + log
    
    print(lines)
    

def make_redundant_bed(insertions, sample_name, out_dir, method="popoolationte"):
    tmp_bed = out_dir+"/tmp.bed"

    insertion_dict = {}
    out_inserts = []
    for insert in insertions:
        insertion_dict[ "_".join([insert.chromosome, str(insert.start-1), str(insert.end), insert.name, "0", insert.strand])] = insert


    with open(tmp_bed, "w") as out:
        for insert in insertions:
            out_line = "\t".join([insert.chromosome, str(insert.start-1), str(insert.end), insert.name, "0", insert.strand])
            out.write(out_line+"\n")
    
    sorted_bed = out_dir+"/sorted.bed"
    command = ["bedtools", "sort", "-i", tmp_bed]
    run_command_stdout(command, sorted_bed)

    redundant_bed = out_dir+"/"+sample_name+"_"+method+"_redundant.bed"
    with open(redundant_bed, "w") as outbed:
        header = 'track name="'+sample_name+'_'+method+'" description="'+sample_name+'_'+method+'"\n'
        outbed.write(header)
        with open(sorted_bed, "r") as inbed:
            for x, line in enumerate(inbed):

                # outputs inserts in sorted order with unique number added to name
                key = line.replace("\t","_")
                key = key.replace("\n","")
                insert = insertion_dict[key]
                insert.name += str(x+1)
                out_inserts.append(insert)

                # write to bed with unique number added to name
                split_line = line.split("\t")
                split_line[3] += str(x+1)
                line = "\t".join(split_line)
                outbed.write(line)
    
    remove(tmp_bed)
    remove(sorted_bed)

    return out_inserts

def make_nonredundant_bed(insertions, sample_name, out_dir, method="popoolationte"):
    uniq_inserts = {}

    for insert in insertions:
        key = "_".join([insert.chromosome, str(insert.end)])
        if key not in uniq_inserts.keys():
            uniq_inserts[key] = insert
        else:
            ## method specific way to determine which duplicate to keep
            if method == "popoolationte":
                if (uniq_inserts[key].popoolationte.f_read_support + uniq_inserts[key].popoolationte.r_read_support) >  (insert.popoolationte.f_read_support + insert.popoolationte.r_read_support):
                    uniq_inserts[key] = insert
            
            elif method == "popoolationte2":
                if (uniq_inserts[key].popoolationte2.frequency) >  (insert.popoolationte2.frequency):
                    uniq_inserts[key] = insert
            
            elif method == "relocate":
                if (uniq_inserts[key].relocate.left_support + uniq_inserts[key].relocate.right_support) < (insert.relocate.left_support + insert.relocate.right_support):
                    uniq_inserts[key] = insert
            
            elif method == "relocate2":
                if (uniq_inserts[key].left_support + uniq_inserts[key].right_support) < (insert.left_support + insert.right_support):
                    uniq_inserts[key] = insert
            
            elif method == "retroseq":
                if uniq_inserts[key].retroseq.read_support >  insert.retroseq.read_support:
                    uniq_inserts[key] = insert
            
            elif method == "te-locate":
                if uniq_inserts[key].telocate.read_pair_support >  insert.telocate.read_pair_support:
                    uniq_inserts[key] = insert
            
            elif method == "temp":
                if uniq_inserts[key].temp.support != "!" and insert.temp.support > uniq_inserts[key].temp.support:
                    uniq_inserts[key] = insert
            
            elif method == "ngs_te_mapper":
                if insert.ngs_te_mapper.support > uniq_inserts[key].ngs_te_mapper.support:
                    uniq_inserts[key] = insert
    
    tmp_bed = out_dir+"/tmp.bed"
    with open(tmp_bed, "w") as outbed:
        for key in uniq_inserts.keys():
            insert = uniq_inserts[key]
            out_line = "\t".join([insert.chromosome, str(insert.start-1), str(insert.end), insert.name, "0", insert.strand])
            outbed.write(out_line+"\n")
    
    sorted_bed = out_dir+"/sorted.bed"
    command = ["bedtools", "sort", "-i", tmp_bed]
    run_command_stdout(command, sorted_bed)

    nonredundant_bed = out_dir+"/"+sample_name+"_"+method+"_nonredundant.bed"
    with open(sorted_bed, "r") as inbed:
        with open(nonredundant_bed, "w") as outbed:
            header = 'track name="'+sample_name+'_'+method+'" description="'+sample_name+'_'+method+'"\n'
            outbed.write(header)
            for line in inbed:
                outbed.write(line)

    remove(tmp_bed)
    remove(sorted_bed)