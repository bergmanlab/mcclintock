import os
import subprocess
import sys
import urllib.request
import hashlib
import socket
import shutil
import errno
import statistics
from datetime import date


INVALID_SYMBOLS = [";","&","(",")","|","*","?","[","]","~","{","}","<","!","^",'"',"'","\\","$","/","+","-","#"," "]

class insertSizeError(Exception):
    def __init__(self, message):
        self.message = message
    pass


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

def get_base_name(path):
    no_path = os.path.basename(path)
    if no_path[-3:] == ".gz":
        no_path = no_path[:-3]
    no_ext = ".".join(no_path.split(".")[:-1])

    return no_ext

def is_empty_file(infile):
    if os.stat(infile).st_size == 0:
        return True
    else:
        return False

def run_command_stdout(cmd_list, out_file, log=None, fatal=False):
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
            if fatal:
                sys.exit(1)
            else:
                return False
    
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
            if fatal:
                sys.exit(1)
            else:
                return False
    
    return True


def run_command(cmd_list, log=None, fatal=False):
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
            if fatal:
                sys.exit(1)
            else:
                return False
    
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
            if fatal:
                sys.exit(1)
            else:
                return False
    
    return True

def writelog(log, msg):
    if log is not None:
        with open(log, "a") as out:
            out.write(msg)


def download(url, out_file, md5=None, timeout=60, _attempt=1, max_attempts=1):
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
    if os.path.exists(infile) or os.path.islink(infile):
        try:
            if os.path.isfile(infile) or os.path.islink(infile):
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

def calc_median_insert_size(insam):
    insert_sizes = []
    with open(insam,"r") as sam:
        for line in sam:
            split_line = line.split("\t")
            if len(split_line) >= 8:
                insert_size = int(split_line[8])
                if insert_size > 0:
                    insert_sizes.append(insert_size)
    
    if len(insert_sizes) < 1:
        raise insertSizeError("Can't calculate median insert size due to lack of valid insert size values")
        return 0

    insert_sizes.sort()
    median = statistics.median(insert_sizes)

    return median

def check_file_exists(infile):
    if os.path.exists(infile):
        return True
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), infile)

def check_status_file(infile):
    succeeded = True
    if os.path.exists(infile):
        with open(infile,"r") as inf:
            for line in inf:
                if "FAILED" in line:
                    succeeded = False
    
    return succeeded

def replace_special_chars(string, encode="_"):
    for char in INVALID_SYMBOLS:
        string = string.replace(char,encode)
    
    return string


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
