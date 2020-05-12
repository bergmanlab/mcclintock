import os
import subprocess
import sys
import urllib.request
import hashlib
import socket
import shutil

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
    if no_path.split(".")[-1] == ".gz":
        # removes .gz from end of file
        no_path = ".".join(no_path.split(".")[:-1])
    no_ext = os.path.splitext(no_path)[0]
    

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


def download(url, out_file, md5=None, timeout=15):
    socket.setdefaulttimeout(timeout)
    print("downloading ", url, "to", out_file)
    try:
        urllib.request.urlretrieve(url, out_file)
    
    except KeyboardInterrupt:
        sys.exit(1)

    except:
        print(sys.exc_info())
        print("download failed...")
        return False

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
                return False
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