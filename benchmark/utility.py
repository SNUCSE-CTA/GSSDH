import os, platform, subprocess, re, git, sys
import psutil
import datetime
import tqdm
from subprocess import Popen, PIPE
import pandas as pd
import numpy as np

def get_processor_name():
    if platform.system() == "Windows":
        return platform.processor()
    elif platform.system() == "Darwin":
        os.environ['PATH'] = os.environ['PATH'] + os.pathsep + '/usr/sbin'
        command ="sysctl -n machdep.cpu.brand_string"
        return subprocess.check_output(command, shell=True).strip().decode("UTF-8")
    elif platform.system() == "Linux":
        command = "cat /proc/cpuinfo"
        all_info = subprocess.check_output(command, shell=True).decode().strip()
        for line in all_info.split("\n"):
            if "model name" in line:
                return re.sub( ".*model name.*:", "", line,1)
    return ""


def get_avail_ram():
    return psutil.virtual_memory().total // (1024 ** 2)

def get_system_info():
    return (f"""Experiments running on {platform.uname().version}
    System CPU : {get_processor_name()}
    System RAM : {get_avail_ram()} MB""")

def get_current_repository_ver():
    repo = git.Repo(search_parent_directories=True)
    sha = repo.head.object.hexsha
    return sha

def get_experiment_date():
    return str(datetime.datetime.today().date())



def execute_binary(args):
    cmd = ' '.join(args)
    print(cmd)
    process = Popen(cmd, shell=True, stdout=PIPE, stderr=sys.stderr)
    (std_output, std_error) = process.communicate()
    process.wait()
    rc = process.returncode
    return rc, std_output, std_error

def to_number(s):
    try:
        return int(s)
    except ValueError:
        return float(s)
