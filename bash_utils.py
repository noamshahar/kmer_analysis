# Implement bash commands using subprocess

import subprocess
import os
script_path = os.getcwd()

def write_cmd(command_str):
    # Input: command as string.
    # output: activate command in bash via subprocess
    
    ps = subprocess.getoutput(command_str)
    return ps
