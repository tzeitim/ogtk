import time 
import hashlib
import subprocess
import os
import re
from rich.console import Console
from rich.live import Live
import subprocess

class SGE_CONF:
    def __init__(self, **dict):
        for k,v in dict.items():
            setattr(self, k, v)
            
class SGE_JOB:
    def __init__(self, job_template_path, wd, sge_conf, id=None, exp=None, console=None):
        self.id = id if id is not None else "test"
        self.job_template_path = job_template_path
        self.job_template = None
        self.sge_conf = sge_conf
        self.read_job_template()
        self.filled_job = self.job_template
        self.hash_size = 7
        self.console = console

        if self.console is not None:
            self.use_console = True

        self.live = Live(console=console) if console else None

        #if exp is not None:
        #    self.generate_hash(exp)
        #elif self.id != "test":
        #    self.generate_hash(self)
        
        self.wd = f'{wd}/{self.id}'
        
        self.fill_job_template({ 'id': f'{self.id}',
                                'stdout': f'{self.wd}/log.std', 
                                'stderr': f'{self.wd}/log.err',
                                'workdir': f'{self.wd}'
                               })
        
        self.find_placeholders()
        
    def set_use_console(self, value):
        self.use_console = value

    def read_job_template(self):
        with open(self.job_template_path, 'r') as file:
            self.job_template = file.read()
        
    def fill_job_template(self, var_map):
        for k, v in var_map.items():
            self.filled_job = self.filled_job.replace('{'+k.upper()+'}', v)

    def generate_hash(self, object):
        if self.id != "test":
            attributes_string = ";".join(f"{key}={value}" for key, value in sorted(vars(object).items()))
            hash_string = hashlib.md5(attributes_string.encode()).hexdigest()[0:self.hash_size]
            self.id = f"{self.id}_{hash_string}"

    def submit_job(self, files):
        sge_conf = self.sge_conf

        # Create workspace if it doesn't exist
        if not os.path.exists(self.wd):
            print(f"created {self.wd}")
            os.makedirs(self.wd)

        # Transfer files to workspace using rsync
        print(f"Transfer began", end='')
        for file in files:
            subprocess.run(["rsync", file, self.wd], check=True)
        print(f"...ended")

        # Write the filled job script to a file in the working directory
        job_script_path = os.path.join(self.wd, f"{self.id}.sh")

        with open(job_script_path, 'w') as job_file:
            job_file.write(self.filled_job)

        # Submit the job using SSH
        ssh_command = ["ssh", f"{self.sge_conf.user}@{self.sge_conf.host}", "qsub", job_script_path]
        ssh_result = subprocess.run(ssh_command, check=True, capture_output=True)
        job_submitted_id = ssh_result.stdout.decode().strip()

        print(job_submitted_id)

        return job_submitted_id


    def find_placeholders(self):
        # Regular expression to find placeholders in format {PLACEHOLDER}
        pattern = r'\{([^}]+)\}'
        self.placeholders = re.findall(pattern, self.filled_job)

    def print_placeholders(self):
        self.print_placeholders()
        ''' Helper function for interactive work to determine which are the currently missing placeholders to replace
        '''
        print(self.id)
        print('{')
        for placeholder in self.placeholders:
            print(f"    '{placeholder}': '', ")
        print('}')

    def custom_print(self, message, live_update=False):
        """
        Print the message using rich.Console if available, otherwise use standard print.
        Optionally, update the Live display with the message.

        Args:
            message (str): The message to be printed.
            live_update (bool): If True, update the Live display instead of a simple print.
        """
        
        if self.console and live_update and self.use_console:
            with self.live:
                self.live.update(message)
        elif self.console and self.use_console:
            self.console.print(message)
        else:
            print(message)

    def qstat(self, times=10, sleep=1):
        ssh_cmd = ["ssh", f"{self.sge_conf.user}@{self.sge_conf.host}", f"for i in {{1..{times}}}; do qstat; sleep {sleep}; done"]
        process = subprocess.Popen(ssh_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                self.custom_print(output.strip(), live_update=True)

        process.wait()

    def send_cmd(self, cmd):
        '''
        '''
        # Submit the job using SSH
        ssh_command = ["ssh", f"{self.sge_conf.user}@{self.sge_conf.host}"]

        for i in cmd.split():
            ssh_command.append(i)

        ssh_result = subprocess.run(ssh_command, check=True, capture_output=True)
        output = ssh_result.stdout.decode().strip()

        print(output)

