import os
import subprocess
import re
from datetime import datetime
import time

DEFAULT_SCRATCH_DIR = os.path.join(os.sep, "home", os.environ.get("USER"), "sge_logs")

MAX_WAIT = 120  # seconds


class SgeJob(object):
    def __init__(self, scratch_dir = DEFAULT_SCRATCH_DIR):
        self.scratch_dir = scratch_dir
        os.makedirs(self.scratch_dir, exist_ok = True)
        self.params = []
    
    def set_params(self, *args, **kwargs):
        self.params = []
        
        for param in args:
            if not param.startswith('-'):
                param = '-' + param
            self.params.append(param)
        
        for param, value in kwargs.items():
            if not param.startswith('-'):
                param = '-' + param
            self.params.append(param)
            self.params.append(value)
    
    def submit(self, command_line, job_name = 'NO_NAME', sync = True):
        if '-q' not in self.params:
            self.params.extend(['-q', 'all.q'])
        
        if '-V' not in self.params:
            self.params.append('-V')
        
        if '-N' not in self.params:
            self.params.extend(['-N', job_name])
        
        if '-e' not in self.params:
            self.params.extend(['-e', os.path.join(self.scratch_dir, job_name + ".err")])
        
        if '-o' not in self.params:
            self.params.extend(['-o', os.path.join(self.scratch_dir, job_name + ".out")])
        
        if '-b' not in self.params:
            self.params.extend(['-b', 'y'])
        
        qsub_command = "qsub -clear %s '%s'" % (" ".join(self.params), command_line.strip().replace("'", '"'))
        print(qsub_command)
        qsub_response = subprocess.check_output(qsub_command, shell = True)
        
        try:
            job_id = re.findall("Your job (\d+) ", qsub_response.decode())[0]
        except IndexError:
            raise Exception("something went wrong with the job submission: %s" % qsub_response)
        
        job = self.qstat(job_id)
        
        if sync:
            wait_for = 2  # seconds
            while job and job['state'] != 'Eqw':
                print("job ID %s (%s) - state: %s (next check in %ssec)" % (job_id, job_name, job['state'], wait_for))
                wait_for *= 2
                wait_for = min([wait_for, MAX_WAIT])
                time.sleep(wait_for)
                job = self.qstat(job_id)
        
        return job
    
    def qstat(self, job_id = None):
        qstat_response = subprocess.check_output("qstat", shell = True)
        
        start_reading = False
        jobs = {}
        
        for line in qstat_response.decode().splitlines():
            if line.startswith("----"):
                start_reading = True
                continue
            
            if not start_reading:
                continue
            
            columns = re.split("\s+", line.strip())
            jobs[columns[0]] = {
                'id': columns[0],
                'priority': columns[1],
                'name': columns[2],
                'state': columns[4],
                'submit_date': datetime.strptime("%s %s" % (columns[5], columns[6]), "%m/%d/%Y %H:%M:%S"),
                'queue': columns[6]
            }
        
        if job_id:
            return jobs.get(job_id, None)
        return jobs