import os
import subprocess
import re
from datetime import datetime
import time
import xml.etree.ElementTree as ET

from bioinfo_tools.utils.log import Log

DEFAULT_SCRATCH_DIR = os.path.join(os.sep, os.environ.get("HOME"), "sge_logs")

MAX_WAIT = 120  # seconds


class SgeJob(Log):
    def __init__(self, scratch_dir = DEFAULT_SCRATCH_DIR, **kwargs):
        super().__init__(**kwargs)
        
        self.scratch_dir = scratch_dir
        os.makedirs(self.scratch_dir, exist_ok = True)
        self.params = []
        self._job_id = None
    
    def get_job_id(self):
        if not self._job_id:
            self.log("No command submitted yet, there's no job ID to return")
        return self._job_id
    
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
    
    def submit(self, command_line, job_name = 'NO_NAME', sync = False):
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
        self.log(qsub_command)
        qsub_response = subprocess.check_output(qsub_command, shell = True)
        
        try:
            self._job_id = re.findall("Your job (\d+) ", qsub_response.decode())[0]
        except IndexError:
            raise Exception("something went wrong with the job submission: %s" % qsub_response)
        
        job = self.qstat(self._job_id)
        
        if sync:
            wait_for = 2  # seconds
            while job and job['state'] != 'Eqw':
                self.log("job ID %s (%s) - state: %s (next check in %ssec)" % (self._job_id, job_name, job['state'], wait_for))
                wait_for *= 2
                wait_for = min([wait_for, MAX_WAIT])
                time.sleep(wait_for)
                job = self.qstat(self._job_id)
        
        return job
    
    def qstat(self, job_id = None):
        xml_string = subprocess.check_output("qstat -xml", shell = True)
        xml_obj = ET.fromstring(xml_string)
        jobs = list()
        
        for job_list in xml_obj.iter("job_list"):
            job = dict()
            for elem in job_list:
                job[elem.tag] = elem.text
            
            if job_id and job['JB_job_number'] != job_id:
                continue
            
            jobs.append(job)
        
        if job_id and len(jobs) == 1:
            return jobs[0]
        
        return jobs
