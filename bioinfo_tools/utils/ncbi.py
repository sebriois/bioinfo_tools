from urllib import request, parse
import os
import time
import json
import subprocess

from bioinfo_tools.utils.log import Log
from bioinfo_tools.utils.sge import SgeJob


def fetch_taxonomy(species_name):
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = parse.urlencode({
        'db': 'taxonomy',
        'term': species_name,
        'retmode': 'json'
    })
    uri = url + '?' + params
    
    with request.urlopen(uri) as r:
        ncbi_response = r.read()
    
    json_response = json.loads(ncbi_response.decode())
    
    ncbi_id = None
    ncbi_xml = None
    
    if int(json_response.get("esearchresult", {}).get("count", 0)) > 0:
        ncbi_id = json_response.get("esearchresult", {}).get("idlist", [])[0]
        if ncbi_id:
            url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=" + ncbi_id
            with request.urlopen(url) as r:
                ncbi_response = r.read()
            ncbi_xml = ncbi_response.decode()
    
    return ncbi_id, ncbi_xml


def makeblastdb(fasta_filepath, dbtype = 'prot'):
    cmd = [
        "makeblastdb",
        "-in", fasta_filepath,
        "-dbtype", dbtype
    ]
    subprocess.check_call(" ".join(cmd), shell = True)


class BlastCommand(Log):
    def __init__(self, *args, **kwargs):
        super().__init__(**kwargs)
        
        self.blast_prog = None
        self.sge_job = None
        
        for arg_name in ['query', 'db', 'out']:
            if arg_name not in kwargs:
                raise Exception("'%s' argument is required" % arg_name)

        self.args = kwargs
        if not os.path.exists(self.args['query']):
            raise Exception("query file not found: " + self.args['query'])

        if isinstance(self.args['db'], list):
            self.args['db'] = "'" + " ".join(self.args['db']) + "'"
    
    def run(self, job_name = None, use_sge = False, async = False):
        """
        run blast command, either locally or through qsub
        """
        blast_commandline = self._create_commandline()

        if use_sge:
            if not job_name:
                job_name = os.path.basename(self.args['out']).split('.')[0]
            self.sge_job = SgeJob()
            self.sge_job.submit(blast_commandline, job_name = '%s_%s' % (self.blast_prog, job_name), sync = not async)
            return self.sge_job.get_job_id()
        else:
            if async:
                # TODO: make this option available locally
                self.log("Can't run in async mode locally. Will run synchronously.")
            
            try:
                return subprocess.check_output(blast_commandline, shell = True)
            except subprocess.CalledProcessError as e:
                raise Exception(e.output)
    
    def run_async(self, job_name = None, use_sge = False):
        return self.run(job_name = job_name, use_sge = use_sge, async = True)
    
    def wait_for_output(self, max_checks = 10):
        """
        wait for output file to be writen on disk
        """
        while not os.path.exists(self.args['out']) and max_checks > 0:
            time.sleep(10)
            max_checks -= 1
    
        if not os.path.exists(self.args['out']):
            raise Exception("blast output does not exist: " + self.args['out'])

    def _create_commandline(self):
        if not self.blast_prog:
            raise Exception("blast_prog must be set")
    
        blast_params = []
        for key, value in self.args.items():
            if not key.startswith('-'):
                key = '-' + key
            blast_params.append(key + '=' + str(value))
        return self.blast_prog + " " + " ".join(blast_params)


class Blastp(BlastCommand):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.blast_prog = 'blastp'


class Blastx(BlastCommand):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.blast_prog = 'blastx'


class Blastn(BlastCommand):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.blast_prog = 'blastn'
