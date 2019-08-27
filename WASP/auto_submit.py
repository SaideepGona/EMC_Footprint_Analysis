'''
Author : Saideep Gona

Script allowing for the quick submission of jobs on the cluster

example:

/project2/lbarreiro/programs/Python-3.7.3/bin/python3 auto_submit.py [pipeline file] [2:4] [sbatch script]

/project2/lbarreiro/programs/Python-3.7.3/bin/python3 
auto_submit.py 
decaluwe.txt.sh 
37:39
/project2/lbarreiro/users/Saideep/pipelines/job_submission/example_batch_scripts/star.sbatch

'''

import sys, os, argparse
from datetime import datetime

cwd = os.getcwd()
print("Python Interpreter: ", sys.executable)
try:
    os.mkdir(cwd+"/submitted_jobs")
except:
    pass

parser = argparse.ArgumentParser()
# parser.add_argument("input_file", help="Input pipeline commmands")
# parser.add_argument("job_lines", help="start:end,start:end,start:end")
# parser.add_argument("sbatch_root", help="sbatch script to be appended to")
args = parser.parse_args()

args.input_file = "wasp.txt.sh"

job_range = [259, 344]
per_job = 2
jobs = [str(str(x) + ":" + str(x + per_job - 1)) for x in range(job_range[0], job_range[1], per_job)]
args.job_lines = ",".join(jobs)
print(args.job_lines)
# sys.exit()
# args.job_lines = "1:10,11:20,21:30,31:40,41:50,51:60,61:70,71:80,81:90,91:100,101:110,111:120,121:130,131:140,141:150,151:160,161:170,171:172"

# args.sbatch_root = "/project2/lbarreiro/users/Saideep/pipelines/job_submission/example_batch_scripts/star.sbatch"
# args.sbatch_root = "/project2/lbarreiro/users/Saideep/pipelines/job_submission/example_batch_scripts/standard.sbatch"
args.sbatch_root = "/project2/lbarreiro/users/Saideep/pipelines/job_submission/example_batch_scripts/samtools.sbatch"

job_split = args.job_lines.split(",")

commands = []

with open(args.input_file, "r") as inp:
    for line in inp:
        commands.append(line.rstrip("\n"))

for job in job_split:

    submit_commands_inds = [int(x) for x in job.split(":")]
    submit_commands = "\n".join(commands[submit_commands_inds[0]-1:submit_commands_inds[1]])

    with open(args.sbatch_root, "r") as sr:
        sbatch = sr.read()

    sbatch += "\n"
    sbatch += submit_commands


    print(sbatch)
    dtime = str(datetime.now())
    dtime = dtime.replace(" ", ",")

    job_script = cwd + "/submitted_jobs/"+ dtime + ".sbatch"
    with open(job_script, "w") as sb:
        sb.write(sbatch)

    os.system("sbatch " + job_script)
