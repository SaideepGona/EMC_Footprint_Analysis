'''
Author : Saideep Gona

This PYTHON3 script is intended for generating a basic paired-end RNA-seq 
pipeline on the midway cluster. This is a variation which performs a STAR
alignment with WASP tagging followed by filtering out reads with biased 
mapping. 

Basic instructions:
Copy this file into a "project directory" of your choosing. Put your raw data into
a directory called "data" next to this file. Should look like:

----project directory/
-------pipeline.py
-------data/
----------raw_data1
----------raw_data2
            .
            .
            .
----------raw_dataN

To run this script, navigate to the "project directory":

Then, run using:

/project2/lbarreiro/programs/Python-3.7.3/bin/python3 <pipeline>.py [output_filename]

This will output pipeline components into the [output_filename]

With these commands, the user can then submit jobs on compute nodes as they fancy. 
Note that any sinteractive/sbatch commands used should look something like these examples:

Example sinteractive command (readable):
    sinteractive
    --time=24:00:00        (Can be modified if runtime is longer than this)
    --cpus-per-task=8     (Number of processors for each task.
                           Default number of tasks is 1. 
                           Most pipeline software should be 
                           faster with more cpus)
    --mem-per-cpu=5000   (10GB memory per cpu)

Example sinteractive command (copy-paste):

    sinteractive --time=8:00:00 --cpus-per-task=8 --mem-per-cpu=5000

Dependencies include:
    TrimmingQC: Trimmomatic
    Aligner: STAR
    Assembler: Stringtie
'''

import sys, os, argparse, glob
from datetime import datetime

parser = argparse.ArgumentParser()
parser.add_argument("output_dir", help="Output text file for pipeline commands")
args = parser.parse_args()

# Boilerplate
cwd = os.getcwd()           # Store current working directory
print("Current Working Directory: ", cwd)
dtime = str(datetime.now())
print("Time of Running ", dtime)
print("Python Interpreter: ", sys.executable)

# GLOBALS
GENOME = "/project2/lbarreiro/REF/hg19/hs_ref_GRCh37.fasta"
SPLIT_GENOME_DIR = "/project2/lbarreiro/REF/hg19/split_genome"
CHROM_INFO = "/project2/lbarreiro/REF/hg19/hs_ref_GRCh37.fasta"
PHASED_DATA_DIR = "/project2/lbarreiro/users/katie/WASP/phasing_outputs/phased/"
STAR_INDEX = "/project2/lbarreiro/REF/hg19/STAR/"
ANNOTATION = "/project2/lbarreiro/REF/hg19/annotation-8-20-19/Homo_sapiens.GRCh37.87.gtf"
# PERSONAL_VCF = "/project2/lbarreiro/users/katie/admixture/vcf_files/all_my_samples_merged_chr.vcf"
PERSONAL_VCF = "/project2/lbarreiro/users/katie/admixture/vcf_files/all_my_samples_merged_chr.vcf"
ORIGINAL_SCRIPT_DIR = "/project2/lbarreiro/programs/sai_raw_pipelines"
EXPERIMENT = "Flu-Afr-Eur"     # Name of experiment for which data is sourced
USER = "Saideep"            # Put in your name

# Pre-run checks
if cwd == ORIGINAL_SCRIPT_DIR: # Check if someone is running in the original script directory rather than their own 
    print("DO NOT USE THIS VERSION OF THE SCRIPT. MAKE A COPY IN A DIFFERENT DIRECTORY")
    sys.exit()

# Directory containing all input files 
path_to_input_data = "/project2/lbarreiro/users/katie/RNA-seq_files/trimmed/paired/"

# Create directories for intermediate files
trimmomatic_output_dir_all = cwd + "/trimmed_data/"
os.system("mkdir" + trimmomatic_output_dir_all)
trimmomatic_output_dir_paired = trimmomatic_output_dir_all + "paired/"
os.system("mkdir" + trimmomatic_output_dir_paired)
trimmomatic_output_dir_unpaired = trimmomatic_output_dir_all + "unpaired/"
os.system("mkdir" + trimmomatic_output_dir_unpaired)

star_output_dir = cwd + "/star_output/"
os.system("mkdir " + star_output_dir)

filter_wasp_output_dir = cwd + "/filter_was_output/"
os.system("mkdir " + filter_wasp_output_dir)

h5_output_dir = cwd + "/h5_output/"
os.system("mkdir " + h5_output_dir)

snp2h5_output_dir = h5_output_dir + "snp2h5/"
os.system("mkdir " + snp2h5_output_dir)

fasta2h5_output_dir = h5_output_dir + "fasta2h5/"
os.system("mkdir " + fasta2h5_output_dir)

bam2h5_output_dir = h5_output_dir + "bam2h5/"
os.system("mkdir " + bam2h5_output_dir)

# Load dependencies, set run paths
os.system("module load java")
trimmomatic_path = "java -jar /project2/lbarreiro/programs/Trimmomatic-0.39/trimmomatic-0.39.jar"

os.system("module load STAR")
star_path = "STAR"

os.system("module load samtools")
samtools_sort_path = "samtools sort"
samtools_index_path = "samtools index"

filter_wasp_path = "python " + cwd + "/filter_WASP.py"

snp2h5_path = "/project2/lbarreiro/programs/WASP/snp2h5/snp2h5"

fasta2h5_path = "/project2/lbarreiro/programs/WASP/snp2h5/fasta2h5"

bam2h5_path = "python /project2/lbarreiro/programs/WASP/CHT/bam2h5"

get_target_regions_path = "python /project2/lbarreiro/programs/WASP/CHT/get_target_regions.py"

# STAR Indexing
# COMMENT OUT IF INDEXING ALREADY DONE 

# star_index = [
#     star_path,
#     "--runThreadN", "8",
#     "--runMode", "genomeGenerate",
#     "--genomeDir", "/project2/lbarreiro/REF/mouse_mm10/Mus_musculus/mouse_GRCm38_ENSEMBL/Mus_musculus/Ensembl/GRCm38/Sequence/STAR_test",
#     "--genomeFastaFiles", "/project2/lbarreiro/REF/mouse_mm10/Mus_musculus/mouse_GRCm38_ENSEMBL/Mus_musculus/Ensembl/GRCm38/Sequence/WholeGenomeFasta/genome.fa",
# ]
# os.system(" ".join(star_index))

class Pipeline():
    '''
    Add steps to the pipeline, and then run the pipeline when ready
    '''
    def __init__(self, name, dtime, user):
        self.name = name
        self.blocks = []

    def run_pipeline(self):
        print("STARTING PIPELINE")
        for step in self.steps:
            print("Running step: "+step.name)
            print("Command: "+ step.command)
            step.run_step()
            print("Completed step: "+step.name)

    def output_commands(self, write_to = False):

        step_count = 1
        block_count = 1

        if not write_to:
            for block in self.blocks:
                print("/"*100)
                print("(Block "+str(block_count)+") "+ block.name)
                print("/"*100)
                for step in block.steps:
                    print("*"*100)
                    print("(Step "+str(step_count)+") "+ step.name)
                    print("-"*100)
                    print(step.command)
                    print("*"*100)
                    step_count += 1
                block_count += 1
                print("/"*100)
        else:
            print("Writing to: " + write_to)
            writeable = ""
            shell_writeable = ""
            step_count = 1
            block_count = 1

            for block in self.blocks:

                writeable += "/"*100 + "\n"
                writeable += "(Block "+str(block_count)+") "+ block.name + "\n"
                writeable += "/"*100 + "\n"

                for step in block.steps:
                    writeable += "*"*100 + "\n"
                    writeable += "(Step "+str(step_count)+") "+ step.name + "\n"
                    writeable += "-"*100 + "\n"
                    writeable += step.command + "\n"
                    shell_writeable += step.command + "\n"
                    writeable += "*"*100 + "\n"
                    step_count += 1
                writeable += "/"*100 + "\n"
                block_count += 1
            with open(write_to, "w") as out:
                out.write(writeable)
            with open(write_to + ".sh", "w") as sh_out:
                sh_out.write(shell_writeable)

class Block():
    '''
    Unit of parallelization
    '''
    def __init__(self, name):
        self.name = name
        self.steps = []

class Step():
    '''
    Unit of individual jobs/commands
    '''
    def __init__(self, name, command):
        self.name = name
        self.command = command

    def run_step(self):
        os.system(self.command)

class Experiment():
    '''
    Describes the overall organization of the experiment, treatments, controls,
    etc. Should be modified on an experiment by experiment basis until standardized
    '''

    def __init__(self, name):
        self.name = name

    def parse_meta(self, meta_file):
        '''
        Parses a metadata file and organizes inputs into specific treatments/groupings
        '''
        

# Initialize pipeline and experiment object

cur_experiment = Experiment(EXPERIMENT)

rnaseq_pipeline = Pipeline("rnaseq", dtime, USER)

# Run Trimmomatic

trimming_block = Block("Trimming - Trimmomatic")

input_files = glob.glob(path_to_input_data + "/*")
input_files.sort()

paired_files = []

trimmomatic_command = [
    trimmomatic_path,
    "PE",
    "PLACEHOLDER INFORWARD(INDEX2)", "PLACEHOLDER INBACKWARD(INDEX3)",
    "PLACEHOLDER OUTFORWARDPAIRED(INDEX4)", "PLACEHOLDER OUTFORWARDUNPAIRED(INDEX5)",
    "PLACEHOLDER OUTREVERSEPAIRED(INDEX6)", "PLACEHOLDER OUTREVERSEUNPAIRED(INDEX7)",
    "ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36"
]

for x in range(int(len(input_files)/2)):

    current_reads = input_files[2*x: 2*x + 2]
    cur_sample = current_reads[0].split("/")[-1].rstrip(".fastq")
    paired_outputs = [(trimmomatic_output_dir_paired + x.split("/")[-1].rstrip(".fastq") + ".paired.fastq") for x in current_reads]

    for p in paired_outputs:
        paired_files.append(p)
    unpaired_outputs = [(trimmomatic_output_dir_unpaired + x.split("/")[-1].rstrip(".fastq") + ".unpaired.fastq") for x in current_reads]

    trimmomatic_command_cur = trimmomatic_command[:]
    trimmomatic_command_cur[2:4] = current_reads[0:2]
    trimmomatic_command_cur[4] = paired_outputs[0]
    trimmomatic_command_cur[5] = unpaired_outputs[0]
    trimmomatic_command_cur[6] = paired_outputs[1]
    trimmomatic_command_cur[7] = unpaired_outputs[1]
    trimmomatic_command_cur = " ".join(trimmomatic_command_cur)
    
    trimming_block.steps.append(Step("Trimmomatic_" + cur_sample, trimmomatic_command_cur))

# rnaseq_pipeline.blocks.append(trimming_block)

# Run STAR

alignment_block = Block("Alignment - STAR")

unique_sample_ids = []
input_files_trimmed = paired_files[:]
input_files_trimmed.sort()

star_out_files = []

star_command = [       
    star_path, 
    "--outFileNamePrefix", "PLACEHOLDER(INDEX2)",
    "--runThreadN", "8",
    "--genomeDir", STAR_INDEX,
    "--outSAMstrandField", "intronMotif",
    "--outSAMtype", "BAM SortedByCoordinate",
    "--outFilterIntronMotifs", "RemoveNoncanonical",
    "--varVCFfile", PERSONAL_VCF,
    "--waspOutputMode", "SAMtag",
    "--outSAMattributes", "vA vG",
    "--readFilesCommand", "zcat",
    "--readFilesIn"
]

for x in range(int(len(input_files_trimmed)/2)):

    current_reads = input_files[2*x: 2*x + 2]

    star_output_dir_cur_suf = current_reads[0].split("/")[-1].rstrip(".fastq")
    star_output_dir_cur = star_output_dir + star_output_dir_cur_suf + "/"
    os.system("mkdir " + star_output_dir_cur)
    unique_sample_ids.append(star_output_dir_cur_suf)

    star_command_cur = star_command[:]
    star_command_cur[2] = star_output_dir_cur
    star_out_files.append(star_output_dir_cur + "/Aligned.sortedByCoord.out.bam")
    star_command_cur = " ".join(star_command_cur) + " " + " ".join(current_reads)
    
    alignment_block.steps.append(Step("STAR_" + star_output_dir_cur_suf, star_command_cur))

rnaseq_pipeline.blocks.append(alignment_block)
# print(unique_sample_ids)
# Index alignment files (primary)

index_block = Block("Index alignments 1 - samtools index")

samtools_index_command = [
    samtools_index_path,
    "PLACEHOLDER INPUT(INDEX1)",
    "PLACEHOLDER OUTPUT(INDEX2)"
]

for fout in star_out_files:

    ind_file = fout + ".bai"
    # print(fout)
    samtools_index_command_cur = samtools_index_command[:]
    samtools_index_command_cur[1] = fout
    samtools_index_command_cur[2] = ind_file

    samtools_index_command_cur = " ".join(samtools_index_command_cur)

    index_block.steps.append(Step("samtools_index1", samtools_index_command_cur))

rnaseq_pipeline.blocks.append(index_block)

# Filter WASP alignments

filter_wasp_block = Block("Filter - filter_wasp.py")
filtered_wasp_outputs = []

filter_wasp_command = [
    filter_wasp_path,
    "-in", "PLACEHOLDER INPUT(INDEX2)",
    "-out", "PLACEHOLDER INPUT(INDEX4)"
]

for uid in unique_sample_ids:

    cur_file = star_output_dir + uid + "/Aligned.sortedByCoord.out.bam"
    filter_wasp_command_cur = filter_wasp_command[:]
    filter_wasp_command_cur[2] = cur_file
    filter_wasp_command_cur[4] = filter_wasp_output_dir + "/" + uid + ".filter_wasp.bam"
    filtered_wasp_outputs.append(filter_wasp_output_dir + "/" + uid + ".filter_wasp.bam")
    filter_wasp_command_cur = " ".join(filter_wasp_command_cur)

    filter_wasp_block.steps.append(Step("filter_wasp_" + uid, filter_wasp_command_cur))

rnaseq_pipeline.blocks.append(filter_wasp_block)

# Index alignment files

index_block = Block("Index alignments 2 - samtools index")

samtools_index_command = [
    samtools_index_path,
    "PLACEHOLDER INPUT(INDEX1)",
    "PLACEHOLDER OUTPUT(INDEX2)"
]

for fout in filtered_wasp_outputs:
    ind_file = fout + ".bai"

    samtools_index_command_cur = samtools_index_command[:]
    samtools_index_command_cur[1] = fout
    samtools_index_command_cur[2] = ind_file

    samtools_index_command_cur = " ".join(samtools_index_command_cur)

    index_block.steps.append(Step("samtools_index2", samtools_index_command_cur))

rnaseq_pipeline.blocks.append(index_block)

# Convert to HDF5

convert_hdf5_block = Block("Convert input files to hdf5 format")

# Convert SNP to HDF5

snp2h5_outputs = {
    "geno_probs": snp2h5_outputs + "/geno_probs.h5",
    "snp_index": snp2h5_outputs + "/snp_index.h5",
    "snp_tab": snp2h5_outputs + "/snp_tab.h5",
    "haplotype": snp2h5_outputs + "/haplotype.h5"
}

snp2h5_command = [
    snp2h5_path,
    "--chrom", CHROM_INFO,
    "--format", "impute",
    "--geno_prob", snp2h5_outputs["geno_probs"],
    "--snp_index", snp2h5_outputs["snp_index"],
    "--snp_tab", snp2h5_outputs["snp_tab"],
    "--haplotype", snp2h5_outputs["haplotype"],
    PHASED_DATA_DIR + "/chr*.impute2*gz"
]

snp2h5_command = " ".join(snp2h5_command)

convert_hdf5_block.steps.append(Step("snp2h5", snp2h5_command))

# Convert FASTA to HDF5

fasta2h5_command = [
    fasta2h5_path,
    "--chrom", CHROM_INFO,
    "--seq", fasta2h5_output_dir  + "/seq.h5",
    SPLIT_GENOME_DIR + "/*.fasta.gz"
]

fasta2h5_command = " ".join(fasta2h5_command)

convert_hdf5_block.steps.append(Step("fasta2h5", snp2h5_command))

# Convert BAM to HDF5

bam2h5_command = [
    bam2h5_path,
    "--snp_index", snp2h5_outputs["snp_index"],
    "--snp_tab", snp2h5_outputs["snp_tab"],
    "--haplotype", snp2h5_outputs["haplotype"],
    "--individual", "PLACEHOLDER INDIVIDUAL(INDEX4)",
    "--ref_as_counts", "PLACEHOLDER OUTPUTMATCHREF(INDEX5)",
    "--alt_as_counts", "PLACEHOLDER OUTPUTMATCHALT(INDEX6)",
    "--other_as_counts", "PLACEHOLDER OUTPUTMATCHOTHER(INDEX7)",
    "--read_counts", "PLACEHOLDER OUTPUT(INDEX8)",
    "PLACEHOLDER INPUTBAMFILE(INDEX9)"
]

for fout in filtered_wasp_outputs:

    bam2h5_command_cur = bam2h5_command[:]
    bam2h5_command_cur = " ".join(bam2h5_command_cur)

    cur_file = fout.split("/")[-1]

    bam2h5_command_cur[4] = cur_file
    current_out_dir = bam2h5_output_dir + "/" + cur_file 
    os.system("mkdir " + current_out_dir)

    bam2h5_command_cur[5] = current_out_dir+ "/ref_as_counts.h5"
    bam2h5_command_cur[6] = current_out_dir+ "/alt_as_counts.h5"
    bam2h5_command_cur[7] = current_out_dir+ "/other_as_counts.h5"
    bam2h5_command_cur[8] = current_out_dir+ "/read_counts.h5"
    bam2h5_command_cur[9] = fout

    convert_hdf5_block.steps.append(Step("bam2h5", bam2h5_command_cur))

rnaseq_pipeline.blocks.append(convert_hdf5_block)



# 

# Run whole pipeline

rnaseq_pipeline.output_commands(args.output_dir)

# rnaseq_pipeline.run_pipeline()



