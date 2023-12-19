import os
import sys
import gzip
import json
import subprocess
from tempfile import TemporaryDirectory
from Bio import SeqIO
from loguru import logger


def syscall(cmd, stdout=False, stderr=True):
    if stdout:
        stdout_str = subprocess.PIPE
    else:
        stdout_str = None
    if stderr:
        stderr_str = subprocess.PIPE
    else:
        stderr_str = None
    shell, executable = (True, "/bin/bash") if isinstance(cmd, str) else (False, None)
    child_process = subprocess.run(
        cmd, shell=shell, stdout=stdout_str, stderr=stderr_str, universal_newlines=True, executable=executable,
    )
    if child_process.returncode:
        logger.error(f"Command {cmd} is fail")
        logger.error(child_process.stderr)
        sys.exit('Abort')
    return child_process


def validate_medaka_model(model):
    output = syscall(['medaka', 'tools', 'list_models'], stdout=True).stdout
    available_models = set(output.splitlines()[0].replace('Available: ', '').split(', '))
    if not (model in available_models):
        logger.error(f"Medaka model {model} unavailable")
        sys.exit('Abort')


def validate_fastq(file):
    """Checks the input fastq is really a fastq
        :param file: fastq file
    :return: zipped - Boolean whether the input fastq is gzipped.
    """

    # to get extension
    filename, file_extension = os.path.splitext(file)
    if file_extension == ".gz":
        # if gzipped
        with gzip.open(file, "rt") as handle:
            fastq = SeqIO.parse(handle, "fastq")
            if any(fastq):
                logger.info(f"FASTQ {file} checked")
            else:
                logger.error(f"Input file {file} is not in the FASTQ format.")
                sys.exit('Abort')
    else:
        with open(file, "r") as handle:
            fastq = SeqIO.parse(handle, "fastq")
            if any(fastq):
                logger.info(f"FASTQ {file} checked")
            else:
                logger.error(f"Input file {file} is not in the FASTQ format.")
                sys.exit('Abort')


def fastq_scan(file):
    p = syscall(f"nanoq -s -f -j -i {file}", stdout=True)
    return json.loads(p.stdout)


def estimate_genome_size(fastq_file, num_threads):
    with TemporaryDirectory() as tmpdir:
        output = syscall(
            f"kmc -sm -fq -t{num_threads} -k21 -ci10 {fastq_file} {tmpdir}/kmc {tmpdir} | "
            f"grep 'No. of unique counted k-mers' | "
            f"awk '{{print $NF}}'", stdout=True
        ).stdout
    return int(output.strip())


def exclude_target_from_single_end(input_reads, output_reads, target, threads):
    cmd = (f"minimap2 --secondary=no -L -t {threads} -ax map-ont {target} {input_reads} | "
           f"samtools sort -@ {threads} -O BAM - | "
           f"samtools view -@ {threads} -f 4 -O BAM - | "
           f"samtools fastq -@ {threads} - | "
           f"pigz -9 -p {threads} > {output_reads}")
    syscall(cmd)


def exclude_target_from_paired_end(paired_1, paired_2, output_1, output_2, target, threads):
    cmd = f"minimap2 -t {threads} -ax sr {target} {paired_1} {paired_2} | " \
          f"samtools sort -@ {threads} -n -O BAM - | " \
          f"samtools view -@ {threads} -f 12 -O BAM - | " \
          f"samtools fastq -@ {threads} -1 {output_1} -2 {output_2} -0 /dev/null -s /dev/null -n -"
    syscall(cmd)
