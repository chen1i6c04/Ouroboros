import sys
import subprocess
from tempfile import TemporaryDirectory
from loguru import logger


def syscall(cmd):
    shell = True if isinstance(cmd, str) else False
    child_process = subprocess.run(
        cmd, shell=shell, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
    )
    if child_process.returncode:
        logger.error(f"Command {cmd} is fail")
        sys.exit('Abort')
    return child_process


def medaka_model_check(model):
    output = syscall(['medaka', 'tools', 'list_models']).stdout
    available_models = set(output.splitlines()[0].replace('Available: ', '').split(', '))
    return model in available_models


def estimate_genome_size(fastq_file, num_threads):
    with TemporaryDirectory() as tmpdir:
        output = syscall(
            f"kmc -sm -fm -t{num_threads} -k21 -ci10 {fastq_file} {tmpdir}/kmc {tmpdir} | "
            f"grep 'No. of unique counted k-mers' | "
            f"awk '{{print $NF}}'",
        ).stdout
    return int(output.strip())


def read_alignments(assembly, short_reads, alignments, num_threads):
    syscall(f"bwa-mem2 index {assembly}")
    syscall(f"bwa-mem2 mem -t {num_threads} -a {assembly} {short_reads} > {alignments}")
