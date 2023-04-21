import sys
import subprocess
from tempfile import TemporaryDirectory


def syscall(cmd):
    shell = True if isinstance(cmd, str) else False
    child_process = subprocess.run(
        cmd, shell=shell, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
    )
    if child_process.returncode:
        sys.exit(f"Command {cmd} is fail")
    return child_process


def medaka_model_check(model):
    output = syscall(['medaka', 'tools', 'list_models']).stdout
    model_list = set(output.splitlines()[0].replace('Available: ', '').split(', '))
    return model in model_list


def estimate_genome_size(fastq_file, num_threads):
    with TemporaryDirectory() as tmpdir:
        output = syscall(
            f"kmc -sm -fm -t{num_threads} -k21 -ci10 {fastq_file} {tmpdir}/kmc {tmpdir} | "
            f"grep 'No. of unique counted k-mers' | "
            f"awk '{{print $NF}}'",
        ).stdout
    return int(output.strip().split()[-1])

