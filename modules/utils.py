import re
import sys
import gzip
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
    def is_gzip(f):
        with open(f, 'rb') as f:
            signature = f.read(2)
        return signature == b'\x1f\x8b'
    # to get extension
    open_function = gzip.open if is_gzip(file) else open
    with open_function(file, 'rt') as handle:
        if any(SeqIO.parse(handle, 'fastq')):
            logger.info(f"FASTQ {file} checked")
        else:
            logger.error(f"Input file {file} is not in the FASTQ format.")
            sys.exit('Abort')


def parse_genome_size(pattern):
    unit_map = {'K': 1e3, 'M': 1e6, 'G': 1e9}
    result = re.fullmatch(r'^([\d.]+)([KMG])', pattern)
    if result is None:
        return
    else:
        value, unit = result.groups()
        return int(float(value) * unit_map[unit])


def estimate_genome_size(fastq_file, num_threads):
    with TemporaryDirectory() as tmpdir:
        output = syscall(
            f"kmc -sm -fq -t{num_threads} -k21 -ci10 {fastq_file} {tmpdir}/kmc {tmpdir} | "
            f"grep 'No. of unique counted k-mers' | "
            f"awk '{{print $NF}}'", stdout=True
        ).stdout
    genome_size = int(output.strip())
    logger.info(f"Estimated genome size was {genome_size}bp.")
    return genome_size


def exclude_target_from_single_end(input_reads, output_reads, target, threads):
    cmd = (f"minimap2 --secondary=no -L -t {threads} -ax map-ont {target} {input_reads} | "
           f"samtools sort -@ {threads} -O BAM - | "
           f"samtools view -@ {threads} -f 4 -O BAM - | "
           f"samtools fastq -@ {threads} - | "
           f"pigz -9 -p {threads} > {output_reads}")
    syscall(cmd)


def exclude_target_from_paired_end(paired_1, paired_2, output_1, output_2, target, threads):
    cmd = (
        f"minimap2 -t {threads} -ax sr {target} {paired_1} {paired_2} | "
        f"samtools sort -@ {threads} -n -O BAM - | "
        f"samtools view -@ {threads} -f 12 -O BAM - | "
        f"samtools fastq -@ {threads} -1 {output_1} -2 {output_2} -0 /dev/null -s /dev/null -n -"
    )
    syscall(cmd)
