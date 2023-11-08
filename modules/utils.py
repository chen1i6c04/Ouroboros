import sys
import subprocess
from tempfile import TemporaryDirectory
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
    shell = True if isinstance(cmd, str) else False
    executable = "/bin/bash" if shell else None
    child_process = subprocess.run(
        cmd, shell=shell, stdout=stdout_str, stderr=stderr_str, universal_newlines=True, executable=executable,
    )
    if child_process.returncode:
        logger.error(f"Command {cmd} is fail")
        logger.error(child_process.stderr)
        sys.exit('Abort')
    return child_process


def medaka_model_check(model):
    output = syscall(['medaka', 'tools', 'list_models'], stdout=True).stdout
    available_models = set(output.splitlines()[0].replace('Available: ', '').split(', '))
    return model in available_models


def estimate_genome_size(fastq_file, num_threads):
    with TemporaryDirectory() as tmpdir:
        output = syscall(
            f"kmc -sm -fq -t{num_threads} -k21 -ci10 {fastq_file} {tmpdir}/kmc {tmpdir} | "
            f"grep 'No. of unique counted k-mers' | "
            f"awk '{{print $NF}}'", stdout=True
        ).stdout
    return int(output.strip())


def read_alignments(assembly, short_reads, alignments, num_threads):
    syscall(f"bwa-mem2 index {assembly}")
    syscall(f"bwa-mem2 mem -t {num_threads} -a {assembly} {short_reads} > {alignments}")


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
