import os
import re
import csv
import sys
import gzip
import subprocess
from tempfile import TemporaryDirectory
from Bio import SeqIO
from loguru import logger


def syscall(cmd, stdout=False, stderr=False):
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
        cmd = (f"kmc -sm -t{num_threads} -k21 -ci10 {fastq_file} {tmpdir}/kmc {tmpdir} | "
               f"grep -oP 'No. of unique counted k-mers\s+:\s+\K\d+'")
        child_process = syscall(cmd, stdout=True)
    return int(child_process.stdout)


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


def annotate_assembly(assembly, assembly_info):
    coverage_dict = dict()
    circular_dict = dict()
    with open(assembly_info) as handle:
        spamreader = csv.reader(handle, delimiter='\t')
        next(spamreader)  # ignore header
        for row in spamreader:
            coverage_dict[row[0]] = row[2]
            circular_dict[row[0]] = row[3]
    records = []
    for record in SeqIO.parse(assembly, 'fasta'):
        length = len(record)
        coverage = coverage_dict[record.id]
        circular = circular_dict[record.id]
        record.description = f"length={length} coverage={coverage}{' circular=true' if circular == 'Y' else ''}"
        records.append(record)
    records = sorted(records, key=lambda record: len(record.seq), reverse=True)
    return records


def merge_chromosome_amd_plasmid(dirpath):
    assembly = os.path.join(dirpath, '2_medaka.fasta')
    assembly_info = os.path.join(dirpath, 'flye_info.txt')
    records = annotate_assembly(assembly, assembly_info)
    records = list(filter(lambda rec: len(rec.seq)>=1e6, records))
    if len(records) == 1:
        records[0].id = records[0].name = 'chromosome'
    else:
        for idx, record in enumerate(records, 1):
            record.id = record.name = f'chromosome_{idx}'
    plasmids = os.path.join(dirpath, 'plassembler', 'plassembler_plasmids.fasta')
    with open(os.path.join(dirpath, '3_plassembler.fasta'), 'w') as handle:
        SeqIO.write(records, handle, 'fasta')
        for record in SeqIO.parse(plasmids, 'fasta'):
            SeqIO.write(record, handle, 'fasta')


def collate(query, subject, result):
    describe = {}
    for record in SeqIO.parse(subject, 'fasta'):
        describe[record.id] = record.description
    records = []
    for record in SeqIO.parse(query, 'fasta'):
        description = describe[record.id]
        description = re.sub('length=\d+', f'length={len(record)}', description)
        record.description = description
        records.append(record)
    SeqIO.write(records, result, 'fasta')
