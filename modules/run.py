import os
import re
import csv
from tempfile import NamedTemporaryFile
from loguru import logger
from .utils import syscall

def paired_mapping(paired_1, paired_2, alignment_1, alignment_2, reference, num_threads):
    syscall(f"bwa index {reference}")
    syscall(f"bwa mem -t {num_threads} -a {reference} {paired_1} > {alignment_1}")
    syscall(f"bwa mem -t {num_threads} -a {reference} {paired_2} > {alignment_2}")
    syscall(f"rm {reference}.*")


def run_polypolish(assembly, short_1, short_2, output_dir, num_threads):
    alignments_1 = os.path.join(output_dir, 'alignments_1.sam')
    alignments_2 = os.path.join(output_dir, 'alignments_2.sam')
    filtered_1 = os.path.join(output_dir, 'filtered_1.sam')
    filtered_2 = os.path.join(output_dir, 'filtered_2.sam')
    report = os.path.join(output_dir, 'polypolish.report')
    output = os.path.join(output_dir, '4_polypolish.fasta')
    paired_mapping(short_1, short_2, alignments_1, alignments_2, assembly, num_threads)

    with open(report, 'w') as handle, NamedTemporaryFile('w') as tmpfile:
        child_process = syscall(
            f'polypolish filter --in1 {alignments_1} --in2 {alignments_2} --out1 {filtered_1} --out2 {filtered_2}',
            stderr=True
        )
        handle.write(re.sub(r'\x1b\[[0-9;]*m', '', child_process.stderr))
        child_process = syscall(
            f"polypolish polish {assembly} {filtered_1} {filtered_2}", stdout=True, stderr=True
        )
        handle.write(re.sub(r'\x1b\[[0-9;]*m', '', child_process.stderr))
        tmpfile.write(re.sub('polypolish', '', child_process.stdout))
        tmpfile.flush()
        syscall(f'seqkit sort -r -l {tmpfile.name} -o {output}')
    for f in (alignments_1, alignments_2, filtered_1, filtered_2):
        os.remove(f)
    return output


@logger.catch
def run_dnaapler(flye_output, outdir, threads):
    assembly = os.path.join(flye_output, 'assembly.fasta')
    assembly_info = os.path.join(flye_output, 'assembly_info.txt')
    with open(assembly_info) as handle, NamedTemporaryFile('w') as tmpfile:
        spamreader = csv.reader(handle, delimiter='\t')
        next(spamreader)  # ignore header
        for row in spamreader:
            if row[3] == 'N':
                tmpfile.write(row[0] + '\n')
        tmpfile.flush()
        cmd = f"dnaapler all -i {assembly} -o {outdir} -t {threads} --ignore {tmpfile.name}"
        logger.info(f"Running : {cmd}")
        syscall(cmd)


@logger.catch
def run_pypolca(assembly, short_1, short_2, output_dir, num_threads):
    """
    POLCA is a polishing tool in MaSuRCA (Maryland Super Read Cabog Assembler)
    https://github.com/alekseyzimin/masurca#polca
    """
    cmd = f"pypolca run -a {assembly} -1 {short_1} -2 {short_2} -t {num_threads} -o {output_dir} --careful"
    logger.info(f"Running : {cmd}")
    syscall(cmd)


def run_flye(reads, outdir, threads):
    cmd = f'flye --nano-hq {reads} -o {outdir} -t {threads}'
    logger.info(f"Running : {cmd}")
    syscall(cmd)


def run_medaka(assembly, reads, outdir, model, threads):
    cmd = f"medaka_consensus -i {reads} -d {assembly} -o {outdir} -m {model} -t {threads} -f"
    logger.info(f"Running : {cmd}")
    syscall(cmd)
    os.remove(assembly + '.fai')
    os.remove(assembly + '.map-ont.mmi')
