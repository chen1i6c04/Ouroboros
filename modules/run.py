import os
import csv
from tempfile import NamedTemporaryFile
from loguru import logger
from .utils import syscall, fasta_scan, fastq_scan


def run_bwa_mem2(assembly, fastq, output, num_threads):
    syscall(f"bwa-mem2 index {assembly}")
    syscall(f"bwa-mem2 mem -t {num_threads} -a {assembly} {fastq} > {output}")
    syscall(f"rm {assembly}.*")


def run_polypolish(assembly, short_1, short_2, output_dir, num_threads):
    align_1 = os.path.join(output_dir, 'alignments_1.sam')
    align_2 = os.path.join(output_dir, 'alignments_2.sam')
    filtered_1 = os.path.join(output_dir, 'filtered_1.sam')
    filtered_2 = os.path.join(output_dir, 'filtered_2.sam')
    report = os.path.join(output_dir, 'polypolish.report')
    output = os.path.join(output_dir, '4_polypolish.fasta')
    run_bwa_mem2(assembly, short_1, align_1, num_threads)
    run_bwa_mem2(assembly, short_2, align_2, num_threads)

    handle = open(report, 'w')
    p = syscall(
        f'polypolish filter --in1 {align_1} --in2 {align_2} --out1 {filtered_1} --out2 {filtered_2}', stderr=True
    )
    handle.write(p.stderr)
    p = syscall(
        f"polypolish polish {assembly} {filtered_1} {filtered_2} | "
        f"sed 's/polypolish//g' | "
        f"seqkit sort -l -r -o {output}", stderr=True
    )
    handle.write(p.stderr)
    handle.close()

    syscall(f"sed -r 's/\x1b\[[0-9;]*m//g' -i {report}")
    for f in (align_1, align_2, filtered_1, filtered_2):
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
def run_pypolca(assembly, short_1, short_2, output_dir, depth, num_threads):
    """
    POLCA is a polishing tool in MaSuRCA (Maryland Super Read Cabog Assembler)
    https://github.com/alekseyzimin/masurca#polca
    """
    cmd = f"pypolca run -a {assembly} -1 {short_1} -2 {short_2} -t {num_threads} -o {output_dir}"
    if depth < 25:
        cmd += " --careful"
    logger.info(f"Running : {cmd}")
    syscall(cmd)


def run_flye(reads, outdir, threads, high_quality):
    input_type = '--nano-hq' if high_quality else '--nano-raw'
    cmd = f'flye {input_type} {reads} -o {outdir} -t {threads}'
    logger.info(f"Running : {cmd}")
    syscall(cmd)


def run_medaka(assembly, reads, outdir, model, threads):
    cmd = f"medaka_consensus -i {reads} -d {assembly} -o {outdir} -m {model} -t {threads} -f"
    logger.info(f"Running : {cmd}")
    syscall(cmd)
    os.remove(assembly + '.fai')
    os.remove(assembly + '.map-ont.mmi')
