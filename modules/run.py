import os
import csv
from tempfile import NamedTemporaryFile
from loguru import logger
from .utils import syscall


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
        f"polypolish_insert_filter.py --in1 {align_1} --in2 {align_2} --out1 {filtered_1} --out2 {filtered_2}",
        stderr=True
    )
    handle.write(p.stderr)
    p = syscall(
        f"polypolish {assembly} {filtered_1} {filtered_2} | sed 's/_polypolish//g' > {output}",
        stderr=True
    )
    handle.write(p.stderr)
    handle.close()

    syscall(f"sed -r 's/\x1b\[[0-9;]*m//g' -i {report}")
    for f in (align_1, align_2, filtered_1, filtered_2):
        os.remove(f)
    return output


@logger.catch
def run_dnaapler(assembly, outdir, asm_info, threads):
    with open(asm_info) as handle, NamedTemporaryFile('w') as tmpfile:
        spamreader = csv.reader(handle, delimiter='\t')
        next(spamreader)  # ignore header
        for row in spamreader:
            if row[3] == 'N':
                tmpfile.write(row[0] + '\n')
        tmpfile.flush()
        cmd = ['dnaapler', 'all', '-i', assembly, '-o', outdir, '-t', str(threads), '--ignore', tmpfile.name]
        syscall(cmd)
    return os.path.join(outdir, 'dnaapler_reoriented.fasta')


@logger.catch
def run_polca(assembly, short_1, short_2, output_dir, num_threads):
    """
    POLCA is a polishing tool in MaSuRCA (Maryland Super Read Cabog Assembler)
    https://github.com/alekseyzimin/masurca#polca
    """
    cmd = ['pypolca', 'run', '-a', assembly, '-1', short_1, '-2', short_2, '-t', str(num_threads), '-o', output_dir]
    syscall(cmd)


def run_flye(reads, outdir, threads, high_quality):
    input_type = '--nano-hq' if high_quality else '--nano-raw'
    cmd = ['flye', input_type, reads, '-o', outdir, '-t', str(threads)]
    logger.info(f"Flye command: {' '.join(cmd)}")
    syscall(cmd)
    log = os.path.join(outdir, 'flye.log')
    assembly = os.path.join(outdir, 'assembly.fasta')
    assembly_info = os.path.join(outdir, 'assembly_info.txt')
    assembly_graph = os.path.join(outdir, 'assembly_graph.gfa')
    return assembly, assembly_info, assembly_graph, log


def run_medaka(assembly, reads, outdir, model, threads):
    cmd = f"medaka_consensus -i {reads} -d {assembly} -o {outdir} -m {model} -t {threads}"
    logger.info("Medaka command: " + cmd)
    syscall(cmd)
    return os.path.join(outdir, 'consensus.fasta')


def run_plassembler(long, outdir, database, flye_asm, flye_info, threads, **kwargs):
    short_1, short_2 = kwargs.get("short_1", False), kwargs.get("short_2", False)
    mode = "run" if short_1 and short_2 else "long"
    cmd = [
            'plassembler', mode, '--skip_qc', '--keep_chromosome', '-t', str(threads), '-d', database, '-l', long,
            '-o', outdir, '--flye_assembly', flye_asm, '--flye_info', flye_info
        ]
    if short_1 and short_2:
        cmd += ['-1', short_1, '-2', short_2]
    syscall(cmd)
    chromosome = os.path.join(outdir, 'chromosome.fasta')
    plasmids = os.path.join(outdir, 'plassembler_plasmids.fasta')
    return chromosome, plasmids
