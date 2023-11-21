import os
import re
import sys
import json
import shutil
import argparse
from tempfile import NamedTemporaryFile

import pandas as pd
from loguru import logger
from modules.utils import (syscall,
                           medaka_model_check,
                           estimate_genome_size,
                           read_alignments,
                           exclude_target_from_single_end)


def check_dependency():
    version = {
        'rasusa': 'rasusa --version',
        'nanoq': 'nanoq --version',
        'filtlong': 'filtlong --version',
        'flye': 'flye -v',
        'medaka': 'medaka --version',
        "KMC": "kmc | grep K-Mer",
        'polypolish': 'polypolish --version',
        'bwa-mem2': 'bwa-mem2 version 2> /dev/null',
        'pypolca': 'pypolca -V',
        'dnaapler': 'dnaapler -V',
    }
    for program_name, cmd in version.items():
        p = syscall(cmd, stdout=True)
        if p.returncode:
            logger.error(f"Could not determine version of {program_name}")
            sys.exit("Abort")
        else:
            version = p.stdout.strip()
            logger.info(f"Using {program_name:11} | {version}")


def print_used_values(args):
    command = []
    args = vars(args)
    for argument, value in args.items():
        if value:
            command.append("--" + argument.replace('_', '-'))
            if isinstance(value, bool) is False:
                command.append(str(value))
    command.insert(0, sys.argv[0])
    logger.info("You ran: " + ' '.join(command))


def initialize(args):
    os.makedirs(args.outdir, exist_ok=True)
    logfile = os.path.join(args.outdir, 'ouroboros.log')
    fmt = "{time:YYYY-MM-DD HH:mm:ss} [{level}] {message}"
    logger.add(logfile, format=fmt, level='INFO')
    logger.add(sys.stderr, format=fmt, level='ERROR')
    print_used_values(args)
    check_dependency()


def parse_genome_size(pattern):
    unit_map = {'K': 1e3, 'M': 1e6, 'G': 1e9}
    result = re.fullmatch(r'^([\d.]+)([KMG])', pattern)
    if result is None:
        logger.error(f"Couldn't parse {pattern}")
        sys.exit("Aborted")
    else:
        value, unit = result.groups()
        return int(float(value) * unit_map[unit])


def long_reads_trimming(input_file, output_file, threads):
    cmd = f"porechop -i {input_file} -o {output_file} --format fastq.gz -t {threads}"
    syscall(cmd)


def long_reads_filter(input_file, output_dir, min_length=1, min_quality=1, keep_percent=90):
    logger.info(f"Filter out reads length less than {min_length} and average quality score less {min_quality}")
    first_filter = os.path.join(output_dir, 'READS_1st_filt.fastq.gz')
    cmd = f'nanoq -i {input_file} -l {min_length} -q {min_quality} -O g -o {first_filter}'
    syscall(cmd)
    logger.info(f"Keep only {keep_percent} percentage of the best reads")
    second_filter = os.path.join(output_dir, 'READS_2nd_filt.fastq.gz')
    cmd = f'filtlong --keep_percent {keep_percent} {first_filter} | pigz  > {second_filter}'
    syscall(cmd)
    return second_filter


def short_reads_trimming(input_1, input_2, output_1, output_2, threads):
    cmd = [
        'fastp', '-i', input_1, '-I', input_2, '-o', output_1, '-O', output_2,
        '--thread', str(threads), '--detect_adapter_for_pe', '-j', '/dev/null', '-h', '/dev/null', '-t', '1',
    ]
    syscall(cmd)


def subsampling(infile, outfile, gsize, depth):
    bases = gsize * depth
    cmd = ['rasusa', '-b', str(bases), '-i', infile, '-o', outfile]
    syscall(cmd)


def run_flye(reads, outdir, threads, meta, high_quality):
    input_type = '--nano-hq' if high_quality else '--nano-raw'
    cmd = ['flye', input_type, reads, '-o', outdir, '-t', str(threads)]
    if meta:
        cmd += ['--meta']
    logger.info(f"Flye command: {' '.join(cmd)}")
    syscall(cmd)
    log = os.path.join(outdir, 'flye.log')
    assembly = os.path.join(outdir, 'assembly.fasta')
    assembly_info = os.path.join(outdir, 'assembly_info.txt')
    assembly_graph = os.path.join(outdir, 'assembly_graph.gfa')
    return assembly, assembly_info, assembly_graph, log


@logger.catch
def run_polca(assembly, short_1, short_2, output_dir, num_threads):
    """
    POLCA is a polishing tool in MaSuRCA (Maryland Super Read Cabog Assembler)
    https://github.com/alekseyzimin/masurca#polca
    """
    cmd = ['pypolca', 'run', '-a', assembly, '-1', short_1, '-2', short_2, '-t', str(num_threads), '-o', output_dir]
    syscall(cmd)


def reorient(assembly, outdir, asm_info, threads):
    df = pd.read_csv(asm_info, sep='\t')
    s = df[df['circ.'] == 'N']['#seq_name']
    with NamedTemporaryFile('w') as tmpfile:
        tmpfile.write(s.to_csv(index=False, header=False))
        tmpfile.flush()
        cmd = ['dnaapler', 'all', '-i', assembly, '-o', outdir, '-t', str(threads), '--ignore', tmpfile.name]
        syscall(cmd)
    return os.path.join(outdir, 'dnaapler_reoriented.fasta')


def run_medaka(assembly, reads, outdir, model, options, threads):
    cmd = ['medaka_consensus', '-i', reads, '-d', assembly, '-o', outdir, '-m', model,
           '-t', str(threads)]
    if options:
        cmd += options
    syscall(cmd)
    return os.path.join(outdir, 'consensus.fasta')


def run_polypolish(assembly, short_reads_1, short_reads_2, output_dir, num_threads):
    alignments_1 = os.path.join(output_dir, 'alignments_1.sam')
    alignments_2 = os.path.join(output_dir, 'alignments_2.sam')
    filtered_1 = os.path.join(output_dir, 'filtered_1.sam')
    filtered_2 = os.path.join(output_dir, 'filtered_2.sam')
    read_alignments(assembly, short_reads_1, alignments_1, num_threads)
    read_alignments(assembly, short_reads_2, alignments_2, num_threads)
    polypolish_report = os.path.join(output_dir, 'polypolish.report')
    handle = open(polypolish_report, 'w')
    p = syscall(
        f"polypolish_insert_filter.py --in1 {alignments_1} --in2 {alignments_2} --out1 {filtered_1} --out2 {filtered_2}",
        stderr=True
    )
    handle.write(p.stderr)
    polypolish_output = os.path.join(output_dir, '3_polypolish.fasta')
    p = syscall(
        f"polypolish {assembly} {filtered_1} {filtered_2} | sed 's/_polypolish//g' > {polypolish_output}",
        stderr=True
    )
    handle.write(p.stderr)
    handle.close()
    syscall(f"sed -r 's/\x1b\[[0-9;]*m//g' -i {polypolish_report}")
    for f in (alignments_1, alignments_2, filtered_1, filtered_2):
        os.remove(f)
    return polypolish_output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile',
                        required=True,
                        help='Input Nanopore FASTQ')
    parser.add_argument('-o', '--outdir',
                        required=True,
                        help='Specify directory in which output has to be created.')
    parser.add_argument('-t', '--num-threads', default=1, type=int,
                        help='Set the allowed number of threads to be used by the script (default: 1)')
    parser.add_argument('-l', '--min-length', metavar='', default=1, type=int)
    parser.add_argument('-q', '--min-quality', metavar='', default=1, type=int)
    parser.add_argument('--nofilter', action='store_true', default=False,
                        help='Disable read length filtering')
    parser.add_argument('--notrim', action='store_true', default=False,
                        help='Disable read trimming.')
    parser.add_argument('-1', '--short_1',
                        help='Read 1 FASTQ to use for polishing')
    parser.add_argument('-2', '--short_2',
                        help='Read 2 FASTQ to use for polishing')
    parser.add_argument('-x', '--depth',
                        default=50, type=int,
                        help='Sub-sample reads to this depth. Disable with --depth 0 (default: 50)')
    parser.add_argument('-g', '--gsize', default=None, metavar='',
                        help='Estimated genome size eg. 3.2M <blank=AUTO> (default: "")')
    parser.add_argument('--meta', action='store_true',
                        help='Metagenome / uneven coverage')
    parser.add_argument('--medaka_model', default='r1041_e82_400bps_sup_g615',
                        help='The model to be used by Medaka (default: r1041_e82_400bps_sup_g615)')
    parser.add_argument('--medaka_opt', nargs='+',
                        help='Additional options to be given to Medaka')
    parser.add_argument('--hq', action='store_true',
                        help="Flye will use '--nano-hq' instead of --nano-raw")
    parser.add_argument('--contaminants',
                        help='Contaminants FASTA file or Minimap2 index to map long reads against to filter out.')
    args = parser.parse_args()

    short_reads_polishing = False
    if (args.short_1 and not args.short_2) or (not args.short_1 and args.short_2):
        parser.error("--sr1 and --sr2 must be used together")
    elif args.short_1 and args.short_2:
        short_reads_polishing = True
    else:
        pass

    initialize(args)

    if medaka_model_check(args.medaka_model) is False:
        logger.error(f"Medaka model {args.medaka_model} unavailable")
        sys.exit("Aborted")

    reads = ''
    if not os.access(args.infile, os.R_OK):
        logger.error(f"{args.infile} is not a readable file")
        sys.exit("Aborted")
    if args.infile.endswith('.fastq.gz'):
        reads = os.path.join(args.outdir, 'READS.fastq.gz')
    elif args.infile.endswith('.fastq'):
        reads = os.path.join(args.outdir, 'READS.fastq')
    else:
        sys.exit()

    os.symlink(args.infile, reads)
    interm_reads = reads
    if not args.notrim:
        trimmed_reads = os.path.join(args.outdir, 'READS_trim.fastq.gz')
        long_reads_trimming(interm_reads, trimmed_reads, args.num_threads)
        interm_reads = trimmed_reads
    if not args.nofilter:
        filtered_reads = long_reads_filter(interm_reads, args.outdir, args.min_length, args.min_quality)
        interm_reads = filtered_reads
    if args.contaminants:
        logger.info("Clean contaminants.")
        cleaned_reads = os.path.join(args.outdir, 'READS_clean.fastq.gz')
        exclude_target_from_single_end(
            input_reads=interm_reads,
            output_reads=cleaned_reads,
            target=args.contaminants,
            threads=args.num_threads
        )
        interm_reads = cleaned_reads
    cmd = f'nanoq -s -j -i {interm_reads}'
    p = syscall(cmd, stdout=True)
    fastq_stats = json.loads(p.stdout)

    if not args.meta:
        total_bases = fastq_stats['bases']
        if args.gsize:
            gsize = parse_genome_size(args.gsize)
            logger.info(f"Using genome size was {gsize}bp.")
        else:
            gsize = estimate_genome_size(interm_reads, args.num_threads)
            logger.info(f"Estimated genome size was {gsize}bp.")
        origin_depth = total_bases / gsize
        logger.info(f'Estimated sequencing depth: {origin_depth:.0f} x')
        if args.depth:
            if origin_depth > args.depth:
                logger.info(f"Subsampling reads from {origin_depth:.0f}x to {args.depth}x.")
                sub_reads = os.path.join(args.outdir, 'READS.sub.fastq')
                subsampling(interm_reads, sub_reads, gsize, args.depth)
                final_reads = sub_reads
            else:
                logger.info("No read depth reduction requested or necessary.")
                final_reads = interm_reads
        else:
            final_reads = interm_reads
    else:
        logger.info("Flag '--meta' was be set, won't subsampling reads.")
        final_reads = interm_reads

    logger.info("Assembling reads with Flye")
    flye_dir = os.path.join(args.outdir, 'flye')
    flye_asm, flye_info, flye_graph, flye_log = run_flye(final_reads, flye_dir, args.num_threads, args.meta, args.hq)

    shutil.copyfile(
        flye_asm,
        os.path.join(args.outdir, '1_flye.fasta')
    )
    shutil.copyfile(
        flye_info,
        os.path.join(args.outdir, 'flye_info.txt')
    )
    shutil.copyfile(
        flye_graph,
        os.path.join(args.outdir, 'flye-unpolished.gfa')
    )
    shutil.copyfile(
        flye_log,
        os.path.join(args.outdir, 'flye.log')
    )
    logger.info('Reorients complete circular sequence.')
    dnaapler_dir = os.path.join(args.outdir, 'dnaapler')
    reoriented_asm = reorient(flye_asm, dnaapler_dir, flye_info, args.num_threads)

    logger.info("Polishing with medaka.")
    medaka_dir = os.path.join(args.outdir, 'medaka')
    medaka_asm = run_medaka(reoriented_asm, interm_reads, medaka_dir, args.medaka_model, args.medaka_opt, args.num_threads)
    shutil.copyfile(medaka_asm, os.path.join(args.outdir, '2_medaka.fasta'))

    if short_reads_polishing:
        logger.info("Trimming short reads")
        trim_1, trim_2 = os.path.join(args.outdir, 'R1.trim.fastq'), os.path.join(args.outdir, 'R2.trim.fastq')
        short_reads_trimming(args.short_1, args.short_2, trim_1, trim_2, args.num_threads)
        logger.info("Polishing with short reads")
        polypolish_output = run_polypolish(medaka_asm, trim_1, trim_2, args.outdir, args.num_threads)
        polca_dir = os.path.join(args.outdir, 'polca')
        run_polca(polypolish_output, trim_1, trim_2, polca_dir, args.num_threads)
        shutil.copyfile(
            os.path.join(polca_dir, 'polca_corrected.fasta'),
            os.path.join(args.outdir, '4_polca.fasta')
        )
        shutil.copy(
            os.path.join(polca_dir, 'polca.report'),
            args.outdir
        )
        shutil.rmtree(polca_dir)
    for d in (flye_dir, medaka_dir):
        shutil.rmtree(d)
    syscall(f"rm {os.path.join(args.outdir, 'READS*')}")
    logger.info("Done")


if __name__ == '__main__':
    main()
