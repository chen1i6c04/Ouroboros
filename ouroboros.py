import os
import re
import sys
import shutil
import argparse
from loguru import logger
from modules.utils import (syscall,
                           validate_medaka_model,
                           estimate_genome_size,
                           exclude_target_from_single_end,
                           validate_fastq,
                           fastq_scan)
from modules.run import run_polypolish, run_dnaapler, run_polca, run_flye, run_medaka, run_plassembler

LOCATION = os.path.dirname(os.path.abspath(__file__))


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
        'plassembler': 'plassembler -V',
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
            command.append("--" + argument)
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
    command = ' '.join(sys.argv)
    logger.info("You ran: " + command)
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
    logger.info("Long reads trimming.")
    cmd = f"porechop -i {input_file} -o {output_file} --format fastq.gz -t {threads}"
    syscall(cmd)


def long_reads_filter(input_file, output_dir, threads, min_length=1, min_quality=1, keep_percent=90):
    logger.info(f"Filter out reads length less than {min_length:,}bp and average quality score less {min_quality}")
    first_filter = os.path.join(output_dir, 'READS_1st_filt.fastq.gz')
    cmd = f'nanoq -i {input_file} -l {min_length} -q {min_quality} | pigz -6 -p {threads} > {first_filter}'
    syscall(cmd)
    logger.info(f"Keep only {keep_percent} percentage of the best reads")
    second_filter = os.path.join(output_dir, 'READS_2nd_filt.fastq.gz')
    cmd = f'filtlong --keep_percent {keep_percent} {first_filter} | pigz -6 -p {threads} > {second_filter}'
    syscall(cmd)
    return second_filter


def short_reads_trimming(input_1, input_2, output_1, output_2, threads):
    logger.info("Trimming short reads.")
    cmd = [
        'fastp', '-i', input_1, '-I', input_2, '-o', output_1, '-O', output_2,
        '--thread', str(threads), '--detect_adapter_for_pe', '-j', '/dev/null', '-h', '/dev/null', '-t', '1',
    ]
    syscall(cmd)


def short_reads_polish(assembly, short_1, short_2, outdir, threads):
    polca_dir = os.path.join(outdir, 'polca')
    polypolish_output = run_polypolish(assembly, short_1, short_2, outdir, threads)
    run_polca(polypolish_output, short_1, short_2, polca_dir, threads)
    shutil.copyfile(
        os.path.join(polca_dir, 'polca_corrected.fasta'),
        os.path.join(outdir, '5_polca.fasta')
    )
    shutil.copy(
        os.path.join(polca_dir, 'polca.report'),
        outdir
    )
    shutil.rmtree(polca_dir)


def subsampling(infile, outfile, gsize, depth):
    bases = gsize * depth
    cmd = ['rasusa', '-b', str(bases), '-i', infile, '-o', outfile]
    syscall(cmd)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile',
                        required=True,
                        help='Input Nanopore FASTQ')
    parser.add_argument('-1', '--short_1', default='',
                        help='Read 1 FASTQ to use for polishing')
    parser.add_argument('-2', '--short_2', default='',
                        help='Read 2 FASTQ to use for polishing')
    parser.add_argument('-o', '--outdir',
                        required=True,
                        help='Specify directory in which output has to be created.')
    parser.add_argument('-t', '--num-threads', default=1, type=int,
                        help='Set the allowed number of threads to be used by the script (default: 1)')
    parser.add_argument('-l', '--min-length', metavar='', default=1, type=int)
    parser.add_argument('-q', '--min-quality', metavar='', default=1, type=int)
    parser.add_argument('--nofilter', action='store_true', default=False,
                        help='Disable long reads length filtering')
    parser.add_argument('--notrim', action='store_true', default=False,
                        help='Disable long reads trimming.')
    parser.add_argument('-x', '--depth',
                        default=50, type=int,
                        help='Sub-sample reads to this depth. Disable with --depth 0 (default: 50)')
    parser.add_argument('-g', '--gsize', default=None, metavar='',
                        help='Estimated genome size eg. 3.2M <blank=AUTO> (default: "")')
    parser.add_argument('--meta', action='store_true',
                        help='Metagenome / uneven coverage')
    parser.add_argument('--medaka_model', default='r1041_e82_400bps_sup_v4.3.0',
                        help='The model to be used by Medaka (default: r1041_e82_400bps_sup_v4.3.0)')
    parser.add_argument('--medaka_opt',
                        help='Additional options to be given to Medaka')
    parser.add_argument('--hq', action='store_true',
                        help="Flye will use '--nano-hq' instead of --nano-raw")
    parser.add_argument('--contaminants',
                        help='Contaminants FASTA file or Minimap2 index to map long reads against to filter out.')
    parser.add_argument('--plasmids',
                        help='Run Plassembler assemble plasmids.')
    args = parser.parse_args()

    if (args.short_1 and not args.short_2) or (not args.short_1 and args.short_2):
        parser.error("-1 and -2 must be used together")

    validate_medaka_model(args.medaka_model)

    initialize(args)

    validate_fastq(args.infile)
    if args.short_1 and args.short_2:
        validate_fastq(args.short_1)
        validate_fastq(args.short_2)

    reads = args.infile
    if not args.notrim:
        trimmed_reads = os.path.join(args.outdir, 'READS_trim.fastq.gz')
        long_reads_trimming(reads, trimmed_reads, args.num_threads)
        reads = trimmed_reads
    if not args.nofilter:
        filtered_reads = long_reads_filter(reads, args.outdir, args.num_threads, args.min_length, args.min_quality)
        reads = filtered_reads
    if args.contaminants:
        logger.info("Clean contaminants.")
        cleaned_reads = os.path.join(args.outdir, 'READS_clean.fastq.gz')
        exclude_target_from_single_end(
            input_reads=reads,
            output_reads=cleaned_reads,
            target=args.contaminants,
            threads=args.num_threads
        )
        reads = cleaned_reads
    fastq_stats = fastq_scan(reads)

    if not args.meta:
        total_bases = fastq_stats['bases']
        if args.gsize:
            gsize = parse_genome_size(args.gsize)
            logger.info(f"Using genome size was {gsize}bp.")
        else:
            gsize = estimate_genome_size(reads, args.num_threads)
            logger.info(f"Estimated genome size was {gsize}bp.")
        origin_depth = total_bases / gsize
        logger.info(f'Estimated sequencing depth: {origin_depth:.0f} x')
        if args.depth:
            if origin_depth > args.depth:
                logger.info(f"Subsampling reads from {origin_depth:.0f}x to {args.depth}x.")
                sub_reads = os.path.join(args.outdir, 'READS.sub.fastq')
                subsampling(reads, sub_reads, gsize, args.depth)
            else:
                logger.info("No read depth reduction requested or necessary.")
                sub_reads = reads
        else:
            sub_reads = reads
    else:
        logger.info("Flag '--meta' was be set, won't subsampling reads.")
        sub_reads = reads

    logger.info("Assembling reads with Flye")
    flye_dir = os.path.join(args.outdir, 'flye')
    flye_asm, flye_info, flye_graph, flye_log = run_flye(sub_reads, flye_dir, args.num_threads, args.meta, args.hq)

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
    reoriented_asm = run_dnaapler(flye_asm, dnaapler_dir, flye_info, args.num_threads)

    logger.info("Polishing with medaka.")
    medaka_dir = os.path.join(args.outdir, 'medaka')
    medaka_asm = run_medaka(
        assembly=reoriented_asm,
        reads=reads,
        outdir=medaka_dir,
        model=args.medaka_model,
        options=args.medaka_opt,
        threads=args.num_threads
    )
    shutil.copyfile(medaka_asm, os.path.join(args.outdir, '2_medaka.fasta'))

    if args.short_1 and args.short_2:
        trim_1, trim_2 = os.path.join(args.outdir, 'READS_1.trim.fastq'), os.path.join(args.outdir, 'READS_2.trim.fastq')
        short_reads_trimming(args.short_1, args.short_2, trim_1, trim_2, args.num_threads)

        logger.info("Plasmid assembly.")
        plassembler_dir = os.path.join(args.outdir, 'plassembler')
        chromosome, plasmids = run_plassembler(
            long=sub_reads, short_1=trim_1, short_2=trim_2, outdir=plassembler_dir, flye_asm=medaka_asm,
            flye_info=flye_info, threads=args.num_threads, database=os.path.join(LOCATION, 'plassembler_db')
        )
        plassembler_asm = os.path.join(args.outdir, '3_plassembler.fasta')
        syscall(f"cat {chromosome} {plasmids} > {plassembler_asm}")
        logger.info("Polishing with short reads")
        short_reads_polish(plassembler_asm, trim_1, trim_2, args.outdir, args.num_threads)
    else:
        logger.info("Plasmid assembly.")
        plassembler_dir = os.path.join(args.outdir, 'plassembler')
        chromosome, plasmids = run_plassembler(
            long=sub_reads, outdir=plassembler_dir, flye_asm=medaka_asm, flye_info=flye_info, threads=args.num_threads,
            database=os.path.join(LOCATION, 'plassembler_db')
        )
        plassembler_asm = os.path.join(args.outdir, '3_plassembler.fasta')
        syscall(f"cat {chromosome} {plasmids} > {plassembler_asm}")

    for d in (flye_dir, medaka_dir, dnaapler_dir):
        shutil.rmtree(d)
    syscall(f"rm {os.path.join(args.outdir, 'READS*')}")
    logger.info("Done")


if __name__ == '__main__':
    main()
