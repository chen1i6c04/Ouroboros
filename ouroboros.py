import os
import sys
import shutil
import argparse
from loguru import logger
from modules.utils import (syscall,
                           validate_medaka_model,
                           estimate_genome_size,
                           exclude_target_from_single_end,
                           validate_fastq,
                           fastq_scan,
                           parse_genome_size)
from modules.filter import length_filter, quality_filter
from modules.trim import trim_short_reads, trim_long_reads
from modules.run import run_polypolish, run_dnaapler, run_pypolca, run_flye, run_medaka
from modules.plassembler import run_plassembler

__location__ = os.path.dirname(os.path.abspath(__file__))
__version__ = 'v1.0.0'


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


def begin(outdir):
    os.makedirs(outdir, exist_ok=True)
    log_format = "{time:YYYY-MM-DD HH:mm:ss} [{level}] {message}"
    logfile = os.path.join(outdir, 'ouroboros.log')
    logger.add(logfile, format=log_format, level='INFO')
    logger.add(sys.stderr, format=log_format, level='ERROR')
    logger.info(f"You are using Ouroboros version {__version__}")
    logger.info("You ran: " + ' '.join(sys.argv))


def short_reads_polish(assembly, short_1, short_2, outdir, depth, threads):
    polypolish_asm = run_polypolish(assembly, short_1, short_2, outdir, threads)

    polca_dir = os.path.join(outdir, 'pypolca')
    polca_asm = os.path.join(polca_dir, 'pypolca_corrected.fasta')
    run_pypolca(polypolish_asm, short_1, short_2, polca_dir, depth, threads)

    syscall(f"seqkit sort -r -l -o {os.path.join(outdir, '5_pypolca.fasta')} {polca_asm}")
    shutil.copy(os.path.join(polca_dir, 'pypolca.report'), outdir)
    shutil.rmtree(polca_dir)


@logger.catch
def check_plassembler_db_installation(db_dir):
    syscall(f"plassembler download -f -d {db_dir}")


def main():
    parser = argparse.ArgumentParser(
        prog="ouroboros.py",
        description='Ouroboros is a pipeline of hybrid assembly for bacteria.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    required = parser.add_argument_group("Required")
    required.add_argument('-i', '--infile', required=True, help='Input Nanopore FASTQ')
    required.add_argument('-1', '--short_1', required=True, help='Read 1 FASTQ to use for polishing')
    required.add_argument('-2', '--short_2', required=True, help='Read 2 FASTQ to use for polishing')
    required.add_argument('-o', '--outdir', required=True, help='Specify directory in which output has to be created.')

    optional = parser.add_argument_group("Optional")
    optional.add_argument(
        '-t', '--num-threads', default=1, type=int, help='Set the allowed number of threads to be used by the script'
    )
    optional.add_argument('-l', '--min-length', metavar='', default=1, type=int, help='Minimum length')
    optional.add_argument('-q', '--min-quality', metavar='', default=1, type=int, help='Minimum mean quality')
    optional.add_argument(
        '-p', '--keep-percent', metavar='', default=90, type=int, help='keep only this percentage of the best reads'
    )
    optional.add_argument('--nofilter', action='store_true', default=False, help='Disable long reads length filtering')
    optional.add_argument('--notrim', action='store_true', default=False, help='Disable long reads trimming.')
    optional.add_argument(
        '-x', '--depth', default=100, type=int, help='Sub-sample reads to this depth. Disable with --depth 0'
    )
    optional.add_argument('-g', '--gsize', default=None, metavar='', help='Estimated genome size eg. 3.2M <blank=AUTO>')
    optional.add_argument(
        '--medaka_model', default='r1041_e82_400bps_sup_v4.3.0', help='The model to be used by Medaka'
    )
    optional.add_argument('--hq', action='store_true', help="Flye will use '--nano-hq' instead of --nano-raw")
    optional.add_argument(
        '--contaminants', help='Contaminants FASTA file or Minimap2 index to map long reads against to filter out.'
    )
    optional.add_argument('-v', '--version', action='version', version=__version__)
    args = parser.parse_args()

    if (args.short_1 and not args.short_2) or (not args.short_1 and args.short_2):
        parser.error("-1 and -2 must be used together")

    begin(args.outdir)
    check_dependency()
    plassembler_db = os.path.join(__location__, 'plassembler_db')
    check_plassembler_db_installation(plassembler_db)
    validate_medaka_model(args.medaka_model)
    validate_fastq(args.infile)
    validate_fastq(args.short_1)
    validate_fastq(args.short_2)

    reads = args.infile
    if not args.notrim:
        logger.info("Remove long reads adapter and barcode sequence")
        trimmed_reads = os.path.join(args.outdir, 'READS_trim.fastq.gz')
        trim_long_reads(reads, trimmed_reads, args.num_threads)
        reads = trimmed_reads
    if not args.nofilter:
        first_filter = os.path.join(args.outdir, 'READS.qual.fastq.gz')
        quality_filter(
            reads, first_filter, keep_percent=args.keep_percent, min_quality=args.min_quality, threads=args.num_threads
        )
        reads = first_filter

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
    total_bases = fastq_stats['bases']

    if args.gsize:
        gsize = parse_genome_size(args.gsize)
        if gsize:
            logger.info(f"Using genome size was {gsize}bp.")
        else:
            logger.warning(f"Couldn't parse string {args.gsize}, will auto estimate.")
            gsize = estimate_genome_size(reads, args.num_threads)
    else:
        gsize = estimate_genome_size(reads, args.num_threads)

    origin_depth = total_bases / gsize
    logger.info(f'Estimated long sequencing depth: {origin_depth:.0f} x')
    if args.depth:
        if origin_depth > args.depth:
            logger.info(f"Subsampling reads from {origin_depth:.0f}x to {args.depth}x.")
            sub_reads = os.path.join(args.outdir, 'READS.sub.fastq')
            syscall(f"rasusa reads -b {gsize * args.depth} {reads} -o {sub_reads}")
        else:
            logger.info("No read depth reduction requested or necessary.")
            sub_reads = reads
    else:
        sub_reads = reads

    second_filter = length_filter(
        sub_reads, args.outdir, args.min_length, args.num_threads
    )

    logger.info('Trimming short reads.')
    trimmed_one = os.path.join(args.outdir, 'READS_1.fastq.gz')
    trimmed_two = os.path.join(args.outdir, 'READS_2.fastq.gz')
    trim_short_reads(args.short_1, args.short_2, trimmed_one, trimmed_two, args.num_threads)
    total_bases = fastq_scan(trimmed_one)['bases'] + fastq_scan(trimmed_one)['bases']

    origin_depth = total_bases / gsize
    logger.info(f'Estimated short sequencing depth: {origin_depth:.0f} x')

    logger.info("Assembling reads with Flye")
    flye_dir = os.path.join(args.outdir, 'flye')
    run_flye(second_filter, flye_dir, args.num_threads, args.hq)
    shutil.copyfile(os.path.join(flye_dir, 'flye.log'), os.path.join(args.outdir, 'flye.log'))
    shutil.copyfile(os.path.join(flye_dir, 'assembly_info.txt'), os.path.join(args.outdir, 'flye_info.txt'))
    shutil.copyfile(os.path.join(flye_dir, 'assembly_graph.gfa'), os.path.join(args.outdir, 'flye-unpolished.gfa'))

    logger.info('Reorients complete circular sequence.')
    dnaapler_dir = os.path.join(args.outdir, 'dnaapler')
    dnaapler_asm = os.path.join(dnaapler_dir, 'dnaapler_reoriented.fasta')
    run_dnaapler(flye_output=flye_dir, outdir=dnaapler_dir, threads=args.num_threads)
    syscall(f"seqkit sort -l -r {dnaapler_asm} -o {os.path.join(args.outdir, '1_flye.fasta')}")

    logger.info("Polishing with medaka.")
    medaka_dir = os.path.join(args.outdir, 'medaka')
    medaka_asm = os.path.join(args.outdir, '2_medaka.fasta')
    run_medaka(
        assembly=dnaapler_asm, reads=reads, outdir=medaka_dir, model=args.medaka_model, threads=args.num_threads
    )
    syscall(f"seqkit sort -l -r {os.path.join(medaka_dir, 'consensus.fasta')} -o {medaka_asm}")

    logger.info("Plasmid assembly.")
    plassembler_dir = os.path.join(args.outdir, 'plassembler')
    chromosome, plasmids = run_plassembler(
        reads=sub_reads,
        assembly=medaka_asm, assembly_info=os.path.join(args.outdir, 'flye_info.txt'),
        outdir=plassembler_dir, database=plassembler_db, short_one=trimmed_one, short_two=trimmed_two,
        threads=args.num_threads
    )
    plassembler_asm = os.path.join(args.outdir, '3_plassembler.fasta')
    syscall(f"cat {chromosome} {plasmids} > {plassembler_asm}")



    short_reads_polish(
        plassembler_asm, trimmed_one, trimmed_two, args.outdir, origin_depth, args.num_threads
    )
    syscall(f"seqkit sort -l -r  -o {os.path.join(args.outdir, 'assembly.fasta')} "
            f"{os.path.join(args.outdir, '5_pypolca.fasta')}")

    for d in (flye_dir, medaka_dir, dnaapler_dir):
        shutil.rmtree(d)
    syscall(f"rm {os.path.join(args.outdir, 'READS*')}")
    logger.info("Done")


if __name__ == '__main__':
    main()
