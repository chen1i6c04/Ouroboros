import os
import re
import sys
import json
import shutil
import argparse
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

def reads_filter(input_file, output_file, min_length=0, min_quality=0, keep_percent=90):
    if min_length or min_quality:
        logger.info(f"Filter out reads length less than {min_length} and average quality score less {min_quality}")
        cmd = f'filtlong --min_length {min_length} --min_mean_q {min_quality} {input_file}'
    else:
        logger.info(f"Keep only {keep_percent} percentage of the best reads")
        cmd = f'filtlong --keep_percent {keep_percent} {input_file}'
    cmd += f' | pigz > {output_file}'
    syscall(cmd)


def subsampling(infile, outfile, gsize, depth):
    bases = gsize * depth
    cmd = ['rasusa', '-b', str(bases), '-i', infile, '-o', outfile]
    syscall(cmd)


@logger.catch
def run_polca(assembly, short_1, short_2, output_dir, num_threads):
    """
    POLCA is a polishing tool in MaSuRCA (Maryland Super Read Cabog Assembler)
    https://github.com/alekseyzimin/masurca#polca
    """
    cmd = ['pypolca', 'run', '-a', assembly, '-1', short_1, '-2', short_2, '-t', str(num_threads), '-o', output_dir]
    syscall(cmd)


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
    parser.add_argument('-i', '--infile', metavar='',
                        required=True,
                        help='Input Nanopore FASTQ')
    parser.add_argument('-o', '--outdir', metavar='',
                        required=True,
                        help='Specify directory in which output has to be created.')
    parser.add_argument('-t', '--num-threads', metavar='',
                        default=1, type=int,
                        help='Set the allowed number of threads to be used by the script (default: 1)')
    parser.add_argument('-l', '--min-length', metavar='',
                        default=0, type=int,
                        help='')
    parser.add_argument('-q', '--min-quality', metavar='',
                        default=0, type=int,
                        help='')
    parser.add_argument('--nofilter', action='store_true', default=False,
                        help='Disable read length filtering')
    parser.add_argument('-1', '--short_1', metavar='',
                        help='Read 1 FASTQ to use for polishing')
    parser.add_argument('-2', '--short_2', metavar='',
                        help='Read 2 FASTQ to use for polishing')
    parser.add_argument('-x', '--depth', metavar='',
                        default=50, type=int,
                        help='Sub-sample reads to this depth. Disable with --depth 0 (default: 50)')
    parser.add_argument('-g', '--gsize', default=None, metavar='',
                        help='Estimated genome size eg. 3.2M <blank=AUTO> (default: "")')
    parser.add_argument('--meta', action='store_true',
                        help='Metagenome / uneven coverage')
    parser.add_argument('--medaka_model', default='r1041_e82_400bps_sup_g615', metavar='',
                        help='The model to be used by Medaka (default: r1041_e82_400bps_sup_g615)')
    parser.add_argument('--medaka_opt', metavar='', nargs='+',
                        help='Additional options to be given to Medaka')
    parser.add_argument('--hq', action='store_true',
                        help="Flye will use '--nano-hq' instead of --nano-raw")
    parser.add_argument('--contaminants', metavar='',
                        help=' Contaminants FASTA file or Minimap2 index to map long reads against to filter out.')
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

    if not os.access(args.infile, os.R_OK):
        logger.error(f"{args.infile} is not a readable file")
        sys.exit("Aborted")
    if args.infile.endswith('.fastq.gz'):
        raw_reads = os.path.join(args.outdir, 'READS.fastq.gz')
    elif args.infile.endswith('.fastq'):
        raw_reads = os.path.join(args.outdir, 'READS.fastq')
    else:
        sys.exit()

    os.symlink(args.infile, raw_reads)
    interm_reads = raw_reads
    if not args.nofilter:
        filtered_reads = os.path.join(args.outdir, 'READS_filter.fastq.gz')
        reads_filter(interm_reads, filtered_reads, args.min_length, args.min_quality)
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
    input_type = '--nano-hq' if args.hq else '--nano-raw'
    flye_dir = os.path.join(args.outdir, 'flye')
    flye_asm = os.path.join(flye_dir, 'assembly.fasta')
    flye_info = os.path.join(flye_dir, 'assembly_info.txt')
    flye_cmd = ['flye', input_type, final_reads, '-o', flye_dir, '-t', str(args.num_threads)]
    if args.meta:
        flye_cmd += ['--meta']
    logger.info(f"Flye command: {' '.join(flye_cmd)}")
    p = syscall(flye_cmd)
    if p.returncode:
        logger.error("Assembly failed.")
        sys.exit('Aborted')

    shutil.copyfile(
        flye_asm,
        os.path.join(args.outdir, '1_flye.fasta')
    )
    shutil.copyfile(
        flye_info,
        os.path.join(args.outdir, 'flye_info.txt')
    )
    shutil.copyfile(
        os.path.join(flye_dir, 'assembly_graph.gfa'),
        os.path.join(args.outdir, 'flye-unpolished.gfa')
    )
    shutil.copyfile(
        os.path.join(flye_dir, 'flye.log'),
        os.path.join(args.outdir, 'flye.log')
    )

    logger.info("Polishing with medaka.")
    medaka_dir = os.path.join(args.outdir, 'medaka')
    medaka_asm = os.path.join(medaka_dir, 'consensus.fasta')
    cmd = ['medaka_consensus', '-i', interm_reads, '-d', flye_asm, '-o', medaka_dir, '-m', args.medaka_model, '-t', str(args.num_threads)]
    if args.medaka_opt:
        cmd += args.medaka_opt
    logger.info(' '.join(cmd))
    syscall(cmd)
    shutil.copyfile(medaka_asm, os.path.join(args.outdir, '2_medaka.fasta'))

    if short_reads_polishing:
        logger.info("Polishing with short reads")
        polypolish_output = run_polypolish(medaka_asm, args.short_1, args.short_2, args.outdir, args.num_threads)
        polca_dir = os.path.join(args.outdir, 'polca')
        run_polca(polypolish_output, args.short_1, args.short_2, polca_dir, args.num_threads)
        shutil.copyfile(
            os.path.join(polca_dir, 'polca_corrected.fasta'),
            os.path.join(args.outdir, '4_polca.fasta')
        )
        shutil.copy(
            os.path.join(polca_dir, 'polca.report'),
            args.outdir
        )
        shutil.rmtree(polca_dir)
    shutil.rmtree(flye_dir)
    shutil.rmtree(medaka_dir)
    syscall(f"rm {os.path.join(args.outdir, 'READS*')}")
    logger.info("Done")


if __name__ == '__main__':
    main()
