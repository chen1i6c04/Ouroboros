import os
import re
import sys
import json
import shutil
import argparse
from tabulate import tabulate
from loguru import logger
from modules.utils import syscall, medaka_model_check, estimate_genome_size, read_alignments
from modules.assembly_process import reorient_assembly


current_location = os.path.dirname(os.path.abspath(__file__))
bin_path = os.path.join(current_location, 'bin')


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
        'masurca': f'{os.path.join(bin_path, "masurca")} --version',
        'bwa': 'bwa 2>&1 | grep Version:',
    }
    for program_name, cmd in version.items():
        p = syscall(cmd)
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


def reads_filter(raw_reads, filtered_reads, min_length=0, min_quality=0):
    if min_length or min_quality:
        logger.info(f"Filter out reads length less than {min_length} and average quality score less {min_quality}")
        cmd = f'nanoq -l {min_length} -q {min_quality} -i {raw_reads} > {filtered_reads}'
    else:
        logger.info("Keep only 90 percentage of the best reads")
        cmd = f'filtlong --keep_percent 90 --mean_q_weight 10 {raw_reads} > {filtered_reads}'
    syscall(cmd)


def subsampling(infile, outfile, gsize, depth):
    bases = gsize * depth
    cmd = ['rasusa', '-b', str(bases), '-i', infile, '-o', outfile]
    syscall(cmd)


def nanoq_stats(seqfile):
    cmd = ['nanoq', '-s', '-j', '-i', seqfile]
    p = syscall(cmd)
    result = json.loads(p.stdout)
    return result


@logger.catch
def run_polca(assembly, short_reads_1, short_reads_2, output_dir, num_threads):
    """
    POLCA is a polishing tool in MaSuRCA (Maryland Super Read Cabog Assembler)
    https://github.com/alekseyzimin/masurca#polca
    """
    assembly = os.path.abspath(assembly)
    program = os.path.join(bin_path, 'polca.sh')
    os.makedirs(output_dir, exist_ok=True)
    os.chdir(output_dir)
    cmd = f"{program} -a {assembly} -r '{short_reads_1} {short_reads_2}' -t {num_threads}"
    syscall(cmd)
    os.remove(assembly + '.fai')


def run_polypolish(assembly, short_reads_1, short_reads_2, output_dir, num_threads):
    alignments_1 = os.path.join(output_dir, 'alignments_1.sam')
    alignments_2 = os.path.join(output_dir, 'alignments_2.sam')
    filtered_1 = os.path.join(output_dir, 'filtered_1.sam')
    filtered_2 = os.path.join(output_dir, 'filtered_2.sam')
    polypolish_output = os.path.join(output_dir, '3_polypolish.fasta')
    read_alignments(assembly, short_reads_1, alignments_1, num_threads)
    read_alignments(assembly, short_reads_2, alignments_2, num_threads)
    cmd = f"polypolish_insert_filter.py " \
          f"--in1 {alignments_1} --in2 {alignments_2} --out1 {filtered_1} --out2 {filtered_2}"
    syscall(cmd)
    cmd = f"polypolish {assembly} {filtered_1} {filtered_2} > {polypolish_output}"
    syscall(cmd)
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
    parser.add_argument('--sr1', metavar='',
                        help='Read 1 FASTQ to use for polishing')
    parser.add_argument('--sr2', metavar='',
                        help='Read 2 FASTQ to use for polishing')
    parser.add_argument('--depth', metavar='',
                        default=100, type=int,
                        help='Sub-sample reads to this depth. Disable with --depth 0 (default: 100)')
    parser.add_argument('--gsize', default=None, metavar='',
                        help='Estimated genome size eg. 3.2M <blank=AUTO> (default: "")')
    parser.add_argument('--meta', action='store_true',
                        help='Metagenome / uneven coverage')
    parser.add_argument('--model', default='r941_min_hac_g507', metavar='',
                        help='The model to be used by Medaka (default: r941_min_hac_g507)')
    parser.add_argument('--hq', action='store_true',
                        help="Flye will use '--nano-hq' instead of --nano-raw")
    args = parser.parse_args()

    short_reads_polishing = False
    if (args.sr1 and not args.sr2) or (not args.sr1 and args.sr2):
        parser.error("--sr1 and --sr2 must be used together")
    elif args.sr1 and args.sr2:
        short_reads_polishing = True
    else:
        pass

    initialize(args)

    if medaka_model_check(args.model) is False:
        logger.error(f"Medaka model {args.model} unavailable")
        sys.exit("Aborted")

    if args.infile.endswith('.fastq.gz'):
        reads = os.path.join(args.outdir, 'READS.fastq.gz')
    elif args.infile.endswith('.fastq'):
        reads = os.path.join(args.outdir, 'READS.fastq')
    else:
        sys.exit()
    os.symlink(args.infile, reads)
    stats_before_filter = nanoq_stats(reads)
    filt_reads = os.path.join(args.outdir, 'READS.flit.fastq')
    reads_filter(reads, filt_reads, args.min_length, args.min_quality)
    reads = filt_reads
    stats_after_filter = nanoq_stats(reads)

    index = ['reads', 'bases', 'n50', 'mean_length', 'median_length', 'mean_quality', 'median_quality']
    data = [[i, stats_before_filter[i], stats_after_filter[i]] for i in index]
    comparison = tabulate(data, headers=['', 'before filter', 'after filter'], floatfmt='.0f')
    logger.info(f"""\n\n{comparison}\n""")

    if not args.meta:
        total_bases = stats_after_filter['bases']
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
                reads = sub_reads
            else:
                logger.info("No read depth reduction requested or necessary.")
    else:
        logger.info("Flag '--meta' was be set, won't subsampling reads.")

    logger.info("Assembling reads with Flye")
    input_type = '--nano-hq' if args.hq else '--nano-raw'
    flye_dir = os.path.join(args.outdir, 'flye')
    flye_asm = os.path.join(flye_dir, 'assembly.fasta')
    flye_cmd = ['flye', input_type, reads, '-o', flye_dir, '-t', str(args.num_threads)]
    if args.meta:
        flye_cmd += ['--meta']
    logger.info(f"Flye command: {' '.join(flye_cmd)}")
    p = syscall(flye_cmd)
    if p.returncode:
        logger.error("Assembly failed.")
        sys.exit('Aborted')

    logger.info("Rotate circular sequences.")
    reorient_dir = os.path.join(args.outdir, 'reorient')
    reoriented_assembly = os.path.join(reorient_dir, 'assembly.fasta')
    reorient_assembly(
        flye_asm,
        os.path.join(flye_dir, 'assembly_info.txt'),
        reorient_dir,
        args.num_threads
    )

    shutil.copyfile(
        reoriented_assembly,
        os.path.join(args.outdir, '1_flye.fasta')
    )
    shutil.move(
        os.path.join(flye_dir, 'assembly_info.txt'),
        os.path.join(args.outdir, 'flye_info.txt')
    )
    shutil.move(
        os.path.join(flye_dir, 'assembly_graph.gfa'),
        os.path.join(args.outdir, 'flye-unpolished.gfa')
    )
    shutil.move(
        os.path.join(flye_dir, 'flye.log'),
        os.path.join(args.outdir, 'flye.log')
    )

    logger.info("Polishing with medaka.")
    medaka_dir = os.path.join(args.outdir, 'medaka')
    medaka_asm = os.path.join(medaka_dir, 'consensus.fasta')
    cmd = ['medaka_consensus', '-i', filt_reads, '-d', reoriented_assembly, '-o', medaka_dir, '-m', args.model, '-t', str(args.num_threads)]
    syscall(cmd)

    shutil.copyfile(medaka_asm, os.path.join(args.outdir, '2_medaka.fasta'))

    if short_reads_polishing:
        logger.info("Polishing with short reads")
        polypolish_output = run_polypolish(medaka_asm, args.sr1, args.sr2, args.outdir, args.num_threads)
        polca_dir = os.path.join(args.outdir, 'polca')
        run_polca(polypolish_output, args.sr1, args.sr2, polca_dir, args.num_threads)
        shutil.copyfile(
            os.path.join(polca_dir, '3_polypolish.fasta.PolcaCorrected.fa'),
            os.path.join(args.outdir, '4_polca.fasta')
        )
        shutil.rmtree(polca_dir)
    shutil.rmtree(flye_dir)
    shutil.rmtree(medaka_dir)
    shutil.rmtree(reorient_dir)
    syscall(f"rm {os.path.join(args.outdir, 'READS.*')}")
    logger.info("Done")


if __name__ == '__main__':
    main()
