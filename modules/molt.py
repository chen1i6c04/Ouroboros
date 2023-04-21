import os
import argparse
import shutil
from tempfile import TemporaryDirectory
from modules.utils import syscall


def polca_polish(assembly, short_reads_1, short_reads_2, output_dir, num_threads):
    """
    POLCA is a polishing tool in MaSuRCA (Maryland Super Read Cabog Assembler)
    https://github.com/alekseyzimin/masurca#polca
    """
    assembly = os.path.abspath(assembly)
    os.chdir(output_dir)
    cmd = f"polca.sh -a {assembly} -r '{short_reads_1} {short_reads_2}' -t {num_threads} 2>&1 | tee /dev/tty"
    syscall(cmd)


def polypolish_polish(assembly, short_reads_1, short_reads_2, output_dir, num_threads):
    alignments_1, alignments_2 = read_alignments(assembly, short_reads_1, short_reads_2, output_dir, num_threads)
    filtered_1 = os.path.join(output_dir, 'filtered_1.sam')
    filtered_2 = os.path.join(output_dir, 'filtered_2.sam')
    polished_assembly = os.path.join(output_dir, 'polypolish.fasta')
    cmd = f"polypolish_insert_filter.py --in1 {alignments_1} --in2 {alignments_2} --out1 {filtered_1} --out2 {filtered_2} " \
          f"2>&1 | tee /dev/tty"
    syscall(cmd)
    cmd = f"polypolish {assembly} {filtered_1} {filtered_2} 2>&1 1> {polished_assembly} | " \
          f"tee /dev/tty"
    syscall(cmd)


def read_alignments(assembly, short_reads_1, short_reads_2, output_dir, num_threads):
    alignments_1 = os.path.join(output_dir, 'alignments_1.sam')
    alignments_2 = os.path.join(output_dir, 'alignments_2.sam')
    syscall(f"bwa-mem2 index {assembly}")
    syscall(f"bwa-mem2 mem -t {num_threads} -a {assembly} {short_reads_1} > {alignments_1}")
    syscall(f"bwa-mem2 mem -t {num_threads} -a {assembly} {short_reads_2} > {alignments_2}")
    return alignments_1, alignments_2


def molt(assembly, short_reads_1, short_reads_2, output_file, num_threads):
    with TemporaryDirectory() as tmp_dir:
        polypolish_assembly = os.path.join(tmp_dir, 'polypolish.fasta')
        alignments_1, alignments_2 = read_alignments(
            assembly, short_reads_1, short_reads_2, tmp_dir, num_threads
        )
        polypolish_polish(assembly, alignments_1, alignments_2, tmp_dir, num_threads)

        polca_assembly = os.path.join(tmp_dir, 'polypolish.fasta.PolcaCorrected.fa')
        polca_polish(polypolish_assembly, short_reads_1, short_reads_2, tmp_dir, num_threads)
        shutil.copyfile(polca_assembly, output_file)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('assembly')
    parser.add_argument('--sr1', required=True)
    parser.add_argument('--sr2', required=True)
    parser.add_argument('-o', '--output-file', required=True)
    parser.add_argument('-t', '--num-threads', default=1, type=int)
    parser.add_argument('--tmp-dir', default='/tmp')
    args = parser.parse_args()

    with TemporaryDirectory(dir=args.tmp_dir) as tmp_dir:
        polypolish_assembly = os.path.join(tmp_dir, 'polypolish.fasta')
        alignments_1, alignments_2 = read_alignments(
            args.assembly, args.sr1, args.sr2, tmp_dir, args.num_threads
        )
        polypolish_polish(args.assembly, alignments_1, alignments_2, tmp_dir)

        polca_assembly = os.path.join(tmp_dir, 'polypolish.fasta.PolcaCorrected.fa')
        polca_polish(polypolish_assembly, args.sr1, args.sr2, tmp_dir, args.num_threads)
        shutil.copyfile(polca_assembly, args.output_file)


if __name__ == '__main__':
    main()
