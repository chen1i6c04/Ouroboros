from .utils import syscall


def trim_long_reads(infile, outfile, threads):
    cmd = f"porechop -i {infile} -o {outfile} --check_reads 1000 --discard_middle --format fastq.gz -t {threads}"
    syscall(cmd)


def trim_short_reads(input_1, input_2, output_1, output_2, threads):
    cmd = (f"fastp -i {input_1} -I {input_2} -o {output_1} -O {output_2} -l 50 -t 1 -5 -3 -w {threads} "
           f"--detect_adapter_for_pe -j /dev/null -h /dev/null")
    syscall(cmd)
