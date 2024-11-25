import os
from loguru import logger
from .utils import syscall


def length_filter(infile, outdir, min_length=1, threads=1):
    logger.info(f"Filter out reads length less than {min_length:,}bp")
    outfile = os.path.join(outdir, 'READS.len.fastq.gz')
    if min_length == 1:
        outfile = infile
    else:
        syscall(f'nanoq -i {infile} -f -l {min_length} | pigz -6 -p {threads} > {outfile}')
    return outfile


def quality_filter(infile, outfile, keep_percent=90, min_quality=1, threads=1):
    logger.info(f"Keep only {keep_percent} percentage of the best reads and average quality score great than {min_quality}")
    if min_quality == 1:
        cmd = f'filtlong --keep_percent {keep_percent} {infile} | pigz -6 -p {threads} > {outfile}'
    else:
        cmd = f'filtlong --keep_percent {keep_percent} {infile} | nanoq -q {min_quality} | pigz -6 -p {threads} > {outfile}'
    syscall(cmd)
