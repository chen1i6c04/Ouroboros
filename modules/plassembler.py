import os
from tempfile import gettempdir
from loguru import logger
from .utils import syscall


@logger.catch
def run_plassembler(reads, assembly, assembly_info, outdir, database, short_one, short_two, threads):
    cmd = (f"plassembler run --skip_qc --keep_chromosome -t {threads} -d {database} -l {reads} --depth_filter 0.8 "
           f"-1 {short_one} -2 {short_two} "
           f"-o {outdir} --flye_assembly {assembly} --flye_info {assembly_info} "
           f"--spades_options '--tmp-dir {gettempdir()}' --unicycler_options '--keep 0'")
    logger.info(f"Running : {cmd}")
    syscall(cmd)
    chromosome = os.path.join(outdir, 'chromosome.fasta')
    plasmids = os.path.join(outdir, 'plassembler_plasmids.fasta')
    return chromosome, plasmids
