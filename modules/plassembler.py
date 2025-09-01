import os
from tempfile import gettempdir
from loguru import logger
from .utils import syscall


@logger.catch
def run_plassembler(reads, flye_directory, outdir, database, short_one, short_two, threads):
    cmd = (f"plassembler run --skip_qc --skip_mash --depth_filter 0.5 "
           f"-t {threads} -d {database} -l {reads} "
           f"-1 {short_one} -2 {short_two} "
           f"-o {outdir} --flye_directory {flye_directory} "
           f"--spades_options '--tmp-dir {gettempdir()}' --unicycler_options '--keep 0'")
    logger.info(f"Running : {cmd}")
    syscall(cmd)


@logger.catch
def check_plassembler_db_installation(db_dir):
    mash_db_names = ['plsdb_2023_11_03_v2.msh', 'plsdb_2023_11_03_v2.tsv']
    f1 = os.path.join(db_dir, mash_db_names[0])
    f2 = os.path.join(db_dir, mash_db_names[1])
    if os.path.exists(db_dir) is False:
        syscall(f"plassembler download -d {db_dir}")
    if os.path.exists(f1) and os.path.exists(f2):
        logger.info(f"PLSDB Database at {db_dir} has already been downloaded")
    else:
        syscall(f"plassembler download -f -d {db_dir}")
