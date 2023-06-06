from pathlib import Path
from io import StringIO
import pandas as pd
from loguru import logger
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastxCommandline


current_location = Path(__file__)
gene_data = current_location.parent.parent/'gene_data'
gene_dict = {
    'dnaA': gene_data/'dnaA_db',
    'repA': gene_data/'repA_db'
}


def run_blastx(query, db, num_threads):
    columns = [
        "qseqid", "qlen", "sseqid", "slen", "length", "qstart", "qend", "sstart", "send",
        "pident", "nident", "gaps", "mismatch", "evalue", "bitscore", "qseq", "sseq",
    ]
    cline = NcbiblastxCommandline(
        query=query, db=db, evalue=1e-10, num_threads=num_threads, outfmt=f"6 {' '.join(columns)}"
    )
    stdout, stderr = cline()
    return pd.read_csv(StringIO(stdout), sep='\t', nrows=1, names=columns)


def reorient_sequence(query, anchor, strand='forward'):
    record = SeqIO.read(query, 'fasta')
    length = len(record)
    if strand == "forward":
        edge = anchor - 1
        left = record.seq[edge:length]
        right = record.seq[0:edge]
        record.seq = left + right
    if strand == "reverse":
        record.seq = record.seq.reverse_complement()
        edge = length - anchor
        left = record.seq[edge:length]
        right = record.seq[0:edge]
        record.seq = left + right
    SeqIO.write(record, query, 'fasta')


def rotate_circular_sequences(query, threads):
    for gene_name, db in gene_dict.items():
        blast_out = run_blastx(query, db, threads)
        if blast_out.empty:
            pass
        else:
            blast_hit = blast_out.iloc[0].to_dict()
            if (blast_hit['qseq'][0] not in ("M", "V", "L")) or blast_hit['sstart'] != 1:
                pass
            elif blast_hit['qstart'] == 1:
                pass
            else:
                logger.info(f"{gene_name} identified in {query.stem}. It starts at coordinate {blast_hit['qstart']}.")
                strand = "reverse" if blast_hit['qstart'] > blast_hit['qend'] else "forward"
                reorient_sequence(query, blast_hit['qstart'], strand)
                break


def parse_flye_info(flye_info):
    df = pd.read_csv(flye_info, sep='\t', index_col='#seq_name')
    return df.to_dict('index')


def reorient_assembly(assembly, assembly_info, outdir, threads):
    outdir = Path(outdir)
    outdir.mkdir()

    assembly_info = parse_flye_info(assembly_info)

    for record in SeqIO.parse(assembly, 'fasta'):
        outfile = outdir/(record.id + '.fa')
        SeqIO.write(record, outfile, 'fasta')
    for contig in outdir.iterdir():
        if assembly_info[contig.stem]['circ.'] == 'Y':
            rotate_circular_sequences(contig, threads)
    records = [SeqIO.read(contig, 'fasta') for contig in outdir.iterdir()]
    SeqIO.write(records, outdir/'assembly.fasta', 'fasta')

