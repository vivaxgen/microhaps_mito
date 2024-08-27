from ngs_pipeline import cerr, cexit, check_NGSENV_BASEDIR, snakeutils
import pysam
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def init_argparser():
    p = snakeutils.init_argparser("Extracts the merged query sequence of paired reads from a BAM file")
    p.add_argument("--bam", help="Input BAM file")
    p.add_argument("--bed", help="BED file with the region of interest")
    p.add_argument("--fasta", help="Output FASTA file")
    return p

# follow bed format start = 0-based, end = 1-based
# start is exclusive, end is inclusive
def get_paired_merged_query_sequence(start, end, read1, read2):
    result = []
    reads = [read1, read2]
    reads_pos_base = [{
        rpos: (r_1_2.query_sequence[qpos], r_1_2.query_qualities[qpos])
            for qpos, rpos in r_1_2.get_aligned_pairs(matches_only=True)} for r_1_2 in reads]
    
    for read_i, read in enumerate(reads):
        missing_rpos = set(range(read.reference_start, read.reference_end)) - set(reads_pos_base[read_i].keys())
        for pos in missing_rpos:
            reads_pos_base[read_i][pos] = ("-", 1) 
        
    for pos in range(start, end):
        temp = ('', 0)
        if pos in reads_pos_base[0].keys():
            temp = reads_pos_base[0][pos]
        if pos in reads_pos_base[1].keys():
            if reads_pos_base[1][pos][1] > temp[1]:
                temp = reads_pos_base[1][pos]
        if temp[0] == '':
            temp = ('N', 0)
        result.append(temp[0])
    return(''.join(result).upper())

def main(args):
    bamfile = pysam.AlignmentFile(args.bam, "rb")
    bed = pd.read_csv(args.bed, sep="\t", header=None)
    chrom = bed.iloc[0,0]
    roi_start = bed.iloc[0,1]
    roi_end = bed.iloc[0,2]
    read_names = []
    results = []
    for read in bamfile.fetch(chrom, roi_start, roi_end):
        # Check if read is paired
        if not read.is_paired:
            continue
        # Already processed
        if read.query_name in read_names:
            continue
        
        read_names.append(read.query_name)
        try:
            mate = bamfile.mate(read)
        except:
            cerr(f"Read {read.query_name} does not have a mate")
            continue
        results.append(
            SeqRecord(
                Seq(
                    get_paired_merged_query_sequence(roi_start, roi_end, read, mate)
                ),
                id=read.query_name,
                description=""
            )
        )
    SeqIO.write(results, args.fasta, "fasta")