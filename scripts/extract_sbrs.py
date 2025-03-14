import pysam
from tqdm import tqdm


def reverse_complement(dna_sequence: str) -> str:
    """
    Generate the reverse complement of a DNA sequence.

    Args:
        dna_sequence (str): The original DNA sequence.

    Returns:
        str: The reverse complement of the DNA sequence.
    """
    # Define the complement mapping
    complement = str.maketrans("ATCG", "TAGC")

    # Get the complement sequence
    complement_sequence = dna_sequence.translate(complement)

    # Reverse the complement sequence
    reverse_complement_sequence = complement_sequence[::-1]

    return reverse_complement_sequence


def extact_channel_sbrs(bam_file, channel):

    seqs = []
    rc = False
    with pysam.AlignmentFile(bam_file, "rb", check_sq=False, threads=40) as bam:
        for read in tqdm(bam.fetch(until_eof=True)):
            if read.has_tag("ch"):
                if read.get_tag("ch") == channel:
                    seq = read.query_sequence
                    if rc:
                        seq = reverse_complement(seq)
                    rc = not rc
                    seqs.append(seq)
    print("\n\n".join(seqs))


def main():
    channel = 332939
    bam_file = "/data/ccs_data/case-study/20250310-lowQ30/Output.subreads.bam"
    extact_channel_sbrs(bam_file, channel)


if __name__ == "__main__":
    main()
