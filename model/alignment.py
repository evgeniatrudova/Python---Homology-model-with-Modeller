from Bio import AlignIO
from io import StringIO

def main():
    # Define the file paths
    alignment_output_file_path = 'uniprot_alignment_fasta_proteins.txt'
    
    # Instructions:
    # 1. Perform the alignment of your FASTA protein sequences using UniProt's alignment tool.
    # 2. Save the resulting alignment in FASTA format.
    # 3. Replace the content of the `uniprot_alignment_fasta_proteins` variable with the FASTA alignment result.
    # 4. Run this script to save the formatted alignment to `uniprot_alignment_fasta_proteins.txt`.
    
    # Replace the content of uniprot_alignment_fasta_proteins with your alignment
    uniprot_alignment_fasta_proteins = """
>sp|P41968|MC3R_HUMAN
MNASCCLPSVQPTLPNGSEHLQAPFFSNQSSSAFCEQVFIKPEVFLSLGIVSLLENILVILAVVRNGNLHSPMYFFLCSLAVADMLVSVSNALETIMIAIVHSDYLTFEDQFIQHMD-NI-FDSMICISLVASICNLLAIAVDRYVTIFYALRYHSIMTVRKALTLIVAIWVCCGVCGVVFI--------------------------VYSE----SKMVIVCLITMFFAMMLLMGTLYVHMFLFARLHVKRIAALPPAD-GVAPQQHSCMKGAVTITILLGVFIFCWAPFFLHLVLIITCPTNPYCICYTAHFNTYLVLIMCNSVIDPLIYAFRSLELRNTFREILCGCNGMNLG----------------------------------------------------------------------------------------------------
>sp|P29274|AA2AR_HUMAN
----------------------MP--------IMGSSVYIT--VELAIAVLAILGNVLVCWAVWLNSNLQNVTNYFVVSLAAADIAVGVLAIPFAITI----------STGFCAACHGCLFIACFVLVLTQSSIFSLLAIAIDRYIAIRIPLRYNGLVTGTRAKGIIAICWVLSFAIGLTPMLGWNNCGQPKEGKNHSQGCGEGQVACLFEDVVPMNYMVYFNFFACVLVPLLLMLGVYLRIFLAARRQLKQMESQPLPGERARSTLQKEVHAAKSLAIIVGLFALCWLPLHIINCFTFFCPDCSHAPLWLM--YLAIVLSHTNSVVNPFIYAYRIREFRQTFRKIIRSHVLRQQEPFKAAGTSARVLAAHGSDGEQVSLRLNGHPPGVWANGSAPHPERRPNGYALGLVSGGSAQESQGNTGLPDVELLSHELKGVCPEPPGLDDPLAQDGAGVS
"""
    
    # Use StringIO to simulate a file object from the string
    alignment_file = StringIO(uniprot_alignment_fasta_proteins)
    
    # Read the FASTA alignment results
    alignment = AlignIO.read(alignment_file, 'fasta')
    
    # Convert the alignment to a readable format
    alignment_str = format_alignment(alignment)
    
    # Save the formatted alignment to a text file
    with open(alignment_output_file_path, 'w') as file:
        file.write(alignment_str)
    
    print(f"Alignment results saved to {alignment_output_file_path}")

def format_alignment(alignment):
    """
    Format the alignment into a readable string with positions and conserved residues.
    """
    aln_str = "CLUSTAL O(1.2.4) multiple sequence alignment\n\n"
    
    num_sequences = len(alignment)
    seq_length = len(alignment[0].seq)

    # Format the alignment into rows
    for i in range(0, seq_length, 60):  # Adjust width as needed
        # Print positions
        aln_str += " _aln.pos"
        for pos in range(i + 1, min(i + 61, seq_length + 1), 10):
            aln_str += f"{pos:>6}"
        aln_str += "\n"
        
        # Print each sequence
        for record in alignment:
            aln_str += f"{record.id:>17} {format_sequence(str(record.seq), i, 60)}\n"
        
        # Print conserved residues
        aln_str += " " * 18  # For alignment of the conservation line
        for pos in range(i, min(i + 60, seq_length)):
            columns = [seq[pos] for seq in alignment]
            if all(col == columns[0] for col in columns):
                aln_str += "*"
            elif any(col == '-' for col in columns):
                aln_str += "-"
            else:
                aln_str += " "
        aln_str += "\n\n"

    return aln_str

def format_sequence(seq, start, width):
    """
    Format the sequence to ensure proper alignment and fit within the given width.
    """
    formatted_seq = ""
    for i in range(start, min(start + width, len(seq)), 10):
        formatted_seq += f"{seq[i:i+10]} "
    return formatted_seq.strip()

if __name__ == "__main__":
    main()
