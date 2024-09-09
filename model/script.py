from Bio import PDB
from Bio.SeqUtils import seq1

# Define file paths
pdb_files = [
    "C:\\Users\\evgen\\OneDrive\\Skrivbord\\model\\2ydv.pdb",
    "C:\\Users\\evgen\\OneDrive\\Skrivbord\\model\\P41968.ent"
]
pdb_info_file = "C:\\Users\\evgen\\OneDrive\\Skrivbord\\model\\pdb_info.txt"
fasta_file = "C:\\Users\\evgen\\OneDrive\\Skrivbord\\model\\fasta_proteins.txt"
python_proteins_file = "C:\\Users\\evgen\\OneDrive\\Skrivbord\\model\\python_proteins.txt"

def parse_pdb(file_path, file_handle):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', file_path)
    
    file_handle.write(f"Parsing file: {file_path}\n")
    file_handle.write(f"Number of models: {len(structure)}\n")
    for model in structure:
        file_handle.write(f"Model ID: {model.id}\n")
        for chain in model:
            file_handle.write(f"  Chain ID: {chain.id}\n")
            for residue in chain:
                file_handle.write(f"    Residue: {residue.get_resname()} {residue.get_id()}\n")
                for atom in residue:
                    file_handle.write(f"      Atom: {atom.get_name()}\n")
    file_handle.write("\n")

def extract_fasta_from_pdb(file_path, fasta_handle, python_handle, protein_set):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('structure', file_path)
    
    for model in structure:
        for chain in model:
            # Extract sequence and convert to FASTA format
            sequence = ""
            for residue in chain:
                if PDB.is_aa(residue):
                    sequence += seq1(residue.get_resname())
            
            if sequence:
                # Create a dummy UniProt ID and description for illustration
                uniprot_id = file_path.split("\\")[-1].split(".")[0]  # Using filename as dummy ID
                description = f"sp|{uniprot_id}|PROT_HUMAN Example protein description OS=Homo sapiens OX=9606 GN=EXAMPLE PE=1 SV=1"
                
                if uniprot_id not in protein_set:
                    # Write to FASTA file
                    fasta_handle.write(f">{description}\n")
                    for i in range(0, len(sequence), 60):
                        fasta_handle.write(f"{sequence[i:i+60]}\n")
                    fasta_handle.write("\n")  # Add space between sequences
                    
                    # Write to python_proteins.txt in alignment format
                    python_handle.write(f">{description}\n")
                    python_handle.write(f"{sequence}\n")
                    python_handle.write("*\n\n")  # Add space between sequences
                    
                    # Add to protein set for the header
                    protein_set.add(uniprot_id)

# Open the files and perform the parsing and extraction
protein_set = set()
with open(pdb_info_file, 'w') as pdb_file_handle, \
     open(fasta_file, 'w') as fasta_file_handle, \
     open(python_proteins_file, 'w') as python_file_handle:

    for pdb_file in pdb_files:
        parse_pdb(pdb_file, pdb_file_handle)
        extract_fasta_from_pdb(pdb_file, fasta_file_handle, python_file_handle, protein_set)

# Adding a header to the FASTA file
with open(fasta_file, 'r+') as file_handle:
    content = file_handle.read()
    file_handle.seek(0, 0)
    file_handle.write("Homology Model\n\nProteins for Alignment\n")
    for i, protein_id in enumerate(sorted(protein_set), start=1):
        file_handle.write(f"{i}. UniProt ID: {protein_id}\n")
    file_handle.write("\n")
    file_handle.write(content)

# Adding a header to the Python-compatible file
with open(python_proteins_file, 'r+') as file_handle:
    content = file_handle.read()
    file_handle.seek(0, 0)
    file_handle.write("Homology Model\n\nProteins for Alignment\n")
    for i, protein_id in enumerate(sorted(protein_set), start=1):
        file_handle.write(f"{i}. UniProt ID: {protein_id}\n")
    file_handle.write("\n")
    file_handle.write(content)

print(f"Information has been written to {pdb_info_file}")
print(f"FASTA sequences have been written to {fasta_file}")
print(f"Python-compatible sequences for alignment have been written to {python_proteins_file}")
