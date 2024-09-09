from Bio import PDB
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.stats import zscore
import os

def extract_ca_coordinates(pdb_file):
    """
    Extract Cα atom coordinates from a PDB file.
    """
    parser = PDB.PDBParser(QUIET=True)
    try:
        structure = parser.get_structure('protein', pdb_file)
    except Exception as e:
        print(f"Error parsing {pdb_file}: {e}")
        return np.array([])

    coordinates = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if PDB.is_aa(residue):
                    try:
                        ca = residue['CA']
                        coordinates.append(ca.get_coord())
                    except KeyError:
                        # Cα atom not present in this residue
                        continue
    return np.array(coordinates)

def calculate_distance_matrix(coords):
    """
    Calculate pairwise distance matrix for given coordinates.
    """
    if len(coords) < 2:
        return np.array([[]])
    return squareform(pdist(coords))

def calculate_sequential_angles(coords):
    """
    Calculate angles for sequential triplets of Cα atoms.
    """
    if len(coords) < 3:
        return np.array([])

    angles = []
    num_atoms = len(coords)

    for i in range(num_atoms - 2):
        a = coords[i]
        b = coords[i + 1]
        c = coords[i + 2]

        # Vectors
        ba = a - b
        bc = c - b

        # Calculate the angle using the dot product formula
        dot_product = np.dot(ba, bc)
        norm_ba = np.linalg.norm(ba)
        norm_bc = np.linalg.norm(bc)

        if norm_ba > 0 and norm_bc > 0:
            cos_angle = dot_product / (norm_ba * norm_bc)
            # Clamp the cosine to the valid range [-1, 1] to avoid numerical issues
            cos_angle = np.clip(cos_angle, -1.0, 1.0)
            angle = np.arccos(cos_angle)
            angles.append(np.degrees(angle))
    
    return np.array(angles)

def approximate_dope_score(coords):
    """
    Approximate the DOPE score based on pairwise distances.
    """
    distances = calculate_distance_matrix(coords)
    if distances.size == 0:
        return None
    mean_distance = np.mean(distances)
    std_distance = np.std(distances)
    return -mean_distance / (std_distance + 1e-6)  # Avoid division by zero

def approximate_ga341_score(coords):
    """
    Approximate the GA341 score based on sequential angles.
    """
    angles = calculate_sequential_angles(coords)
    if angles.size == 0:
        return None
    mean_angle = np.mean(angles)
    std_angle = np.std(angles)
    return -mean_angle / (std_angle + 1e-6)  # Avoid division by zero

def main():
    pdb_files = [
        'C:/Users/evgen/OneDrive/Skrivbord/model/2ydv.pdb',
        'C:/Users/evgen/OneDrive/Skrivbord/model/P41968.ent'
    ]

    # Prepare the output file
    output_file = 'validation_scores.txt'
    with open(output_file, 'w') as f:
        # Write the header
        f.write("PDB_File\tDOPE_Score\tGA341_Score\n")

        for pdb_file in pdb_files:
            print(f"Analyzing structure: {pdb_file}")

            if not os.path.exists(pdb_file):
                print(f"File not found: {pdb_file}")
                f.write(f"{os.path.basename(pdb_file)}\tError: File not found\tError: File not found\n")
                continue

            coords = extract_ca_coordinates(pdb_file)

            if len(coords) == 0:
                print(f"No Cα coordinates found for {pdb_file}.")
                f.write(f"{os.path.basename(pdb_file)}\tError: No Cα coordinates\tError: No Cα coordinates\n")
                continue

            dope_score = approximate_dope_score(coords)
            ga341_score = approximate_ga341_score(coords)

            # Handle cases where scores couldn't be calculated
            dope_score_str = f"{dope_score:.4f}" if dope_score is not None else "N/A"
            ga341_score_str = f"{ga341_score:.4f}" if ga341_score is not None else "N/A"

            print(f"DOPE Score for {pdb_file}: {dope_score_str}")
            print(f"GA341 Score for {pdb_file}: {ga341_score_str}")

            # Write the scores to the file
            f.write(f"{os.path.basename(pdb_file)}\t{dope_score_str}\t{ga341_score_str}\n")

    print(f"\nValidation scores have been written to {output_file}")

if __name__ == '__main__':
    main()
