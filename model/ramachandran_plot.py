"""
Ramachandran Plot Generator

This script generates a Ramachandran plot from protein structure files in PDB format.

Libraries required:
- biopython: For handling and parsing protein structure files.
- matplotlib: For plotting the Ramachandran plot.
- numpy: For numerical operations and angle conversions.

To install the required libraries, run the following commands in your terminal or command prompt:

pip install biopython
pip install matplotlib
pip install numpy

Usage:
1. Save your protein structure files in PDB format to the same directory as this script.
2. Update the 'pdb_files' list in the script with the paths to your PDB files.
3. Run the script using:
   python ramachandran_plot.py

The Ramachandran plot will be saved as 'ramachandran_plot.png' in the same directory as the script.
The log file with calculation details will be saved as 'ramachandran_plot_log.txt'.
"""

from Bio import PDB
import numpy as np
import matplotlib.pyplot as plt
import os

def extract_phi_psi_angles(structure):
    """
    Extract phi and psi angles from a protein structure.
    """
    phi_angles = []
    psi_angles = []

    def calc_phi_psi_angles(chain):
        for polypeptide in PDB.PPBuilder().build_peptides(chain):
            angles = polypeptide.get_phi_psi_list()
            for phi, psi in angles:
                if phi is not None and psi is not None:
                    phi_angles.append(np.degrees(phi))
                    psi_angles.append(np.degrees(psi))

    for model in structure:
        for chain in model:
            calc_phi_psi_angles(chain)
    
    return phi_angles, psi_angles

def plot_ramachandran(phi_angles, psi_angles, protein_names):
    """
    Plot the Ramachandran plot from phi and psi angles.
    """
    plt.figure(figsize=(8, 6))
    plt.scatter(phi_angles, psi_angles, s=10, alpha=0.5, color='blue', label='Residue Angles')
    
    plt.title('Ramachandran Plot')
    plt.suptitle(f'Proteins: {", ".join(protein_names)}', y=0.95)  # Subtitle with protein names
    plt.xlabel('Phi (φ) Angle (degrees)')
    plt.ylabel('Psi (ψ) Angle (degrees)')
    
    plt.axhline(0, color='k', linestyle='--', linewidth=1)
    plt.axvline(0, color='k', linestyle='--', linewidth=1)
    
    plt.xlim(-180, 180)
    plt.ylim(-180, 180)
    plt.grid(True)
    plt.legend()
    
    plt.savefig('ramachandran_plot.png')
    plt.show()

def test_graph(protein_names):
    """
    Test graph to check if plotting works correctly.
    """
    plt.figure(figsize=(8, 6))
    plt.table(cellText=[[name] for name in protein_names], colLabels=['Protein Name'], cellLoc='center', loc='center')
    plt.axis('off')  # Turn off the axis
    plt.title('Protein Names')
    
    plt.savefig('test_graph.png')
    plt.show()

def log_data(pdb_files, phi_angles, psi_angles):
    """
    Log the details of the calculation into a text file.
    """
    with open('ramachandran_plot_log.txt', 'w') as log_file:
        log_file.write("Ramachandran Plot Calculation Log\n")
        log_file.write("Proteins Analyzed:\n")
        for file in pdb_files:
            log_file.write(f"- {os.path.basename(file)}\n")
        log_file.write("\nAngles:\n")
        log_file.write(f"Total Phi Angles: {len(phi_angles)}\n")
        log_file.write(f"Total Psi Angles: {len(psi_angles)}\n")
        log_file.write("\nSample Angles:\n")
        log_file.write(f"First 10 Phi Angles: {phi_angles[:10]}\n")
        log_file.write(f"First 10 Psi Angles: {psi_angles[:10]}\n")

def main():
    # Update this list with the paths to your PDB files
    pdb_files = [
        'C:/Users/evgen/OneDrive/Skrivbord/model/P41968.ent',
        'C:/Users/evgen/OneDrive/Skrivbord/model/2ydv.pdb'
    ]
    
    parser = PDB.PDBParser(QUIET=True)
    all_phi_angles = []
    all_psi_angles = []

    for pdb_file in pdb_files:
        try:
            structure = parser.get_structure('protein', pdb_file)
            
            # Extract phi and psi angles
            phi_angles, psi_angles = extract_phi_psi_angles(structure)
            
            # Append angles to the list
            all_phi_angles.extend(phi_angles)
            all_psi_angles.extend(psi_angles)
        
        except FileNotFoundError:
            print(f"Error: The file '{pdb_file}' was not found. Please check the file path.")
        except Exception as e:
            print(f"An error occurred: {e}")

    if all_phi_angles and all_psi_angles:
        # Plot the Ramachandran plot
        test_proteins = [os.path.basename(pdb_file) for pdb_file in pdb_files]  # Extract file names for the plot subtitle
        plot_ramachandran(all_phi_angles, all_psi_angles, test_proteins)
        
        # Log data
        log_data(pdb_files, all_phi_angles, all_psi_angles)
    else:
        print("No angles were extracted. Please check your PDB files and try again.")
    
    # Test the graph generation
    test_graph(test_proteins)

if __name__ == '__main__':
    main()

