import argparse
import math


def parse_args():
    parser = argparse.ArgumentParser(description="Process PDB file")
    parser.add_argument("input_file", help="input PDB file")
    parser.add_argument("output_file", help="output PDB file")
    parser.add_argument("target_chain_id", help="chain ID of target residue")
    parser.add_argument("residue_id", type=int, help="residue ID of target residue")
    parser.add_argument("target_label", help="label for target atoms")
    parser.add_argument("cutoff_distance", type=float, help="cutoff distance in angstroms")
    return parser.parse_args()


def process_pdb(input_file, output_file, target_chain_id, residue_id, target_label, cutoff_distance):
    atoms = []
    target_atoms = []
    target_residue_found = False
    hr_residues = {}
    lr_residues = {}
    hr_atoms = {}
    with open(input_file, "r") as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                record_name = line[0:6].strip()
                atom_serial_number = int(line[6:11])
                atom_name = line[12:16].strip()
                alt_location_indicator = line[16]
                residue_name = line[16:20].strip()
                chain_id = line[21]
                residue_number = int(line[22:26])
                insertion_code = line[26]
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                occupancy = float(line[54:60])
                b_factor = float(line[60:66])
                element_symbol = line[76:78].strip()
                charge = line[78:80].strip()

                if residue_number == residue_id and chain_id == target_chain_id:
                    target_residue_found = True
                    target_atoms.append((record_name, atom_name, residue_name, chain_id, residue_number, x, y, z, occupancy, b_factor, args.target_label))
                else:
                    atoms.append((record_name, atom_name, residue_name, chain_id, residue_number, x, y, z, occupancy, b_factor))

    if not target_residue_found:
        print("Target residue not found in input file.")
        return

    fa_atoms = set()
    for atom in target_atoms:
        for other_atom in atoms:
            dx = atom[5] - other_atom[5]
            dy = atom[6] - other_atom[6]
            dz = atom[7] - other_atom[7]
            distance = math.sqrt(dx*dx + dy*dy + dz*dz)
            if distance <= cutoff_distance:
                fa_atoms.add(other_atom)

    for i, atom in enumerate(atoms):
        if atom in fa_atoms:
            atoms[i] = (atom[0], atom[1], atom[2], atom[3], atom[4], atom[5], atom[6], atom[7], atom[8], atom[9], "FA")
        else:
            atoms[i] = (atom[0], atom[1], atom[2], atom[3], atom[4], atom[5], atom[6], atom[7], atom[8], atom[9], "CG")

    fa_residues = set((atom[3], atom[4]) for atom in fa_atoms)

    for target_atom in target_atoms:
        atoms.append((target_atom[0], target_atom[1], target_atom[2], target_atom[3], target_atom[4], target_atom[5],
                      target_atom[6], target_atom[7], target_atom[8], target_atom[9], target_atom[10]))
    atoms.sort(key=lambda x: (x[3], x[4]))

    with open(output_file, "w") as f:
        prev_residue = None
        for i, atom in enumerate(atoms):
            current_residue = (atom[3], atom[4])
            if current_residue != prev_residue:
                if current_residue in fa_residues:
                    res_label = "FA"
                    hr_atoms.setdefault(atom[3], 0)
                    hr_atoms[atom[3]] += 1
                    hr_residues.setdefault(atom[3], set())
                    hr_residues[atom[3]].add(atom[4])
                elif current_residue == (args.target_chain_id, args.residue_id):
                    res_label = args.target_label
                    hr_atoms.setdefault(atom[3], 0)
                    hr_atoms[atom[3]] += 1
                    hr_residues.setdefault(atom[3], set())
                    hr_residues[atom[3]].add(atom[4])
                else:
                    res_label = "CG"
                    lr_residues.setdefault(atom[3], set())
                    lr_residues[atom[3]].add(atom[4])
            else:
                res_label = "FA"
                hr_atoms.setdefault(atom[3], 0)
                hr_atoms[atom[3]] += 1
                hr_residues.setdefault(atom[3], set())
                hr_residues[atom[3]].add(atom[4])
            line = "{:6s}{:>6d}  {:<5s}{:s}{:3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}  {:7.2f}    {:2s}       \n".format(
                atom[0], i + 1, atom[1], "", atom[2], atom[3], atom[4], atom[5], atom[6], atom[7], atom[9], res_label)

            f.write(line)

    total_nodes = 0
    with open("info.out", "w") as f:
        for chain_id in sorted(set([atom[3] for atom in atoms])):
            if cutoff_distance > 0:
                hr_atoms_chain = hr_atoms.get(chain_id, {})
                lr_residues_chain = lr_residues.get(chain_id, set())
                if hr_atoms_chain:
                    node_c = len(lr_residues_chain) + hr_atoms_chain
                    total_nodes += node_c
                    print(total_nodes)
                else:
                    node_c = len(lr_residues_chain)
                    total_nodes += node_c
                    print(total_nodes)
            elif cutoff_distance == 0.0:
                hr_atoms_chain = hr_atoms.get(chain_id, {})
                lr_residues_chain = lr_residues.get(chain_id, set())
                node_c = len(lr_residues_chain)
                total_nodes += node_c
                print(total_nodes)

            f.write("Chain {}: FA residues={}, CG residues={}, FA atoms={}, Nodes={}\n".format(
                    chain_id, len(hr_residues.get(chain_id, set())), len(lr_residues.get(chain_id, set())), hr_atoms.get(chain_id, {}), node_c))
            f.write("\nFA residues={}\nCG residues={}\n\n".format(hr_residues.get(chain_id, set()), lr_residues_chain))
            f.write("\nTotal number of Nodes={}".format(total_nodes))



if __name__ == "__main__":
    args = parse_args()
    print(args.input_file, args.output_file, args.target_chain_id, args.residue_id,  args.target_label, args.cutoff_distance,)
    process_pdb(args.input_file, args.output_file, args.target_chain_id, args.residue_id, args.target_label, args.cutoff_distance )



