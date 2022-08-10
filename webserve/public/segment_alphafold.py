#!/usr/bin/env python3
"""
Jiffy to go through an mmCIF file or a PDB file for Alphafold entries and split into segments,
based on some algorithm. The B-factor column needs to contain the confidence value from 
AlphaFold or you will get very unexpected results
"""
#
# normal imports
#
from pathlib import Path
import gzip
import socket
from time import sleep
from argparse import ArgumentParser

#
# special imports
#
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


#
# functions
#
def run_dssp(pdb_file):
    """
    get the secondary structure analysis from DSSP and store the
    elements in an array
    """
    import subprocess

    if "bc.ic.ac.uk" in socket.getfqdn():
        # Jurassic DSSP...
        command = f"/bmm/soft/Linux_2.6_64/bin/dssp -na {pdb_file} temp.dssp"
        subprocess.check_output(list(command.split()), text=True)
        with open("temp.dssp", "r") as filein:
            structure = filein.read()
        Path("temp.dssp").unlink()
    else:
        command = f"mkdssp -i {pdb_file}"
        structure = subprocess.check_output(list(command.split()), text=True)
    started = False
    elements = {}
    #
    # this is to remind me of what they are, not actually used anywhere 17.05.2022
    #
    """
    decode = {
        "H": "alpha helix",
        "B": "beta bridge",
        "E": "extended strand in beta ladder",
        "G": "310 helix",
        "I": "5 helix (pi-helix)",
        "T": "hydrogen bonded turn",
        "S": "bend",
        " ": " "
    }
    """
    #
    # some secondary structure elements are more "important" in domains
    #
    factor = {
        "H": 2.5,
        "B": 2.5,
        "E": 2.5,
        "G": 1.0,
        "I": 1.0,
        "T": 1.0,
        "S": 1.0,
        " ": 0.75,
    }
    for line in structure.split("\n"):
        if line.startswith("  #  RESIDUE AA STRUCTURE"):
            started = True
            continue
        if started:
            if len(line) == 0:
                break
            try:
                residue = int(line.split()[1])
            except Exception:
                #
                # "!" in [1] if His-His - this is an extra line, just go to next entry
                #
                continue
            aa = line[14:15]
            ss_element = line[16:17]
            elements[residue] = [aa, ss_element, factor[ss_element]]
            # print(residue, ss_element, decode[ss_element])
        else:
            continue
    return elements


#
# initialisations
#
parser = MMCIFParser()
three2one = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLY": "G",
    "GLN": "Q",
    "GLU": "E",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TYR": "Y",
    "TRP": "W",
    "VAL": "V",
}
#
# copy to one2three...
#
one2three = {}
for key, item in three2one.items():
    one2three[item] = key
#
# functions
#


def product(this_list):
    """
    jiffy to calculate product of list elements
    """
    if len(this_list) == 0:
        return 1.0
    result = 1.0
    for item in this_list:
        result *= item
    return result


class ProteinModel:
    """
    populates a model with atom names, coordinates, etc.
    builds sequence of this bit of the model
    """

    def __init__(self, model={}, sequence=""):
        """
        take the various bits of the "_atom_site" group of an mmCIF
        and make a model suitable for writing to a short-form PDB file,
        and extract useful things like sequences
        """
        self.pdb_entry = Path()
        self.model = {}
        self.sequence = ""

    def make_model_from_mmcif(self, mmcif_dict):
        """
        make a model object from the protein model in an MMCIF file
        """
        for key in mmcif_dict["_atom_site.id"]:
            ikey = int(key) - 1
            try:
                if mmcif_dict["_atom_site.group_PDB"][ikey] == "ATOM":
                    self.model[key] = {
                        "label": mmcif_dict["_atom_site.group_PDB"][ikey],
                        "atom_number": int(mmcif_dict["_atom_site.id"][ikey]),
                        "residue_number": int(
                            mmcif_dict["_atom_site.auth_seq_id"][ikey]
                        ),
                        "residue_type": mmcif_dict["_atom_site.auth_comp_id"][ikey],
                        "residue_type_1": three2one[
                            mmcif_dict["_atom_site.auth_comp_id"][ikey]
                        ],
                        "atom_name": mmcif_dict["_atom_site.auth_atom_id"][ikey],
                        "atom_type": mmcif_dict["_atom_site.type_symbol"][ikey],
                        "chain": mmcif_dict["_atom_site.auth_asym_id"][ikey],
                        "coord_x": float(mmcif_dict["_atom_site.Cartn_x"][ikey]),
                        "coord_y": float(mmcif_dict["_atom_site.Cartn_y"][ikey]),
                        "coord_z": float(mmcif_dict["_atom_site.Cartn_z"][ikey]),
                        "occupancy": float(mmcif_dict["_atom_site.occupancy"][ikey]),
                        "B_factor": float(
                            mmcif_dict["_atom_site.B_iso_or_equiv"][ikey]
                        ),
                    }
            except IndexError as error:
                print(key, ikey)
                print(mmcif_dict["_atom_site.id"])
                print(f"Error line is {error}")
                exit()
        minimum_confidence = 999
        maximum_confidence = -1
        for key, entry in self.model.items():
            if entry["B_factor"] >= maximum_confidence:
                maximum_confidence = entry["B_factor"]
            if entry["B_factor"] <= minimum_confidence:
                minimum_confidence = entry["B_factor"]
        print(
            f"minimum confidence = {minimum_confidence}, maximum = {maximum_confidence}"
        )

    def make_model_from_pdb(self, protein):
        """
        make a model object from the protein model in a PDB file
        """
        for line in protein:
            if line.startswith("ATOM"):
                #
                # column spacing is very, very important in a PDB file!
                #           1         2         3         4         5         6         7         8
                # 0123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 #
                # ATOM      1  N   MET A   1     -24.555 -26.246   8.302  1.00 36.79            N
                # aaaaaaNNNNN aaaa aaa cNNNN    nnnn.nnnNNNN.NNNnnnn.nnn  N.NNnnn.nn           AA
                key = int(line[6:11].strip())
                self.model[key] = {
                    "label": line[0:6].strip(),
                    "atom_number": int(line[6:11]),
                    "atom_name": line[12:16].strip(),
                    "residue_type": line[17:20],
                    "residue_type_1": three2one[line[17:20].strip()],
                    "chain": line[21:22],
                    "residue_number": int(line[22:26].strip()),
                    "coord_x": float(line[30:38].strip()),
                    "coord_y": float(line[38:46].strip()),
                    "coord_z": float(line[46:54].strip()),
                    "occupancy": float(line[56:60].strip()),
                    "B_factor": float(line[61:66].strip()),
                    "atom_type": line[77:78],
                }

    def make_sequence(self):
        """
        make a single-letter amino acid sequence from the C-alphas in a protein model
        """
        for key in self.model:
            if self.model[key]["atom_name"] == "CA":
                self.sequence += three2one[self.model[key]["residue_type"]]

    def make_calpha(self):
        """
        extract a C-alpha trace from a protein model object
        """
        calpha_trace = {}
        for key in self.model:
            if self.model[key]["atom_name"] == "CA":
                calpha_trace[key] = self.model[key].copy()
        return calpha_trace

    def make_copy(self):
        """
        make a copy of an existing protein model object
        """
        model_copy = ProteinModel()
        for key in self.model:
            model_copy.model[key] = self.model[key].copy()
        model_copy.sequence = self.sequence
        return model_copy

    def populate_model(self, mmcif_dict):
        """
        run methods to populate parts of a protein model object
        """
        self.make_model_from_mmcif(mmcif_dict)
        self.make_sequence()

    def write_pdb_file(self, filename):
        """
        method to write PDB file from model object
        """
        pdb_out = Path(filename)
        line = f"HEADER file for {pdb_out}"
        padding = 80 - len(line)
        line += (" " * padding) + "\n"
        with open(pdb_out, "w") as fileout:
            fileout.write(line)
            for key in self.model:
                atom = self.model[key]
                #
                # column spacing is very, very important in a PDB file!
                #           1         2         3         4         5         6         7         8
                # 0123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 #
                # ATOM      1  N   MET A   1     -24.555 -26.246   8.302  1.00 36.79            N
                line = f"ATOM  {atom['atom_number']:5d}  "
                line += f"{atom['atom_name']:<4}"
                line += f"{atom['residue_type']} "
                line += f"{atom['chain']}"
                line += f"{atom['residue_number']:4d}    "
                line += f"{atom['coord_x']:8.3f}"
                line += f"{atom['coord_y']:8.3f}"
                line += f"{atom['coord_z']:8.3f}  "
                line += f"{atom['occupancy']:4.2f}"
                line += f"{atom['B_factor']:6.2f}           "
                line += f"{atom['atom_type']}\n"
                fileout.write(line)
            line = f"TER{' '*77}\nEND{' '*77}\n"
            fileout.write(line)


#
# end of class ProteinModel
#


def prune_domain(model, threshold=50):
    # self.model = {key: atom for (key, atom) in self.model.items() if self.model[key]["B_factor"] >= threshold}
    prune = True
    nu_model = {}
    reversed_model = []
    # print("before pruning, domain has", len(model), "atoms")
    for key, atom in model.items():
        reversed_model.append(key)
        if prune:
            if model[key]["B_factor"] >= threshold:
                prune = False
                nu_model[key] = atom
        else:
            nu_model[key] = atom
    prune = True
    for key in reversed_model[::-1]:
        if prune:
            if nu_model[key]["B_factor"] >= threshold:
                prune = False
                nu_model[key] = model[key]
            else:
                del nu_model[key]
    model = nu_model.copy()


def write_pdb_chain(model, pdb_stem):
    """
    standalone function to write PDB file from an atomic model
    """
    counting = False
    for key in model:
        if not counting:
            counting = True
            first_residue = model[key]["residue_number"]
    last_residue = model[key]["residue_number"]
    # for key in model:
    #    domain = model[key]["chain"]
    # print(f"writing domain {domain} {first}-{last} to {pdb_stem}_{domain}.pdb")
    pdb_out = Path(f"{pdb_stem}_{first_residue}_{last_residue}.pdb")
    print(f"writing pdb file {pdb_out}")
    # if pdb_out.exists():
    #    return
    # print("writing", pdb_out, end="")
    if len(pdb_stem) > 4:
        line = f"HEADER Domain from AlphaFold model for UniProt entry {pdb_stem}"
    else:
        line = f"HEADER Domain from genuine PDB structure {pdb_stem}"
    padding = 80 - len(line)
    line += (" " * padding) + "\n"
    start = ""
    with open(pdb_out, "w") as fileout:
        fileout.write(line)
        for key in model:
            atom = model[key]
            if start == "":
                start = atom["residue_number"]
            line = f"ATOM  {atom['atom_number']:5d}  "
            line += f"{atom['atom_name']:<4}"
            line += f"{atom['residue_type']} "
            line += f"{atom['chain']}"
            line += f"{atom['residue_number']:4d}    "
            line += f"{atom['coord_x']:8.3f}"
            line += f"{atom['coord_y']:8.3f}"
            line += f"{atom['coord_z']:8.3f}  "
            line += f"{atom['occupancy']:4.2f}"
            line += f"{atom['B_factor']:6.2f}           "
            line += f"{atom['atom_type']}\n"
            fileout.write(line)
        line = f"TER{' '*77}\nEND{' '*77}\n"
        fileout.write(line)
        debug(f" {start}-{atom['residue_number']}")
    return


def write_pdb_sequence(model, pdb_stem):
    """
    standalone function to write FASTA file from an atomic model
    """
    for key in model:
        #
        # only need first hit
        #
        domain = model[key]["chain"]
        break
    fasta_out = Path(f"{pdb_stem}_{domain}.fasta")
    header = f">{pdb_stem}_{domain}\n"
    sequence = ""
    for key in model:
        if model[key]["atom_name"] == "CA":
            sequence += model[key]["residue_type_1"]
    sequence += "\n"
    # print(f"writing sequence to {fasta_out.name}")
    with open(fasta_out, "w") as fileout:
        fileout.write(header)
        fileout.write(sequence)
    return


def debug(argument, debug_value=False):
    if debug_value:
        print(argument)
    return


"""
main program starts here
"""


def main(arguments):
    # unsuitable = []
    model_file = Path(arguments.input)
    min_domain_length = int(arguments.minimum_domain_length)
    list_length = int(arguments.loop_length)
    # domain_break_value = pow(0.4, list_length)
    min_val = float(arguments.confidence)
    counter = 0
    model_count = 0
    debug("=" * 64)
    debug("templates (for Phyre2) were extracted")
    counter += 1
    # print(counter,":",model_file)
    #
    # get a file stem name for the output
    #
    if model_file.name.startswith("AF-"):
        pdb_stem = model_file.name.split("AF-")[1].split("-F1-model")[0]
    else:
        pdb_stem = model_file.name[:4]
    #
    # have we already written PDB output files for this model? Do it anyway if this is a single entry
    #
    if str(model_file).endswith(".cif"):
        # protein = parser.get_structure("alpha_model", model_file)
        # Don't do multi-model structures
        mmcif_dict = MMCIF2Dict(model_file)
        model = ProteinModel()
        model.populate_model(mmcif_dict)
    elif str(model_file).endswith(".pdb") or (
        str(model_file).startswith("pdb") and str(model_file).endswith(".ent")
    ):
        model = ProteinModel()
        with open(model_file, "r") as filein:
            lines = filein.readlines()
            # model = ProteinModel()
            model.make_model_from_pdb(lines)
            model.make_sequence()
    elif str(model_file).endswith(".pdb.gz"):
        model = ProteinModel()
        with gzip.open(model_file, "r") as filein:
            lines = []
            count = 0
            for line in filein:
                nuline = line.decode()
                if nuline.startswith("TITLE") and "COLLAGEN" in nuline:
                    debug(nuline)
                    continue
                lines.append(nuline)
            model = ProteinModel()
            model.make_model_from_pdb(lines)
            model.make_sequence()
    else:
        print("Don't know the format - file must end with .cif, .pdb, .pdb.gz or .ent")
        print("filename is", str(model_file))
        exit()
    #
    # get the secondary structure elements
    #
    model.write_pdb_file("temp.pdb")
    ss_elements = run_dssp("temp.pdb")
    # Path("temp.pdb").unlink()
    #
    # get a file stem name for the output
    #
    if model_file.name.startswith("AF-"):
        pdb_stem = model_file.name.split("AF-")[1].split("-F1-model")[0]
    else:
        pdb_stem = model_file.name[:4]
    #
    # variables for splitting a chain into domains
    #
    count = 0
    domain_dict = {}
    sequence = {}
    domain_list = list(chr(i) for i in range(65, 91))  # A-Z
    # add lowercase
    domain_list[len(domain_list) :] = [chr(i) for i in range(97, 122)]  # a-z
    # add numerals
    domain_list[len(domain_list) :] = [chr(i) for i in range(48, 57)]  # 0-9
    loop_list = []
    value_list = []
    #
    # try simpler criterion - simple sum
    #
    # this_list = [1.0] * list_length
    this_list = [100.0] * list_length
    domain_break_value = min_val * list_length
    this_domain = 0
    domain_break = False
    started = False
    #
    # divide into domains based on some criteria
    #
    calpha_model = model.make_calpha()
    # calc mean val of confidence values
    meanval = 0
    median_list = []
    #
    # step through this model to initialise some values
    #
    for key in calpha_model:
        atom = calpha_model[key]
        meanval += atom["B_factor"]
        median_list.append(atom["B_factor"])
    # print(f"MEAN VAL = {meanval/divisor:.2f}, MEDIAN = {statistics.median(median_list):.2f}")
    min_val = max(min_val, 0.6 * meanval / len(median_list))

    domain_break_value = min_val * list_length
    # print(domain_break_value, min_val, 0.6 * meanval / len(median_list), list_length, meanval, len(median_list))
    #
    # step through this model
    #
    for key, atom in calpha_model.items():
        #
        # adjust confidence (Aka "B_factor") if it's in a secondary structure elephant
        #
        try:
            atom["B_factor"] *= ss_elements[atom["residue_number"]][2]
        except KeyError as error:
            print(error)
            print(f"key error for {pdb_stem}, {atom['residue_number']}")
            print(atom)
            for line in ss_elements:
                print(line)
            exit()
        if atom["B_factor"] >= min_val:
            started = True
        this_list.pop(0)
        this_list.append(atom["B_factor"])
        value = sum(this_list)
        value_list.append(value)
        if not started:
            if atom["B_factor"] < min_val:
                # add this to a loop, i.e. not a domain
                loop_list.append(atom["residue_number"])
                domain_break = True
            else:
                # start a new domain
                started = True
        elif value < domain_break_value:
            if not domain_break:
                domain_break = True
                this_domain += 1
                loop_list = this_list.copy()
            else:
                loop_list.append(atom["residue_number"])
        else:
            if domain_break:
                domain_break = False
                if len(loop_list) > 0 and len(loop_list) < list_length:
                    if this_domain > 0:
                        this_domain -= 1
                    for item in loop_list:
                        if domain_list[this_domain] in domain_dict:
                            domain_dict[domain_list[this_domain]].append(item)
                        else:
                            domain_dict[domain_list[this_domain]] = [item]
                loop_list = []
            try:
                if domain_list[this_domain] in domain_dict:
                    domain_dict[domain_list[this_domain]].append(atom["residue_number"])
                else:
                    domain_dict[domain_list[this_domain]] = [atom["residue_number"]]
            except IndexError as error:
                print(f"\nfor {model_file}:")
                print(f"current index is {this_domain}")
                print(f'Error message is "{error}"')
                try:
                    print(f"current index is {domain_list[this_domain]}")
                except IndexError as error:
                    print("domain_list is", len(domain_list), "long")
                    print(f'Error message is "{error}"')
                    print(domain_dict)
                exit()
            if domain_list[this_domain] in sequence:
                sequence[domain_list[this_domain]] += atom["residue_type_1"]
            else:
                sequence[domain_list[this_domain]] = atom["residue_type_1"]
        count += 1
    #
    # domain_dict now contains a list of residue numbers for each domain
    #
    #
    # prune the low confidence ends of each domain
    #
    threshold = 50
    threshold = min_val
    for chain in domain_dict:
        #
        # head
        #
        atom_list = []
        prune = True
        for key, atom in calpha_model.items():
            atom_list.append(key)
            if atom["residue_number"] in domain_dict[chain]:
                if atom["B_factor"] >= threshold:
                    prune = False
                if prune:
                    domain_dict[chain].remove(atom["residue_number"])
        #
        # tail
        #
        prune = True
        for key in atom_list[::-1]:
            if calpha_model[key]["residue_number"] in domain_dict[chain]:
                if calpha_model[key]["B_factor"] >= threshold:
                    prune = False
                if prune:
                    domain_dict[chain].remove(calpha_model[key]["residue_number"])
    #
    # now go through the domains and remove any that are shorter than min_domain_length
    #
    pop_list = []
    for key in domain_dict:
        if len(domain_dict[key]) < min_domain_length:
            pop_list.append(key)
            # print(
            #     f"will delete domain {key} because it has length {len(domain_dict[key])}",
            #     f"({domain_dict[key][0]}-{domain_dict[key][-1]})",
            # )
    if len(pop_list) > 0:
        nu_domain_dict = {}
        nu_sequence = {}
        count = 0
        for key in domain_dict:
            if key not in pop_list:
                nu_domain_dict[domain_list[count]] = domain_dict[key]
                nu_sequence[domain_list[count]] = sequence[key]
                count += 1
        domain_dict = nu_domain_dict.copy()
        sequence = nu_sequence.copy()
    #
    # go through the domains and see if any have gaps < list_length; if so, stitch them
    # together again - should this be done before removing short ones or after?
    #
    last = 0
    pop_list = []
    for key in domain_dict:
        # print("key = ", key, "start and end = ", domain_dict[key][0], domain_dict[key][-1])
        if last != 0:
            if domain_dict[key][0] < domain_dict[last][-1] + list_length:
                # print("------- stitching domains", last, "and", key, "together --------")
                domain_dict[last] = [
                    x for x in range(domain_dict[last][0], domain_dict[key][-1] + 1)
                ]
                pop_list.append(key)
        last = key
    if len(pop_list) > 0:
        nu_domain_dict = {}
        nu_sequence = {}
        count = 0
        for key in domain_dict:
            if key not in pop_list:
                nu_domain_dict[domain_list[count]] = domain_dict[key]
                nu_sequence[domain_list[count]] = sequence[key]
                count += 1
        domain_dict = nu_domain_dict.copy()
        sequence = nu_sequence.copy()
    #
    # split model into chains
    #
    dict_domain = {}
    for key in domain_dict:
        for item in domain_dict[key]:
            dict_domain[item] = key
    chain = {}
    for key in model.model:
        atom = model.model[key].copy()
        if atom["residue_number"] in dict_domain:
            this_chain = dict_domain[atom["residue_number"]]
            atom["chain"] = this_chain
            if this_chain not in chain:
                chain[this_chain] = {}
            chain[this_chain][key] = atom
    for key in chain:
        model_count += 1
        write_pdb_chain(chain[key], pdb_stem)
        # write_pdb_sequence(chain[key], pdb_stem)


if __name__ == "__main__":
    #
    # command-line parser...
    #
    parser = ArgumentParser()
    help_line = "name of mmcif or pdb format file from AlphaFold (the "
    help_line += (
        '"B-factor" column must contain the confidence values and not real B-factors)'
    )
    parser.add_argument(
        "input", help=help_line,
    )
    parser.add_argument(
        "-l",
        "--loop_length",
        default=12,
        required=False,
        help="minimum number of consecutive low confidence residues to use when calculating a break",
    )
    parser.add_argument(
        "-m",
        "--minimum_domain_length",
        default=30,
        required=False,
        help='minimum number of residues in a "domain"',
    )
    parser.add_argument(
        "-c",
        "--confidence",
        default=50,
        required=False,
        help="minimum value of confidence for pruning out residues",
    )
    config = parser.parse_args()
    main(config)
    exit()
