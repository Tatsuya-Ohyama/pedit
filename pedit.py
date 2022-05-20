#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
pedit - PDB editor
"""

import sys, signal
sys.dont_write_bytecode = True
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
import parmed

from mods.func_prompt_io import check_exist, check_overwrite



# =============== constant =============== #
RESIDUE_WAT = ["WAT", "HOH", "SOL"]
RESIDUE_AMINO = [
	"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
	"LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
	"HIE", "HID", "HIP", "CYM", "ACE", "NME",
]
RESIDUE_NUC = [
	"DA5", "DA", "DA3",
	"DC5", "DC", "DC3",
	"DG5", "DG", "DG3",
	"DT5", "DT", "DT3",
	"A5", "A", "A3",
	"C5", "C", "C3",
	"G5", "G", "G3",
	"U5", "U", "U3",
]



# =============== functions =============== #
def view(obj_mol):
	prev_chain = None
	n_atom = 0
	n_residue = 0
	n_wat = 0
	n_sol = 0
	n_amino = 0
	n_nuc = 0
	n_other = 0
	list_other = set([])
	for obj_residue in obj_mol.residues:
		if prev_chain != obj_residue.chain:
			if prev_chain is None:
				print("*** Chain {0} ***".format(obj_residue.chain))
			else:
				print("\n*** Chain {0} ***".format(obj_residue.chain))
			print("{0:^8} {1:^5}".format("Residue", "nAtom"))
			print("{0:-^3} {0:-^4} {0:-^5}".format(""))
			prev_chain = obj_residue.chain
		print("{0:3} {1:>4} {2:>5}".format(obj_residue.name, obj_residue.number, len(obj_residue.atoms)))
		n_residue += 1
		n_atom += len(obj_residue.atoms)
		if obj_residue.name in RESIDUE_WAT:
			n_wat += 1
			continue

		n_sol += 1
		if obj_residue.name in RESIDUE_AMINO:
			n_amino += 1
		elif obj_residue.name in RESIDUE_NUC:
			n_nuc += 1
		else:
			n_other += 1
			list_other.add(obj_residue.name)

	print("")
	print("Atoms:        {0:>8}".format(n_atom))
	print("Residues:     {0:>8}".format(n_residue))
	print(" > Solute:    {0:>8}".format(n_sol))
	print("  >> Amino:   {0:>8}".format(n_amino))
	print("  >> Nucleic: {0:>8}".format(n_nuc))
	if len(list_other) == 0:
		print("  >> Other:   {0:>8}".format(n_other))
	else:
		print("  >> Other:   {0:>8} ({1})".format(n_other, ",".join(list_other)))
	print(" > Water:     {0:>8}".format(n_wat))


def edit(obj_mol, mask_del):
	if mask_del is not None:
		obj_mol.strip(mask_del)
	return obj_mol



# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "pedit - PDB editor", formatter_class=argparse.RawTextHelpFormatter)

	subparser = parser.add_subparsers()
	subparser.required = True

	common_argument_group = argparse.ArgumentParser(add_help=False)
	common_argument_group.add_argument("-i", dest="INPUT_FILE", metavar="INPUT.pdb", required=True, help="input (*.pdb)")

	parser_view = subparser.add_parser("view", help="Show PDB in display", parents=[common_argument_group])
	parser_view.set_defaults(handler="view")

	parser_edit = subparser.add_parser("edit", help="Edit structure", parents=[common_argument_group])
	parser_edit.set_defaults(handler="edit")
	parser_edit.add_argument("-o", dest="OUTPUT_FILE", metavar="OUTPUT.pdb", required=True, help="output")
	parser_edit.add_argument("-d", dest="MASK_DEL", metavar="AMBERMASK_DEL", help="Ambermask for delete atoms (`:{residue numlist}`, `@{atom numlist}`, `:{residue namelist}`, `@{atom namelist}`, `@%%{atom type name}`, `@/{atom_element_name}`)")
	parser_edit.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite forcibly")
	parser_edit.add_argument("--no-renumber", dest="FLAG_NO_RENUMBER", action="store_true", default=False, help="do not renumber of atoms and residues")

	args = parser.parse_args()


	check_exist(args.INPUT_FILE, 2)

	obj_mol = parmed.load_file(args.INPUT_FILE)

	# 処理
	if args.handler == "view":
		view(obj_mol)

	elif args.handler == "edit":
		obj_mol = edit(obj_mol, args.MASK_DEL)

		if args.FLAG_OVERWRITE == False:
			check_overwrite(args.OUTPUT_FILE)

		if args.FLAG_NO_RENUMBER:
			obj_mol.save(args.OUTPUT_FILE, overwrite=True, renumber=False)
		else:
			obj_mol.save(args.OUTPUT_FILE, overwrite=True)

	else:
		sys.stderr.write("ERROR: undefined sub-command.\n")



