#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
pedit - PDB editor
"""

import sys, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
import tempfile
import shutil
import parmed

from basic_func import check_exist, check_overwrite

# =============== variables =============== #


# =============== functions =============== #


# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "pedit - PDB editor", formatter_class=argparse.RawTextHelpFormatter)

	subparser = parser.add_subparsers(help = "sub-command")
	subparser.required = True

	parser_view = subparser.add_parser("view", help = "Show PDB in display", formatter_class=argparse.RawTextHelpFormatter)
	parser_view.set_defaults(func = "view")
	parser_view.add_argument("input_file", dest = "input_file", metavar = "INPUT.pdb", required = True, help = "input (*.pdb)")

	parser_mod = subparser.add_parser("del", help = "Delete atom", formatter_class=argparse.RawTextHelpFormatter)
	parser_mod.set_defaults(func = "mod")
	parser_mod.add_argument("-i", dest = "input_file", metavar = "INPUT.pdb", required = True, help = "input file")
	parser_mod.add_argument("-o", dest = "output_file", metavar = "OUTPUT.pdb", required = True, help = "output file")
	parser_mod.add_argument("-O", dest = "overwrite", action = "store_true", default = False, help = "overwrite forcibly")


	parser_del.add_argument("-t", dest = "targets", nargs = "+", help = "target objects\n  Ex:\n    R:1-20 (Residues 1-20)\n    R:WAT (Residue name of 'WAT')\n    A:-20 (Atoms 1-20)\n    A:H (Hydrogen atoms)")
	parser_del.add_argument("-C", action = "store_true", dest = "flag_discone", help = "delete connections")


	args = parser.parse_args()


	check_file(args.input_file)

	# 処理
	obj_mol = parmed.load_file(args.input_file)
	pdb = PDBFile(args.input_file)
