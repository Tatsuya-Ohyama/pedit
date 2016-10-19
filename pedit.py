#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
pedit - PDB editor
"""

import sys, os, re, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse

# =============== common variables =============== #
re_header_atom = re.compile(r"^((ATOM)|(HETATM))")
re_header_ter = re.compile(r"^TER")
re_header_end = re.compile(r"END")


# =============== functions =============== #
# check_file (ファイルの有無を確認)
def check_file(file):
	if not os.path.exists(file):
		sys.stderr.write("ERROR: No such file (%s)\n" % file)
		sys.exit(1)
	return file


# check_overwrite (上書き確認)
def check_overwrite(file):
	if os.path.exists(file):
		sys.stderr.write("WARN: %s exists. Overwrite it? (y/N): " % file)
		sys.stderr.flush()
		user = sys.stdin.readline().replace("\n", "")
		if not user == "y": #or not user == "Y":
			sys.exit(0)


# reset_number (原子順序番号と残基番号のリセット)
def reset_number(file, flag_atom, flag_residue, flag_reset):
	atom_number = 1
	residue_number = 1
	prev_residue_info = ""
	atom_label = ""
	residue_label = ""
	datas = []
	flag_ter = 0

	with open(file, "r") as obj_file:
		for line in obj_file:
			line = line.rstrip("\r\n")
			if re_header_atom.search(line):
				if flag_atom == True:
					# 原子順序番号のリセット
					if atom_number <= 99999:
						# 範囲内の原子順序番号
						atom_label = "%5d" % atom_number
						atom_number += 1
					else:
						# 範囲外の場合
						if flag_reset == True:
							# * をつける
							atom_label = "*****"
						else:
							# 番号をリセット
							atom_number = 1
							atom_label = "%5d" % atom_number
				else:
					atom_label = line[6:11]

				if flag_residue == True:
					residue_info = line[17:26]
					if residue_info != prev_residue_info or flag_ter == 1:
						# 前の残基と残基情報が異なる場合
						if residue_number <= 9999:
							# 範囲内の残基番号
							residue_label = "%4d" % residue_number
							residue_number += 1
						else:
							# 範囲外の場合
							if flag_reset == True:
								# * をつける
								residue_label = "****"
							else:
								# 番号をリセット
								residue_number = 1
								residue_label = "%4d" % residue_number
						prev_residue_info = residue_info
				else:
					residue_label = line[22:26]

				line = line[0:6] + "%5s" % atom_label + line[11:22] + "%4s" % residue_label + line[26:]
				datas.append(line)
				flag_ter = 0

			elif re_header_ter.search(line):
				datas.append("TER")
				flag_ter = 1

			elif re_header_end.search(line):
				datas.append("END")
	return datas

# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = "pedit - PDB editor", formatter_class=argparse.RawTextHelpFormatter)

	subparser = parser.add_subparsers(help = "sub-command")
	subparser.required = True

	parser_view = subparser.add_parser("view", help = "Show PDB in display", formatter_class=argparse.RawTextHelpFormatter)
	parser_view.set_defaults(func = "view")
	parser_view.add_argument("-p", metavar = "PDB", dest = "input", required = True, help = "input (*.pdb)")
	parser_view.add_argument("-o", metavar = "PDB", dest = "output", help = "output (Default: STDIN)")
	parser_view.add_argument("-O", action = "store_true", default = False, help = "overwrite forcibly")

	parser_del = subparser.add_parser("del", help = "Delete atom", formatter_class=argparse.RawTextHelpFormatter)
	parser_del.set_defaults(func = "del")
	parser_del.add_argument("-p", metavar = "PDB", dest = "input", required = True, help = "input (*.pdb)")
	parser_del.add_argument("-o", metavar = "PDB", dest = "output", help = "output (Default: STDIN)")
	parser_del.add_argument("-O", action = "store_true", default = False, help = "overwrite forcibly")

	parser_reset = subparser.add_parser("reset", help = "Delete atom", formatter_class=argparse.RawTextHelpFormatter)
	parser_reset.set_defaults(func = "reset")
	parser_reset.add_argument("-p", metavar = "PDB", dest = "input", required = True, help = "input (*.pdb)")
	parser_reset.add_argument("-R", dest = "flag_residue", action = "store_true", default = False, help = "renumbering residue number")
	parser_reset.add_argument("-A", dest = "flag_atom", action = "store_true", default = False, help = "renumbering atom order")
	parser_reset.add_argument("-S", dest = "flag_reset", action = "store_true", default = False, help = "number set to '*' when number overflow\n(Default: reset number to 1)")
	parser_reset.add_argument("-o", metavar = "PDB", dest = "output", help = "output (Default: STDIN)")
	parser_reset.add_argument("-O", action = "store_true", default = False, help = "overwrite forcibly")

	args = parser.parse_args()


	if args.func == "reset":
		check_file(args.input)
		output = reset_number(args.input, args.flag_atom, args.flag_residue, args.flag_reset)

		if args.output == None:
			# 標準出力に出力
			for data in output:
				print(data)
		else:
			# ファイルに出力
			check_overwrite(args.output)
			with open(args.output, "w") as obj_output:
				for data in output:
					obj_output.write(data + "\n")


