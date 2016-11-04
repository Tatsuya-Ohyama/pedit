#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
pedit - PDB editor
"""

import sys, os, re, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
import tempfile
import shutil

# =============== common variables =============== #
re_header_atom = re.compile(r"^((ATOM)|(HETATM))")
re_header_ter = re.compile(r"^TER")
re_header_end = re.compile(r"END")
re_header_connect = re.compile(r"^CONECT")
re_wsp = re.compile(r"^[\s\t]*$")
re_notatom = re.compile(r"[\d'\s]")
re_range = re.compile(r"^((\d+)|(\d+-\d+)|(\d+-)|(-\d+))$")


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


# target_range (範囲指定の解析)
def range_analyze(data, min_value, max_value):
	names = []
	names_condition = []
	numbers = []
	numbers_conditon = []

	datas = data.split(":", 2)
	if re_range.search(datas[1]):
		# 番号による範囲指定
		datas = datas[1].split("-", 2)
		if len(datas) == 1:
			# 単体の場合
			numbers.append([int(datas[0]), int(datas[0])])
		else:
			# 範囲の場合
			if datas[0] == "":
				# 始端が指定されていない場合
				datas[0] = min_value
			if datas[1] == "":
				# 終端が指定されていない場合
				datas[1] = max_value
			datas = list(map(lambda x : int(x), datas))
			if datas[1] < datas[0]:
				# 始端が終端より大きい場合
				sys.stderr.write("ERROR: Invalid target range (%s). End value is small than start.\n" % data)
				sys.exit(1)
			if not (min_value <= datas[0] <= datas[1] <= max_value):
				sys.stderr.write("ERROR: Invalid target range (%s). Range size overflow or underflow from in PDB format\n" % data)
				sys.exit(1)
			numbers.append(datas)

	else:
		# 名前による指定
		names.append(datas[1])

	return names, numbers


# def check_line(line, term1, term2 = []):
# 	if term1[0] == "A":
# 		# 原子
# 		if type(term1[1]) == int:
# 			# 範囲
# 		else:
# 			# 名前
# 	elif


# class ChekTerm:
# 	def __init__(self):
# 		self.targets = []

# 	def getTarget(self, targets):
# 		# 条件取り込み
# 		for target in targets:
# 			queries = target.split("@", 2)
# 			category = 0
# 			if re_target_atom.search(queries[1]):
# 				category = 1
# 				(tmp_names, tmp_numbers) = range_analyze(queries[0], 1, 99999)
# 				targets_atom_name.extend(tmp_names)
# 				targets_atom_num.extend(tmp_numbers)
# 				if len(queries) == 2:
# 					# 条件がある場合
# 					category =
# 					(tmp_names, tmp_numbers) = range_analyze(queries[1], 1, 99999)
# 					targets_atom_name_cnd.extend(tmp_names)
# 					targets_atom_num_cnd.extend(tmp_numbers)
# 			elif re_target_residue.search(queries[1]):
# 				category = 3
# 				(tmp_names, tmp_numbers) = range_analyze(queries[0], 1, 9999)
# 				targets_residue_name.extend(tmp_names)
# 				targets_residue_num.extend(tmp_numbers)
# 				if len(queries) == 2:
# 					# 条件がある場合
# 					category = 4
# 					(tmp_names, tmp_numbers) = range_analyze(queries[1], 1, 99999)
# 					targets_residue_name_cnd.extend(tmp_names)
# 					targets_residue_num_cnd.extend(tmp_numbers)
# 			else:
# 				sys.stderr.write("ERROR: target type is not specified (%s)\n" % target)
# 				sys.exit(1)




# 			self.condition.append([category, term1, term2])

	# def checkTerm(self, line):
	# 	# 該当行を削除するか判断
	# 	atom_name = re_notatom.sub("", line[12:14].strip())
	# 	atom_number = int(line[6:11])
	# 	residue_name = line[17:20]
	# 	residue_number = int(line[22:26])

	# 	flag_del = 0
	# 	for term in self.targets:
	# 		if term[0] == "A" and len(term) == 2:
	# 			# 第一条件は原子
	# 			if isinstance(term[1], list):
	# 				# リストの場合、範囲による指定
	# 				if len(list(filter(lambda array : array[0] <= int(line[6:11]) <= array[1], term[]))) != 0:

	# 				break
	# 			else:
	# 				# 名前による指定
	# 				break

	# 		elif term[0] == "A" and len(term) == 4:
	# 			if isinstance(term[1], list):
	# 				# リストの場合、範囲による指定
	# 				break
	# 			else:
	# 				# 名前による指定
	# 				break


	# 		elif term[0] == "R" and len(term) == 2:
	# 			if isinstance(term[1], list):
	# 				# リストの場合、範囲による指定
	# 				break
	# 			else:
	# 				# 名前による指定
	# 				break



	# 		elif term[0] == "R" and len(term) == 4:
	# 			if isinstance(term[1], list):
	# 				# リストの場合、範囲による指定
	# 				break
	# 			else:
	# 				# 名前による指定
	# 				break


	# 	if flag_del == 0:
	# 		return True
	# 	else:
	# 		return False


	# def printTarget(self):
	# 	# 条件を出力
	# 	return self.condition


# del_object (対象の削除)
def del_object(file, targets, flag_discone):
	re_target_atom = re.compile(r"^A:", re.IGNORECASE)
	re_target_residue = re.compile(r"^R:", re.IGNORECASE)
	targets_atom_name = []
	targets_atom_name_cnd = []
	targets_atom_num = []
	targets_atom_num_cnd = []
	targets_residue_name = []
	targets_residue_name_cnd = []
	targets_residue_num = []
	targets_residue_num_cnd = []
	for target in targets:
		queries = target.split("@", 2)
		print(target, queries)
		if re_target_atom.search(target):
			(tmp_names, tmp_numbers) = range_analyze(queries[0], 1, 99999)
			targets_atom_name.extend(tmp_names)
			targets_atom_num.extend(tmp_numbers)
			if len(queries) == 2:
				# 条件がある場合
				(tmp_names, tmp_numbers) = range_analyze(queries[1], 1, 99999)
				targets_atom_name_cnd.extend(tmp_names)
				targets_atom_num_cnd.extend(tmp_numbers)
		elif re_target_residue.search(target):
			(tmp_names, tmp_numbers) = range_analyze(queries[0], 1, 9999)
			targets_residue_name.extend(tmp_names)
			targets_residue_num.extend(tmp_numbers)
			if len(queries) == 2:
				# 条件がある場合
				(tmp_names, tmp_numbers) = range_analyze(queries[1], 1, 99999)
				targets_residue_name_cnd.extend(tmp_names)
				targets_residue_num_cnd.extend(tmp_numbers)
		else:
			sys.stderr.write("ERROR: target type is not specified (%s)\n" % target)
			sys.exit(1)


	print("targets_atom_name", targets_atom_name)
	print("targets_atom_name_cnd", targets_atom_name_cnd)
	print("targets_atom_num", targets_atom_num)
	print("targets_atom_num_cnd", targets_atom_num_cnd)
	print("targets_residue_name", targets_residue_name)
	print("targets_residue_name_cnd", targets_residue_name_cnd)
	print("targets_residue_num", targets_residue_num)
	print("targets_residue_num_cnd", targets_residue_num_cnd)


	flag_ter = 0
	datas = []
	with open(file, "r") as obj_file:
		for line in obj_file:
			line = line.rstrip("\r\n")
			if re_header_atom.search(line):
				if len(targets_atom_num) != 0 and not re_wsp.search(line[6:11]):
					# 原子順序番号で指定
					if len(list(filter(lambda array : array[0] <= int(line[6:11]) <= array[1], targets_atom_num))) != 0:
						continue

				if len(targets_residue_num) != 0 and not re_wsp.search(line[22:26]):
					# 残基番号で指定
					if len(list(filter(lambda array : array[0] <= int(line[22:26]) <= array[1], targets_residue_num))) != 0:
						continue

				if len(targets_atom_name) != 0:
					# 原子名で指定
					if re_notatom.sub("", line[12:14]) in targets_atom_name:
						continue

				if len(targets_residue_name) != 0:
					# 残基名で指定
					if line[17:20].strip() in targets_residue_name:
						continue
				flag_ter = 0
				datas.append(line)

			elif re_header_ter.search(line):
				if(flag_ter == 1):
					# 前に TER を扱っている場合 (TER の連続対策)
					continue
				datas.append("TER")
				flag_ter = 1

			elif re_header_end.search(line):
				datas.append("END")

			elif re_header_connect.search(line):
				if not flag_discone == True:
					datas.append(line)

	return datas



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
	parser_view.add_argument("input", metavar = "PDB", help = "input (*.pdb)")
	parser_view.add_argument("-o", metavar = "PDB", dest = "output", help = "output (Default: STDIN)")
	parser_view.add_argument("-O", dest = "overwrite", action = "store_true", default = False, help = "overwrite forcibly")

	parser_del = subparser.add_parser("del", help = "Delete atom", formatter_class=argparse.RawTextHelpFormatter)
	parser_del.set_defaults(func = "del")
	parser_del.add_argument("input", metavar = "PDB", help = "input (*.pdb)")
	parser_del.add_argument("-t", dest = "targets", nargs = "+", help = "target objects\n  Ex:\n    R:1-20 (Residues 1-20)\n    R:WAT (Residue name of 'WAT')\n    A:-20 (Atoms 1-20)\n    A:H (Hydrogen atoms)")
	parser_del.add_argument("-C", action = "store_true", dest = "flag_discone", help = "delete connections")
	parser_del.add_argument("-o", metavar = "PDB", dest = "output", help = "output (Default: STDIN)")
	parser_del.add_argument("-O", dest = "overwrite", action = "store_true", default = False, help = "overwrite forcibly")

	parser_reset = subparser.add_parser("reset", help = "Delete atom", formatter_class=argparse.RawTextHelpFormatter)
	parser_reset.set_defaults(func = "reset")
	parser_reset.add_argument("input", metavar = "PDB", help = "input (*.pdb)")
	parser_reset.add_argument("-R", dest = "flag_residue", action = "store_true", default = False, help = "renumbering residue number")
	parser_reset.add_argument("-A", dest = "flag_atom", action = "store_true", default = False, help = "renumbering atom order")
	parser_reset.add_argument("-S", dest = "flag_reset", action = "store_true", default = False, help = "number set to '*' when number overflow\n(Default: reset number to 1)")
	parser_reset.add_argument("-o", metavar = "PDB", dest = "output", help = "output (Default: STDIN)")
	parser_reset.add_argument("-O", dest = "overwrite", action = "store_true", default = False, help = "overwrite forcibly")

	args = parser.parse_args()


	check_file(args.input)

	# 処理
	if args.func == "reset":
		output = reset_number(args.input, args.flag_atom, args.flag_residue, args.flag_reset)

	elif args.func == "del":
		output = del_object(args.input, args.targets, args.flag_discone)


	# 出力
	if args.output == None:
		# 標準出力に出力
		for data in output:
			print(data)
	else:
		# ファイルに出力
		if args.overwrite == False:
			check_overwrite(args.output)
		tempfile_name = ""
		with tempfile.NamedTemporaryFile(mode = "w", prefix = "", delete = False) as obj_output:
			tempfile_name = obj_output.name
			for data in output:
				obj_output.write(data + "\n")

		shutil.move(tempfile_name, args.output)

