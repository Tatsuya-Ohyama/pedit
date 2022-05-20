# pedit.py

## 概要
PDB ファイルを簡易的に編集するプログラム




## 使用方法
### View モード
構造の概要を表示するモード

```sh
$ pedit.py view [-h] -i INPUT.pdb
```

* `-h`, `--help`
	: このヘルプメッセージを表示して終了する。
* `-i INPUT.pdb`
	: 入力 PDB ファイル


### Edit モード
構造を編集するモード

```sh
$ pedit.py edit [-h] -i INPUT.pdb -o OUTPUT.pdb [-d AMBERMASK_DEL] [-O] [--no-renumber]
```

* `-h`, `--help`
	: このヘルプメッセージを表示して終了する。
* `-i INPUT.pdb`
	: 入力 PDB ファイル
* `-o OUTPUT.pdb`
	: 出力 PDB ファイル
* `-d AMBERMASK_DEL`
	: 削除する原子の Ambermask

	* `:{residue numlist}`
	* `@{atom numlist}`
	* `:{residue namelist}`
	* `@{atom namelist}`
	* `@%{atom type name}`
	* `@/{atom_element_name}`

* `-O`
	: 上書きプロンプトを出さずに上書きする。
* `--no-renumber`
	: 原子および残基インデックスを振り直さない


## 動作要件
* Python3
	* parmed


## License
The MIT License (MIT)

Copyright (c) 2016 Tatsuya Ohyama


## Authors
* Tatsuya Ohyama


## ChangeLog
### Ver. 1.0 (2022-05-20)
* 統合版を追加した。
* 従来版を廃止した。
