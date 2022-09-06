# -*- coding: utf-8 -*-
# @author: Wenzhi Gao
# @date: 2022/09/03
# @email: gwz@163.shufe.edu.cn

from typing import List
from random import normalvariate as randn

def gen_csv_file(b: List, fname: str) -> None:

	if (len(b) != 35):
		raise RuntimeError("Invalid size of b")

	b_str = ",".join([str(b[i]) for i in range(len(b))])

	with open(fname, "w") as f:
		f.write(b_str)

	return

if __name__ == '__main__':

	n_samples = 100

	b = [-1.20771696000000,-0.324534720000000,-0.814544800000000,-0.470015680000000,-1.12810432000000,
		  0.651398400000000,1.61020640000000,1.48922560000000,0.191568320000000,-0.170299840000000,
		  0.395972480000000,0.234927040000000,-2.05527616000000,-0.596973120000000,-2.58581840000000,
		  0.674370240000000,-0.600968640000000,0.821288960000000,-1.26244832000000,0.117156800000000,
		  -0.315636159999999,-0.491426560000000,-0.405572160000000,0.796003520000000,-0.139378240000000,
		  0.161203040000000,0.621456960000000,-1.16646112000000,-0.744688320000000,1.02231744000000,
		  -1.06181024000000,-1.00564352000000,0.812587520000000,-0.718058560000000,-2.11628192000000]

	gen_csv_file(b, "rot_example.csv")

	for i in range(n_samples):
		b_sample = [b[i] + randn(0.0, 0.01) for i in range(len(b))]
		gen_csv_file(b_sample, "rot_example_rand_{0}.csv".format(i + 1))