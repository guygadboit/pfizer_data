# Estimate VE using seropositivity for anti-N abs
import csv
import gzip
from collections import namedtuple
from math import *
from pdb import set_trace as brk

Datum = namedtuple("Datum", "row_num subj_id arm when what result phase")

def convert_col(col):
	"convert a column index like AM back to decimal"
	ret = 0
	mul = 1
	offset = 0

	for d in reversed(col):
		ret += mul * (ord(d) - ord('A') + offset)
		mul *= 26
		offset = 1

	return ret

def parse_row(row, num):
	cols = (
		('D', 'subj_id'),
		('F', 'arm'),
		('O', 'when'),
		('S', 'what'),
		('W', 'result'),
		('DE', 'phase'),
		)
	fields = {"row_num": num}
	for index, field in cols:
		fields[field] = row[convert_col(index)]
	return Datum(**fields)

def load_data(filters):
	"Generate datums from the rows that match all the filters"
	with gzip.open("adva.csv.gz", "rt") as fp:
		r = csv.reader(fp)
		for i, row in enumerate(r):
			datum = parse_row(row, i+1)

			for attr, val in filters.items():
				if getattr(datum, attr) != val:
					break
			else:
				yield datum

def ci(a, b, c, d):
	"""a, b: events and non-events in one arm, c,d: the same in the other arm.
	Return the 95CI and the point estimate of the RR"""
	a = float(a)
	c = float(c)

	def rr():
		return (a/(a+b)) / (c/(c+d))

	r = rr()
	f = sqrt((1/a) + (1/c) - (1/(a+b)) - (1/(c+d)))
	e = exp(1.96*f)

	return r/e, r, r*e

def to_ve(rr):
	return (1-rr) * 100

def calc_ve():
	filters = {
			"what": "N-binding antibody - N-binding Antibody Assay",
			"when": "V1_DAY1_VAX1_L",
			"result": "NEG",
			"phase": "Phase 3",
			}

	# Find the people who were negative to start with, and the totals in each
	# arm
	group = set()
	placebo_total, vax_total = 0, 0

	for datum in load_data(filters):
		if datum.subj_id in group: continue
		group.add(datum.subj_id)
		if datum.arm == "Placebo":
			placebo_total += 1
		else:
			vax_total += 1

	print("Found {} interesting people".format(len(group)))

	# Now see how many in each arm became positive afterwards
	placebo, vax = 0, 0
	filters.update({"when": "V3_MONTH1_POSTVAX2_L", "result": "POS"})
	for datum in load_data(filters):

		if datum.subj_id not in group:
			continue

		if datum.arm == "Placebo":
			placebo += 1
		else:
			vax += 1

	print("{} / {} in vax arm went from NEG to POS".format(vax, vax_total))
	print("{} / {} in placebo arm went from NEG to POS".format(placebo,
		placebo_total))

	print("({} / {}) / ({} / {})".format(vax,
		vax_total, placebo, placebo_total))

	results = ci(vax, vax_total - vax, placebo, placebo_total - placebo)
	print([to_ve(r) for r in results])

def main():
	calc_ve()

if __name__ == "__main__":
	main()
