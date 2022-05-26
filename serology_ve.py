# Estimate VE using seropositivity for anti-N abs
import csv
import gzip
from collections import namedtuple
from math import *
from datetime import *
from pdb import set_trace as brk

Datum = namedtuple("Datum", "row_num subj_id arm when what result phase date")

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

def convert_date(s):
	try:
		return datetime.strptime(s, "%Y-%m-%d")
	except ValueError:
		return datetime(2020, 1, 1)

def parse_row(row, num):
	cols = (
		('D', 'subj_id', None),
		('H', 'arm', None),
		('O', 'when', None),
		('S', 'what', None),
		('W', 'result', None),
		('DE', 'phase', None),
		('J', 'date', convert_date),
		)
	fields = {"row_num": num}
	for index, field, convert in cols:
		d = row[convert_col(index)]
		if convert: d = convert(d)
		fields[field] = d
	return Datum(**fields)

def load_data(filters, date=None):
	"""Generate datums from the rows that match all the filters where date is <
	date"""
	with gzip.open("adva.csv.gz", "rt") as fp:
		r = csv.reader(fp)
		for i, row in enumerate(r):
			datum = parse_row(row, i+1)

			if date and datum.date >= date:
				continue

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

def calc_ve(date=None):
	filters = {
			"what": "N-binding antibody - N-binding Antibody Assay",
			"when": "V1_DAY1_VAX1_L",
			"result": "NEG",
			}

	# Find the people who were negative to start with, and the totals in each
	# arm
	placebo_arm, vax_arm = set(), set()

	for datum in load_data(filters, date):
		if datum.arm == "Placebo":
			placebo_arm.add(datum.subj_id)
		else:
			vax_arm.add(datum.subj_id)

	# Now see how many in each arm became positive afterwards
	placebo, vax = 0, 0
	filters.update({"when": "V3_MONTH1_POSTVAX2_L", "result": "POS"})
	for datum in load_data(filters, date):

		if datum.subj_id in placebo_arm:
			placebo += 1

		elif datum.subj_id in vax_arm:
			vax += 1

	# And how many were still negative (don't just assume anyone we don't have
	# a result for at the second visit was negative)
	placebo_neg, vax_neg = 0, 0
	filters.update({"when": "V3_MONTH1_POSTVAX2_L", "result": "NEG"})
	for datum in load_data(filters, date):

		if datum.subj_id in placebo_arm:
			placebo_neg += 1

		elif datum.subj_id in vax_arm:
			vax_neg += 1

	if date:
		print("Considering only rows before {}".format(
			date.strftime("%Y-%m-%d")))

	vax_total = vax + vax_neg
	placebo_total = placebo + placebo_neg

	print("Vax arm: {} went NEG->POS. {} stayed NEG.".format(vax, vax_neg))
	print("Placebo arm: {} went NEG->POS. {} stayed NEG.".format(placebo,
		placebo_neg))

	print("({} / {}) / ({} / {})".format(vax,
		vax_total, placebo, placebo_total))

	results = ci(vax, vax_neg, placebo, placebo_neg)
	print([to_ve(r) for r in results])

def main():
	calc_ve()
	calc_ve(datetime(2020, 11, 15))

if __name__ == "__main__":
	main()
