import csv
import gzip
from collections import namedtuple
from math import *
from datetime import *
from pdb import set_trace as brk

# Datum = namedtuple("Datum", "row_num subj_id arm when what result phase date")

class Datum:
	ATTRS = ("row_num",
			"subj_id",
			"arm",
			"when",
			"what",
			"phase",
			"date",
			"result")

	def __init__(self, *args, **kwargs):
		for attr, arg in zip(self.ATTRS, args):
			setattr(self, attr, arg)

		self.__dict__.update(kwargs)

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
		('DE', 'phase', None),
		('J', 'date', convert_date),
		('W', 'result', None),
		)

	fields = {"row_num": num}

	for index, field, convert in cols:
		d = row[convert_col(index)]
		if convert: d = convert(d)
		fields[field] = d

	return Datum(**fields)

def load_data(filename, filters, date=None):
	"""Generate datums from the rows that match all the filters where date is <
	date"""
	with gzip.open(filename, "rt") as fp:
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
