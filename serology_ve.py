# Estimate VE using seropositivity for anti-N abs
from pdb import set_trace as brk
from common import *

# For each test the cutoff for the Placebo arm and for the Vax arm
QUANTITATIVE_TESTS = (
	("SARS-CoV-2 serum neutralizing titer 50 (titer) - "
		"Virus Neutralization Assay", 20, 10000),
	("SARS-CoV-2 serum neutralizing titer 90 (titer) - "
		"Virus Neutralization Assay", 20, 3000),
# 	"SARS-CoV-2 serum neutralizing titer 90 to COVID-19 S1 IgG",
# 	("SARS-CoV-2 serum neutralizing titer 50 to COVID-19 S1 IgG", 10, 10),
)

class VECalculator:
	def __init__(self, filename, what, date):
		self.filename = filename
		self.what = what
		self.date = date

		self.data = None
		self.placebo_neg, self.placebo = 0, 0
		self.vax_neg, self.vax = 0, 0
		self._calculated = False

	def _filters(self, before=True):
		return {
				"what": self.what,
				"when": ("V1_DAY1_VAX1_L", "V3_MONTH1_POSTVAX2_L")[not before],
				}

	def load_before(self):
		raise NotImplementedError

	def load_after(self):
		raise NotImplementedError

	def _count_event(self, datum, event):
		if datum.arm == "Placebo":
			self.placebo_neg += int(not event)
			self.placebo += int(event)
		else:
			self.vax_neg += int(not event)
			self.vax += int(event)

	def count(self):
		raise NotImplementedError

	def calculate(self):
		assert not self._calculated

		self.load_before()
		self.load_after()
		self.count()

		self._calculated = True

	def summarize(self):
		print("Looking at {}".format(self.what))

		if self.date:
			print("Considering only rows before {}".format(
				self.date.strftime("%Y-%m-%d")))

		placebo, placebo_neg = self.placebo, self.placebo_neg
		vax, vax_neg = self.vax, self.vax_neg

		vax_total = vax + vax_neg
		placebo_total = placebo + placebo_neg

		print("Vax arm: {} went NEG->POS. {} stayed NEG.".format(vax, vax_neg))
		print("Placebo arm: {} went NEG->POS. {} stayed NEG.".format(placebo,
			placebo_neg))

		print("({} / {}) / ({} / {})".format(vax,
			vax_total, placebo, placebo_total))

		if vax_total == 0 or placebo_total == 0:
			print("Insufficient data")
			return

		try:
			results = ci(vax, vax_neg, placebo, placebo_neg)
			print([to_ve(r) for r in results])
		except ValueError:
			print("Bad data")

class NVECalculator(VECalculator):
	def __init__(self, filename, date):
		super(NVECalculator, self).__init__(filename,
			"N-binding antibody - N-binding Antibody Assay", date)

	def load_before(self):
		filters = self._filters(True)
		filters["result"] = "NEG"

		# Load everyone who was neg to start with
		data = {}
		for datum in load_data(self.filename, None, filters, self.date):
			data[datum.subj_id] = datum

		self.data = data

	def load_after(self):
		# Now update with everyone who was still negative at the end
		filters = self._filters(False)
		filters["result"] = "NEG"
		update_data(self.filename, self.data, filters, self.date)

		# Update again with everyone who was positive at the end
		filters["result"] = "POS"
		update_data(self.filename, self.data, filters, self.date)

	def count(self):
		for datum in self.data.values():
			self._count_event(datum, datum.results == ["NEG", "POS"])

class QuantVECalculator(VECalculator):
	"Estimate VE using one of the quantitative tests"
	def __init__(self, filename, what, placebo_cutoff, vax_cutoff, date):
		self.placebo_cutoff = placebo_cutoff
		self.vax_cutoff = vax_cutoff
		super(QuantVECalculator, self).__init__(filename, what, date)

	def load_before(self):
		filters = self._filters(True)

		# Just load everyone's result on the start date
		data = {}
		for datum in load_data(self.filename, None, filters, self.date):
			data[datum.subj_id] = datum

		self.data = data

	def load_after(self):
		filters = self._filters(False)
		update_data(self.filename, self.data, filters, self.date)

	def count(self):
		for datum in self.data.values():
			results = datum.results
			if len(results) != 2: continue

			cutoff = (self.placebo_cutoff if datum.arm == "Placebo" else
					self.vax_cutoff)

			before, after = [float(r) for r in results]
			event = before <= cutoff and after > cutoff
			self._count_event(datum, event)

def calc_ve(date=None):
	calc = NVECalculator("adva.csv.gz", date)
	calc.calculate()
	calc.summarize()

	calc = QuantVECalculator("adva.csv.gz", *QUANTITATIVE_TESTS[0], date)
	calc.calculate()
	calc.summarize()

def main():
	calc_ve(datetime(2020, 11, 15))
	calc_ve()

if __name__ == "__main__":
	main()
