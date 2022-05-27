# Estimate VE using seropositivity for anti-N abs
from pdb import set_trace as brk
from common import *

QUANTITATIVE_TESTS = (
	"COVID-19 RBD IgG (U/mL) - Luminex Immunoassay",
	"COVID-19 S1 IgG (U/mL) - Luminex Immunoassay",
	"SARS-CoV-2 serum neutralizing titer 50 (titer) - Virus Neutralization Assay",
	"SARS-CoV-2 serum neutralizing titer 50 to COVID-19 S1 IgG",
	"SARS-CoV-2 serum neutralizing titer 90 (titer) - Virus Neutralization Assay",
	"SARS-CoV-2 serum neutralizing titer 90 to COVID-19 S1 IgG",
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

		results = ci(vax, vax_neg, placebo, placebo_neg)
		print([to_ve(r) for r in results])

class NVECalculator(VECalculator):
	def __init__(self, filename, date):
		super(NVECalculator, self).__init__(filename,
			"N-binding antibody - N-binding Antibody Assay", date)

	def load_before(self):
		filters = self._filters(True)
		filters["result"] = "NEG"

		# Load everyone who was neg to start with
		data = {}
		for datum in load_data("adva.csv.gz", None, filters, self.date):
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
		# Don't count twice

		for datum in self.data.values():
			non_event = datum.results == ["NEG", "NEG"]
			event = datum.results == ["NEG", "POS"]

			if datum.arm == "Placebo":
				self.placebo_neg += int(non_event)
				self.placebo += int(event)
			else:
				self.vax_neg += int(non_event)
				self.vax += int(event)

def calc_ve(date=None):
	calc = NVECalculator("adva.csv.gz", date)
	calc.calculate()
	calc.summarize()

def main():
	calc_ve()
	calc_ve(datetime(2020, 11, 15))

if __name__ == "__main__":
	main()
