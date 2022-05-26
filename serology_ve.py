# Estimate VE using seropositivity for anti-N abs
from pdb import set_trace as brk
from common import *

def calc_ve(date=None):
	filters = {
			"what": "N-binding antibody - N-binding Antibody Assay",
			"when": "V1_DAY1_VAX1_L",
			"result": "NEG",
			}

	# Find the people who were negative to start with, and the totals in each
	# arm
	placebo_arm, vax_arm = set(), set()

	for datum in load_data("adva.csv.gz", filters, date):
		if datum.arm == "Placebo":
			placebo_arm.add(datum.subj_id)
		else:
			vax_arm.add(datum.subj_id)

	# Now see how many in each arm became positive afterwards
	placebo, vax = 0, 0
	filters.update({"when": "V3_MONTH1_POSTVAX2_L", "result": "POS"})
	for datum in load_data("adva.csv.gz", filters, date):

		if datum.subj_id in placebo_arm:
			placebo += 1

		elif datum.subj_id in vax_arm:
			vax += 1

	# And how many were still negative (don't just assume anyone we don't have
	# a result for at the second visit was negative)
	placebo_neg, vax_neg = 0, 0
	filters.update({"when": "V3_MONTH1_POSTVAX2_L", "result": "NEG"})
	for datum in load_data("adva.csv.gz", filters, date):

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
