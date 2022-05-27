# Estimate VE using seropositivity for anti-N abs
from pdb import set_trace as brk
from common import *

def calc_ve(date=None):
	filters = {
			"what": "N-binding antibody - N-binding Antibody Assay",
			"when": "V1_DAY1_VAX1_L",
			"result": "NEG",
			}

	# Load everyone who was neg to start with
	data = {}
	for datum in load_data("adva.csv.gz", None, filters, date):
		data[datum.subj_id] = datum

	# Now update with everyone who was still negative at the end
	filters.update({"when": "V3_MONTH1_POSTVAX2_L", "result": "NEG"})
	update_data("adva.csv.gz", data, filters, date)

	# Update again with everyone who was positive at the end
	filters.update({"when": "V3_MONTH1_POSTVAX2_L", "result": "POS"})
	update_data("adva.csv.gz", data, filters, date)

	# Now count NEG->NEG and NEG->POS in each arm
	placebo_neg, placebo = 0, 0
	vax_neg, vax = 0, 0

	for datum in data.values():
		non_event = datum.results == ["NEG", "NEG"]
		event = datum.results == ["NEG", "POS"]

		if datum.arm == "Placebo":
			placebo_neg += int(non_event)
			placebo += int(event)
		else:
			vax_neg += int(non_event)
			vax += int(event)

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
