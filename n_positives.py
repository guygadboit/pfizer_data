from common import load_data
from pdb import set_trace as brk

def main():
	filters = {
			"what": "N-binding antibody - N-binding Antibody Assay",
			"when": "V3_MONTH1_POSTVAX2_L",
			"result": "POS"
			}

	people = set()

	for datum in load_data("adva.csv.gz", filters):
		people.add(datum.subj_id)

	print("{} people who were anti-N positive".format(len(people)))

	for test in (
			"COVID-19 RBD IgG (U/mL) - Luminex Immunoassay",
			"COVID-19 S1 IgG (U/mL) - Luminex Immunoassay",
			"SARS-CoV-2 serum neutralizing titer 50 (titer) - Virus Neutralization Assay",
			"SARS-CoV-2 serum neutralizing titer 50 to COVID-19 S1 IgG",
			"SARS-CoV-2 serum neutralizing titer 90 (titer) - Virus Neutralization Assay",
			"SARS-CoV-2 serum neutralizing titer 90 to COVID-19 S1 IgG",
			):

		filters = {
				"what": test,
				"when": "V3_MONTH1_POSTVAX2_L",
				}

		print("Looking at {} in anti-N positive".format(test))

		count = 0
		for datum in load_data("adva.csv.gz", filters):
			if datum.subj_id in people:
				print(datum.arm, datum.result)
				count += 1

		print("Found {} (out of {})".format(count, len(people)))

if __name__ == "__main__":
	main()
