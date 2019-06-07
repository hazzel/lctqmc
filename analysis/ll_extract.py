import glob
import sys
import re
import numpy as np

regex_number = "((?:[-+]?(?:\d*\.\d+|\d+)(?:e[-+]?\d*)?)|(?:[-+]?nan))"
def get_parameter(param, jobfile_string):
	p = re.search(f"\s*{param}\s*=\s*{regex_number}\\n", jobfile_string)
	if p==None:
		return None
	return p.groups()[0]

def get_observable(param, jobfile_string):
	p = re.search(f"\s*Name\s*=\s*{param}\s*Bins\s*=\s*\d*\s*BinLength\s*=\s*\d*\s*Mean\s*=\s*{regex_number}\s*Error\s*=\s*{regex_number}", jobfile_string)
	if p==None:
		return None
	return f"{p.groups()[0]} {p.groups()[1]}"
 
jobname = sys.argv[1]
filenames = glob.glob(jobname + ".task*.out")
filenames.sort()

for filename in filenames:
	with open(filename) as f:
		jobfile_string = f.read()
		output_string = ""
		for i in range(2, len(sys.argv) - 1):
			output_string += f"{get_parameter(sys.argv[i], jobfile_string)} "
		for obs in sys.argv[-1].split():
			output_string += f"{get_observable(obs, jobfile_string)} "
		print(output_string)
