#!/usr/bin/python

"""
Purpose: Filters rMATS output using cutoffs specified for Pvalue, FDR, IncLevelDiff and ReadCounts
Version: 2.0

Usage: python rmats_filter_lw.py [-h] [-d Directory | -f File] [PValue] [FDR] [IncLevelDiff] [Counts]
Options:
	 Directory	Input directory of rMATs output
	 File		Input file of rMATs output
	 PValue		Keeps events with pvalue lower than specified
         FDR		Keeps events that have FDR lower than specified
         IncLevelDiff	Keeps events with levels greater than specified
	 Counts		Keeps events which has atleast one sample with counts greater than specified

	 To avoid a filtering restriction, use none

Example: python rmats_filter_lw.py -d path/to/directory 0.05 0.1 0.1 20
	 Filtering files for events with p-val < 0.05, fdr < 0.1, psi > 0.1 and counts > 20

	 python rmats_filter_lw.py -f file 0.05 0.2 none none
	 Filtering files for p-val < 0.05, fdr < 0.2, all psi and all counts
"""
from __future__ import division
import sys
import os
import license


if __name__ == '__main__':
	# Check license
	license.check_status()
	if sys.argv[1] == "--help" or sys.argv[1] == "-h":
		print(__doc__)
		sys.exit()

	elif sys.argv[1] == '-l':
		print(license.message())
	elif sys.argv[1] == "-d":
		files_list = os.listdir(sys.argv[2])
		files_list.sort(reverse=True)
		for files in files_list:
			if os.path.isfile(sys.argv[2] + files):
				files = sys.argv[2] + files
				if not license.rmats_filter(files, sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6]):
					print("Could not filter %s, not an rMATS output." % files)
					sys.exit()
		print(" ")
		print("Filtered all files with great success. Your welcome! Dont forget to acknowledge the great leweezard ")
		print("on your publication.")

	elif sys.argv[1] == "-f":
		if not license.rmats_filter(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6]):
			print("Could not filter %s, not an rMATS output." % (sys.argv[2]))
			sys.exit()

		print(" ")
		print("Filtered all files with great success. Your welcome! Dont forget to acknowledge the great leweezard")
		print(" on your publication.")

	else:
		print("Usage: python rmats_filter_lw.py [-h] [-d Directory | -f File] [FDR] [IncLevelDiff] [Counts]")
		print("Refer to python rmats_filter_lw.py -h for documentation.")
