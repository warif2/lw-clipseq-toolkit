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
import csv
import sys
import os

def list_avg(list):
	sum = 0
	for elm in list:
		sum += int(elm)
	avg = sum / len(list)
	return avg

def rmats_filter(infile, pvalue_cut, fdr_cut, diff_cut, count_cut):

	summary_out = csv.writer(open("all_events_filtered_Pval-lt-%s_FDR-lt-%s_Diff-gt-%s_Count-gt-%s.csv" % (pvalue_cut, fdr_cut, diff_cut, count_cut), "a"), delimiter=",")

	with open(infile, 'r') as f:
		header = f.readline().strip("\n").split("\t")
		if "FDR" and "IncLevelDifference" not in header:
			print "Missing FDR or IncLevelDiffence in %s" % (infile)
			return False

		if "PValue" not in header:
			print "Missing Pvalue in %s" % (infile)
			return False

		if "IC_SAMPLE_1" not in header:
			if "IJC_SAMPLE_1" not in header:
				print "Missing either IC_SAMPLE_1 or IJC_SAMPLE_2 in %s" % (infile)
				return False

		if "IC_SAMPLE_2" not in header:
			if "IJC_SAMPLE_2" not in header:
				print "Missing either IC_SAMPLE_2 or IJC_SAMPLE_2 in %s" % (infile)
				return False

		if "SC_SAMPLE_1" not in header:
			if "SJC_SAMPLE_1" not in header:
				print "Missing either SC_SAMPLE_1 or SJC_SAMPLE_1 in %s" % (infile)
				return False

		if "SC_SAMPLE_2" not in header:
			if "SJC_SAMPLE_2" not in header:
				print "Missing either SC_SAMPLE_2 or SJC_SAMPLE_2 in %s" % (infile)
				return False

	outfile = infile.split("/")[-1].strip('.txt')
	event = outfile.split(".")[0]
	type = [event]
	output = csv.writer(open("%s_filtered_PVal-lt-%s_FDR-lt-%s_Diff-gt-%s_Count-gt-%s.csv" % (outfile, pvalue_cut, fdr_cut, diff_cut, count_cut), "w"), delimiter="\t")

	if pvalue_cut == "none":
		pvalue_cut = 2
	if fdr_cut == "none":
		fdr_cut = 2
	if diff_cut == "none":
		diff_cut = 0

	with open(infile, 'r') as f:
		header = f.readline().strip("\n").split("\t")

		if event == "SE" or event == "A3SS" or event == "A5SS" or event == "RI":
			output.writerow(header[0:2]+["Event"]+header[2:11]+[""]+[""]+header[11:23])
		elif event == "MXE":
			output.writerow(header[0:2]+["Event"]+header[2:25])
		else:
			output.writerow(header)

		if event == "SE":
			summary_out.writerow(header[0:2]+["Event"]+header[2:11]+[""]+[""]+header[11:23])

		pvalue_col = header.index("PValue")
		fdr_col = header.index("FDR")
		diff_col = header.index("IncLevelDifference")
		if "IC_SAMPLE_1" in header:
			ic_1_col = header.index("IC_SAMPLE_1")
			sc_1_col = header.index("SC_SAMPLE_1")
			ic_2_col = header.index("IC_SAMPLE_2")
			sc_2_col = header.index("SC_SAMPLE_2")

		if "IJC_SAMPLE_1" in header:
			ic_1_col = header.index("IJC_SAMPLE_1")
			sc_1_col = header.index("SJC_SAMPLE_1")
			ic_2_col = header.index("IJC_SAMPLE_2")
			sc_2_col = header.index("SJC_SAMPLE_2")

		for line in csv.reader(f, delimiter='\t'):

			if float(line[pvalue_col]) < float(pvalue_cut) and float(line[fdr_col]) < float(fdr_cut) and abs(float(line[diff_col])) > float(diff_cut):
				copy = 0
				if count_cut == "none":
					if event == "SE" or event == "RI" or event == "A5SS" or event == "A3SS":
						output.writerow(line[0:2]+type+line[2:11]+[""]+[""]+line[11:23])
						summary_out.writerow(line[0:2]+type+line[2:11]+[""]+[""]+line[11:23])
					elif event == "MXE":
						output.writerow(line[0:2]+type+line[2:23])
						summary_out.writerow(line[0:2]+type+line[2:25])
					else:
						output.writerow(line)
						summary_out.writerow(line)
				else:
					inc_count_1 = list_avg(line[ic_1_col].split(","))
					skp_count_1 = list_avg(line[sc_1_col].split(","))
					inc_count_2 = list_avg(line[ic_2_col].split(","))
					skp_count_2 = list_avg(line[sc_2_col].split(","))
					if inc_count_1 < float(count_cut) and skp_count_1 < float(count_cut) and inc_count_2 < float(count_cut) and skp_count_2 < float(count_cut):
						copy = 1

					if copy == 0:
						if event == "SE" or event == "RI" or event == "A5SS" or event == "A3SS":
							output.writerow(line[0:2]+type+line[2:11]+[""]+[""]+line[11:23])
							summary_out.writerow(line[0:2]+type+line[2:11]+[""]+[""]+line[11:23])
						elif event == "MXE":
							output.writerow(line[0:2]+type+line[2:23])
							summary_out.writerow(line[0:2]+type+line[2:25])
						else:
							output.writerow(line)
							summary_out.writerow(line)
	print "Finished with %s" % (infile)
	return True




if __name__ == '__main__':

	print " "
	if sys.argv[1] == "--help" or sys.argv[1] == "-h":
		print __doc__
		sys.exit()

	elif sys.argv[1] == "-d":

		files_list = os.listdir(sys.argv[2])
		files_list.sort(reverse=True)
		for files in files_list:
			if os.path.isfile(sys.argv[2] + files):
				files = sys.argv[2] + files
				if not rmats_filter(files, sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6]):
					print "Could not filter %s, not an rMATS output." % (files)

		print " "
		print "Filtered all files with great success. Your welcome! Dont forget to acknowledge the great leweezard on your publication."

	elif sys.argv[1] == "-f":
		if not rmats_filter(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6]):
			print "Could not filter %s, not an rMATS output." % (sys.argv[2])

		print " "
		print "Filtered all files with great success. Your welcome! Dont forget to acknowledge the great leweezard on your publication."
	else:
		print "Usage: python rmats_filter_lw.py [-h] [-d Directory | -f File] [FDR] [IncLevelDiff] [Counts]"
		print "Refer to python rmats_filter_lw.py -h for documentation."