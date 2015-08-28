#!/bin/env python3
import subprocess
import sys
import os
import pandas as pd

VCFPATH="~/data/NDAR/Daly/VCFs"
METAPATH="~/data/NDAR/Daly/metafiles"

def read_omics_exp():
	''' The second line is skipped as it is explaining the header row'''
	fn="omics_experiments.txt"
	omics_exp = pd.read_csv(os.path.join(METAPATH, fn),quotechar='"',delimiter="\t", header=0, skiprows=1)
	return omics_exp

def read_genomics_subj():
	''' The second line is skipped as it is explaining the header row '''
	fn = "genomics_subject02.txt"
	genomics_subj = pd.read_csv(os.path.join(METAPATH, fn),quotechar='"',delimiter="\t", header=0, skiprows=[1])#, index_col=31)
	return genomics_subj

def read_genomics_sample():
	''' The second line is skipped as it is explaining the header row '''
	fn = "genomics_sample03.txt"
	genomics_samp = pd.read_csv(os.path.join(METAPATH, fn),quotechar='"',delimiter="\t", header=0, skiprows=[1])
	return genomics_samp

def vcf_to_IDs(mypick="broad"):
	if mypick == "broad":
		fn="c1_NIMH_Autism_HGSC_Broad_broad.vcf"
	elif mypick == "baylor":
		fn = "c1_NIMH_Autism_HGSC_BCM_FINAL_hg19.vcf"
	else:
		print("Error: must be broad or baylor")
		sys.exit(1)
	inds=subprocess.check_output( 'grep -m 1 "#CHROM" %s' % os.path.join(VCFPATH, fn), shell=True )
	inds=inds.decode('utf-8').rstrip().split("\t")[9:]
	return inds

def get_ped_for_vcf(mypick="broad"):
	genomics_subj = read_genomics_subj()
	genomics_samp = read_genomics_sample()
	inds = vcf_to_IDs(mypick)
	genomics_subj_for_inds = genomics_subj[genomics_subj['sample_id_biorepository'].isin(inds)]
	genomics_samp_for_inds = genomics_samp[genomics_samp['sample_id_biorepository'].isin(inds)]
	myfields = ['family_user_def_id','sample_id_biorepository', 'subjectkey_father','subjectkey_mother','gender','phenotype','race','interview_age']
	myfields2 = ['experiment_id']
	mytable = genomics_subj_for_inds.ix[:,myfields]
	mytable.to_csv(sys.stdout, index=False, sep="\t")
	## output experimental ID
	mytable2 = genomics_samp_for_inds.ix[:,myfields2]
	mytable2.to_csv(sys.stderr, index=False, sep="\t")

def QC(mypick="broad"):
	genomics_subj = read_genomics_subj()
	inds = vcf_to_IDs(mypick)
	inds_with_meta = set(genomics_subj[genomics_subj['sample_id_biorepository'].isin(inds)]['sample_id_biorepository'])
	inds_all = set(inds)
	inds_wo_meta = inds_all - inds_with_meta
	print(inds_wo_meta)
	print("Totally %d" % len(inds_wo_meta))
	genomics_samp = read_genomics_sample()
	inds_with_samp = set(genomics_samp[genomics_samp['sample_id_biorepository'].isin(inds)]['sample_id_biorepository'])
	inds_wo_samp = inds_all - inds_with_samp
	print(inds_wo_samp)
	print("Totally %d" % len(inds_wo_samp))


def test():
	omics_exp=read_omics_exp()
	print(omics_exp["experiment_id"])
	genomics_subj = read_genomics_subj()
	print(genomics_subj['sample_id_biorepository'])
	broad_ids = vcf_to_IDs()
	print(' '.join(broad_ids[0:10]) + "...")
	

def main():
	#test()
	#QC('broad')
	get_ped_for_vcf('broad')

if __name__ == "__main__":
	main()	
