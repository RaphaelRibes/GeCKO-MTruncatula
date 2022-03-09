import math
import os
import sys
import re

class Snp1G:
	def __init__(self):
		self.acgt_counts=[0, 0, 0, 0]
		self.nbReads=0
		self.A1=0
		self.A2=0
		self.p_val = 0

	def get_allele_counts(self, alleles):
		return [str(self.acgt_counts[allele]) for allele in alleles]

class Snp:
	nucl=['A','C','G','T']
	geno_regex = re.compile('\s([A-Z])([A-Z])\|(\S*)')
	alr_regex = re.compile('\s(\d+)\[(\d+)\/(\d+)\/(\d+)\/(\d+)\]')
	def __init__(self, nb_samples1):
		self.nb_samples=nb_samples1
		self.sample_infos=[ Snp1G() for i in range(self.nb_samples)]
		self.alleles=[]
		self.allele_encode=[]
		self.nb_alleles=0
		self.nb_reads=0


	def set_snp(self, line_alr,line_gen):
		allele_counts=[line_gen.count(nucl) for nucl in Snp.nucl]
		self.nb_alleles =len(Snp.nucl) - allele_counts.count(0)
		if self.nb_alleles > 1:
			alleles_by_freq = [allele[0] for allele in sorted(enumerate(allele_counts), key=lambda x: x[1], reverse=True)]
			self.alleles = alleles_by_freq[0:self.nb_alleles]

			gen_infos=Snp.geno_regex.findall(line_gen);
			for i,gen_info in enumerate(gen_infos):
				self.sample_infos[i].A1 = gen_info[0];
				self.sample_infos[i].A2 = gen_info[1];
				self.sample_infos[i].p_val = float(gen_info[2]);

			alr_infos=Snp.alr_regex.findall(line_alr);
			self.nb_reads = 0
			for i,alr_info in enumerate(alr_infos):
				self.sample_infos[i].acgt_counts[0] = int(alr_info[1])
				self.sample_infos[i].acgt_counts[1] = int(alr_info[2])
				self.sample_infos[i].acgt_counts[2] = int(alr_info[3])
				self.sample_infos[i].acgt_counts[3] = int(alr_info[4])
				self.sample_infos[i].nbReads = int(alr_info[0])
				self.nb_reads =  self.nb_reads + int(alr_info[0])

	@staticmethod
	def encode_nucl(nucl):
		if nucl == 'A':
			return 0
		elif nucl == 'C':
			return 1
		elif nucl == 'G':
			return 2
		elif nucl == 'T':
			return 3
		else:
			return -1

	def encode_allele(self, nucl):
		encode_nucl=Snp.encode_nucl(nucl)
		if(encode_nucl==-1):
			return "."
		else:
			return self.alleles.index(Snp.encode_nucl(nucl))

	def to_vcf(self):
		# no SNP id
		# maj and alt allele
		#no QUAL and all kept position pass the 'two alleles' FILTER
		#DP combined depth across samples
		# for each sample we got its genotype GT, read depth per allele AD, and genotype quality GQ
		maj_allele=self.nucl[self.alleles[0]]
		alt_alleles="/".join([self.nucl[allele] for allele in self.alleles[1:]])
		res= ".\t{}\t{}\t.\tPASS\tDP={snp.nb_reads}\tGT:DP:AD:GQ".format(maj_allele, alt_alleles, snp=self)
		snps_info=[res]
		for snp1g in self.sample_infos:
			# the genotype GT A1/A2
			# the total number of reads at this position for this sample DP
			# and the number of reads for majority and alternative alleles AD
			# and the genotype quality GQ, since read2snp only provide the probability for the best genotype we assume the worst case scenario where the second best is 1-p
			# https://gatk.broadinstitute.org/hc/en-us/articles/360035890451-Calculation-of-PL-and-GQ-by-HaplotypeCaller-and-GenotypeGVCFs
			if(snp1g.p_val==1):
				GQ=99;
			elif(snp1g.p_val==0):
				GQ=0
			else:
				GQ= min ( 99, int(round( (-10*math.log10(1.-snp1g.p_val)) - (-10*math.log10(snp1g.p_val)))))
			AD=",".join(snp1g.get_allele_counts(self.alleles))
			snps_info.append("{}/{}:{}:{}:{}".format(self.encode_allele(snp1g.A1),self.encode_allele(snp1g.A2),snp1g.nbReads,AD,GQ))
		return "\t".join(snps_info)



def get_vcf_header():
	res ="##fileformat=VCFv4.0\n"
	res+="##source=reads2snp\n"
	res+="##phasing=unphased\n"
	res+='##FILTER=<ID=multi,Description="At least two alleles are seen at this position">\n'
	res+= '##INFO=<ID=DP,Number=1,Type=Integer,Description="Number of reads at this position (sum across samples)">\n'
	res+= '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
	res+= '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Number of reads at this position for this sample">\n'
	res+= '##FORMAT=<ID=AD,Number=1,Type=Integer,Description="Number of reads at this position for this sample with the reference and alternate variants">\n'
	res+= '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n'

	return res

def alr_to_vcf(fic_alr_name, fic_gen_name,fic_vcf_name ):
	global line_gen
	### STEP 0 : recovery of all individuals.
	nbr_indiv = 0
	vcf_header_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
	with open(fic_alr_name, 'r') as fic_alr:
		line = fic_alr.readline()
		line = fic_alr.readline().strip()
		infos = line.split()
		for sample_name in infos[2:]:
			vcf_header_line += "\t" + sample_name
		nbr_indiv = len(infos) - 2
	# ------------------------------------------------------------------------------------------------------------------------------------
	print "\n\n--------starting VCF conversion (printing one dot every 10K positions) ------------"
	nb_pos=0
	nb_pos_kept=0
	with open(fic_alr_name, 'r') as fic_alr, open(fic_gen_name, 'r') as fic_gen, open(fic_vcf_name, 'w') as fic_vcf:
		fic_vcf.write(get_vcf_header()+vcf_header_line + "\n")
		line_gen = fic_gen.readline().strip()
		line_alr = fic_alr.readline().strip()
		snp = Snp(nbr_indiv)
		while line_gen:
			if line_alr.startswith('>'):
				contig_name = line_alr[1:]
				fic_gen.readline();
				fic_alr.readline();
			else:
				nb_pos=nb_pos+1
				if(nb_pos%10000==0):
					sys.stdout.write(".")
				if line_alr.split()[1] == "P":
					snp.set_snp(line_alr, line_gen)
					if(snp.nb_alleles > 1):
						nb_pos_kept=nb_pos_kept+1
						pos_info="{}\t{}\t".format(contig_name ,line_gen.split()[0])
						fic_vcf.write(pos_info)
						fic_vcf.write(snp.to_vcf())
						fic_vcf.write("\n")
			line_gen = fic_gen.readline().strip()
			line_alr = fic_alr.readline().strip()
		print "\n\n-------VCF conversion completed ------------"
	return nb_pos,nb_pos_kept



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
	if((len(sys.argv) != 4)):
		sys.stderr.write("\n\nInvalid number of parameters\n"
						 "This script converts read2snp output into vcf,"
						 " it expects three parameters:"
						 "\n\ti) the path to the alr file produced by read2snp "
						 "\n\tii) the path to the corresponding gen file produced by read2snp"
						 "\n\tiii) the name of the output vcf file that will be created by this script"
						 "\n\nusage example:"
						 "\npython "+sys.argv[0]+" wheat_capture_2020.alr wheat_capture_2021.gen wheat_capture_2020.vcf\n\n")
		exit(1)

	fic_alr_name = sys.argv[1]
	fic_gen_name = sys.argv[2]
	fic_vcf_name = sys.argv[3]
	nb_pos,nb_kept=alr_to_vcf(fic_alr_name, fic_gen_name, fic_vcf_name)
	sys.stdout.write ("{} positions have been considered and {} variable positions have been kept.\n\n".format(nb_pos,nb_kept))
