__author__ = 'gidis'

import errno
import os
import pickle
import re
import subprocess
import sys
import glob
from Bio import SeqIO

import snapperdb

class ParsedVcf:
    def __init__(self):
        self.depth = {}
        self.qual = {}
        self.hap_qual = {}
        self.hap_depth = {}
        self.hap_call = {}
        self.hap_var_count = {}
        self.filter_flag = {}
        self.var = {}
        self.ref_base = {}
        self.ref = None
        self.mixed_positions = []
        self.bad_pos = []
        self.good_var = []

# -------------------------------------------------------------------------------------------------


class Vcf:
    def __init__(self):
        self.bad_depth = None
        self.bad_qual = None
        self.bad_var = None
        self.good_var = None
        self.bad_pos = None
        self.sample_name = None
        self.depth_average = 0.0
        self.depth_sd = None
        self.path_to_config = None
        self.reference_genome = None
        self.ref_genome_path = None
        self.snpdb_name = None
        self.tmp_dir = None
        self.vcf_filehandle = None
        self.depth_cutoff = None
        self.mq_cutoff = None
        self.ad_cutoff = None
        self.number_mixed_positions = None
        self.rec_list = []
        self.vcf_max_pos = None
        self.ref = None
        self.contig = None
        self.parsed_vcf_container = []
        self.mapper = None
        self.variant_caller = None
        self.variant_caller_threads = None
        self.mapper_threads = None

    def parse_config_dict(self, config_dict):
        # # we loop through thusly in case not all these things are in the config
        print (config_dict)
        for attr in config_dict:
            if attr == 'reference_genome':
                self.reference_genome = config_dict[attr]
            if attr == 'depth_cutoff':
                self.depth_cutoff = config_dict[attr]
            if attr == 'mq_cutoff':
                self.mq_cutoff = config_dict[attr]
            if attr == 'ad_cutoff':
                self.ad_cutoff = config_dict[attr]
            if attr == 'mapper':
                self.mapper = config_dict[attr]
            if attr == 'variant_caller':
                self.variant_caller = config_dict[attr]
            if attr == 'variant_caller_threads':
                self.variant_caller_threads = config_dict[attr] 
            if attr == 'mapper_threads':
                self.mapper_threads = config_dict[attr]                              
# -------------------------------------------------------------------------------------------------


    def mkdir_p(self, path):
        try:
            os.makedirs(path)
        except OSError as exc:
            if exc.errno == errno.EEXIST and os.path.isdir(path):
                pass
            else:
                raise
# -------------------------------------------------------------------------------------------------

    def make_tmp_dir(self, args):
        try:
            self.tmp_dir = os.path.join(os.path.dirname(args.fastqs[0]), 'snpdb')
            self.mkdir_p(self.tmp_dir)
        except AttributeError:
            # # this is because we also call this from vcf_to_db, where args doesn't contain fastqs, but a vcf.
            self.tmp_dir = os.path.join(os.path.dirname(args.vcf[0]))

# -------------------------------------------------------------------------------------------------

    def create_snippy_vcf(self, ref_seq, strain_seq, snpdb):
        snp_diff = []
        ig_diff = []
        nucs = ['A','T','C','G']


        ref_key = next(iter(ref_seq))
        strain_key = next(iter(strain_seq))

        parsed_vcf = ParsedVcf()
        parsed_vcf.ref = ref_key

        # Ensure both sequences have the same length
        if len(ref_seq[ref_key]) != len(strain_seq[strain_key]):
            raise ValueError("Both lists must have the same length.")

        # Iterate through the lists using enumerate to get the index and value
        for index, (value1, value2) in enumerate(zip(ref_seq[ref_key], strain_seq[strain_key])):
            if value1 != value2:
                # Store the index and both differing values in a tuple
                parsed_vcf.var[index+1] = value2
                parsed_vcf.ref_base[index+1] = value1 
                if value2 in nucs:
                    parsed_vcf.good_var.append(index+1)
                else:
                    parsed_vcf.bad_pos.append(index+1)
        if len(parsed_vcf.bad_pos) > 250000:
            print ("large amount of missing positions - %s" % (len(parsed_vcf.bad_pos)))
            sys.exit()
        if len(parsed_vcf.good_var) > 10000:
            print ("large amount of variants positions - %s" % (len(parsed_vcf.good_var)))
            sys.exit() 
        self.depth_average = 999           
        self.parsed_vcf_container.append(parsed_vcf)


# -------------------------------------------------------------------------------------------------


    def read_vcf_bcftools(self, ref_seq, snpdb):

        try:
            os.path.exists(self.vcf_filehandle)
        except IOError:
            print (self.vcf_filehandle + " not found ... ")

        if self.vcf_filehandle.endswith('.gz'):
            os.system('gunzip {0}'.format(self.vcf_filehandle))
            self.vcf_filehandle = os.path.splitext(self.vcf_filehandle)[0]
        try:
            openfile = open(self.vcf_filehandle, 'r')
        except:
            print (self.vcf_filehandle + " not found ... ")
            sys.exit()

        oref = ''
        parsed_vcf = ""


        for ref_name in ref_seq:
            parsed_vcf = ParsedVcf()
            parsed_vcf.ref = ref_name

            #get good variants
            os.system("bcftools query -i 'GT=\"1/1\" && MQ>%s && FORMAT/DP>%s && CHROM=\"%s\"' -f '%%CHROM\t%%POS\t%%REF\t%%ALT\n' %s > %s/snps.txt" % (snpdb.mq_cutoff, snpdb.depth_cutoff, ref_name, self.vcf_filehandle, self.tmp_dir))
            print("bcftools query -i 'GT=\"1/1\" && MQ>%s && FORMAT/DP>%s && CHROM=\"%s\"' -f '%%CHROM\t%%POS\t%%REF\t%%ALT\n' %s > %s/snps.txt" % (snpdb.mq_cutoff, snpdb.depth_cutoff, ref_name, self.vcf_filehandle, self.tmp_dir))

            openfile = open(self.tmp_dir+'/snps.txt', 'r')
            for line in openfile:
                temp = line.strip().split()
                pos = temp[1]
                ref_call = temp[2]
                var_call = temp[3]

                parsed_vcf.good_var.append(pos)
                parsed_vcf.var[pos] = var_call
                parsed_vcf.ref_base[pos] = ref_call                

            #get bad positions
            os.system("bcftools query -i 'GT=\"0/1\" || MQ<%s || FORMAT/DP<%s && CHROM=\"%s\"' -f '%%CHROM\t%%POS\t%%REF\t%%ALT\n' %s > %s/badpos.txt" % (snpdb.mq_cutoff, snpdb.depth_cutoff, ref_name, self.vcf_filehandle, self.tmp_dir))
            print("bcftools query -i 'GT=\"0/1\" || MQ<%s || FORMAT/DP<%s && CHROM=\"%s\"' -f '%%CHROM\t%%POS\t%%REF\t%%ALT\n' %s > %s/badpos.txt" % (snpdb.mq_cutoff, snpdb.depth_cutoff, ref_name, self.vcf_filehandle, self.tmp_dir))
  
            openfile = open(self.tmp_dir+'/badpos.txt', 'r')
            for line in openfile:
                temp = line.strip().split()
                pos = temp[1]
                ref_call = temp[2]
                var_call = temp[3]

                parsed_vcf.bad_pos.append(pos)
                parsed_vcf.var[pos] = var_call
                parsed_vcf.ref_base[pos] = ref_call 

            #get depth
            os.system("samtools depth -aa -q 13 -Q 30 -d 0 %s/%s.bam | awk '{ sum += $3; n++ } END { if (n > 0) print sum / n; }' > %s/depth.txt" % (self.tmp_dir, self.sample_name, self.tmp_dir))
            openfile = open(self.tmp_dir+'/depth.txt', 'r')
            for line in openfile:
                self.depth_average = line.strip()

            self.parsed_vcf_container.append(parsed_vcf)





# -----------------------------DEPRECIATED--------------------------------------------------------------------

    def read_multi_contig_vcf(self, ref_seq):

        try:
            os.path.exists(self.vcf_filehandle)
        except IOError:
            print (self.vcf_filehandle + " not found ... ")

        if self.vcf_filehandle.endswith('.gz'):
            os.system('gunzip {0}'.format(self.vcf_filehandle))
            self.vcf_filehandle = os.path.splitext(self.vcf_filehandle)[0]
        try:
            openfile = open(self.vcf_filehandle, 'r')
        except:
            print (self.vcf_filehandle + " not found ... ")
            sys.exit()

        oref = ''

        parsed_vcf = ""

        for line in openfile:
            if not re.match(r'^#', line):
                split_line = line.split()
                # get reference
                ref = split_line[0]

                if ref in ref_seq:
                    # if the reference is new
                    if ref != oref:
                        # if the reference is new and this is not the first reference and the parsed vcf object
                        if oref != '':
                            self.parsed_vcf_container.append(parsed_vcf)
                        #create a parsed vcf object
                        parsed_vcf = ParsedVcf()
                        oref = ref
                        parsed_vcf.ref = ref


                    #make some vars so easier to read
                    pos = split_line[1]
                    filter_flag = split_line[6]
                    var_call = split_line[4]
                    ref_call = split_line[3]

                    parsed_vcf.filter_flag[pos] = filter_flag

                    #split into ignore
                    if parsed_vcf.filter_flag[pos] != 'PASS':
                        parsed_vcf.bad_pos.append(pos)
                    else:
                        parsed_vcf.good_var.append(pos)

                    parsed_vcf.var[pos] = var_call
                    parsed_vcf.ref_base[pos] = ref_call

                    # get depth
                    matchObj = re.match(r'.*DP=(\d*)', line)
                    try:
                        parsed_vcf.depth[pos] = matchObj.group(1)
                    except AttributeError:
                        parsed_vcf.depth[pos] = 0
                    # get map quality
                    matchObj = re.match(r'.*MQ=(\d*\.\d*)', line)
                    try:
                        parsed_vcf.qual[pos] = matchObj.group(1)
                    except AttributeError:
                        parsed_vcf.qual[pos] = 0

                    #if a variant get some other things - including with it's a mix
                    if var_call != '.':
                        matchObvar = re.match(r'.*GT:AD:DP:GQ:PL\s+(.*)', line)
                        format_string = matchObvar.group(1).split(':')
                        parsed_vcf.hap_call[pos] = format_string[0]
                        parsed_vcf.hap_depth[pos] = format_string[2]
                        parsed_vcf.hap_qual[pos] = format_string[3]
                        ad_string = format_string[1].split(',')
                        parsed_vcf.hap_var_count[pos] = float(ad_string[1]) / float(parsed_vcf.hap_depth[pos])
                        if parsed_vcf.hap_var_count[pos] < self.ad_cutoff:
                            parsed_vcf.mixed_positions.append(int(pos))
                else:
                    sys.stderr.write('Contig not in reference genome found in config directory.\n')
                    sys.exit()
                                
            else:
                if re.match(r'##coverageMetaData',line):
                    matchObj = re.match(r'.*mean=(\d*\.\d*)',line)
                    try:
                        self.depth_average = matchObj.group(1)
                    except AttributeError:
                        self.depth_average = '0.0'

        #add the last vcf
        self.parsed_vcf_container.append(parsed_vcf)
        #close file
        openfile.close()

        #calculate total number of mixed positions // this has been depreciated
        self.number_mixed_positions = 0
        for p_vcf in self.parsed_vcf_container:
            self.number_mixed_positions += len(p_vcf.mixed_positions)


# -------------------------------------------------------------------------------------------------

    def parse_json_dict(self, json_dict, ref_seq):
        
        parsed_vcf = ""
        for ref_contig in ref_seq:
            parsed_vcf = ParsedVcf()
            parsed_vcf.ref = ref_contig
            for base in json_dict[ref_contig]:
                if base == 'N':
                    for pos in json_dict[ref_contig][base]:
                        parsed_vcf.bad_pos.append(pos)
                elif base != '-':
                    for var in json_dict[ref_contig][base]:
                        pos, ref_call = var.split(".")
                        parsed_vcf.good_var.append(pos)                       
                        parsed_vcf.var[pos] = base
                        parsed_vcf.ref_base[pos] = ref_call    

            self.parsed_vcf_container.append(parsed_vcf)

# -------------------------------------------------------------------------------------------------

    def make_ref_fastqs(self, args):
        try:
            if os.path.exists(os.path.join(snapperdb.__ref_genome_dir__, self.reference_genome + '.R1.fastq.gz')):
                sys.stderr.write('FASTQs found for  %s\n' % self.reference_genome)
            else:
                fastq_path1 = (os.path.join(snapperdb.__ref_genome_dir__, self.reference_genome + '.R1.fastq'))
                fastq_path2 = (os.path.join(snapperdb.__ref_genome_dir__, self.reference_genome + '.R2.fastq'))
                print (fastq_path1, fastq_path2, self.ref_genome_path)
                os.system('wgsim -e 0 -N 3000000 -1 100 -2 100 -r 0 -R 0 -X 0 %s %s %s' % (self.ref_genome_path,fastq_path1,fastq_path2))
                os.system('gzip %s' % (fastq_path1))
                os.system('gzip %s' % (fastq_path2))

        except IOError:
            sys.stderr.write('Error making reference FASTQ')
# -------------------------------------------------------------------------------------------------


    def define_fastq_paths(self,args):
        args.fastqs = []
        args.fastqs.append(os.path.join(snapperdb.__ref_genome_dir__, self.reference_genome + '.R1.fastq.gz'))
        args.fastqs.append(os.path.join(snapperdb.__ref_genome_dir__, self.reference_genome + '.R2.fastq.gz'))

# -------------------------------------------------------------------------------------------------

    def define_class_variables_and_make_output_files(self, args):
        try:
            self.sample_name = os.path.basename(args.vcf[0]).split(os.extsep)[0]
        except AttributeError:
            self.sample_name = os.path.basename(args.fastqs[0]).split(os.extsep)[0]
        try:
            if os.path.exists(os.path.join(snapperdb.__ref_genome_dir__, self.reference_genome + '.fa')):
                self.ref_genome_path = os.path.join(snapperdb.__ref_genome_dir__, self.reference_genome + '.fa')
            elif os.path.exists(os.path.join(snapperdb.__ref_genome_dir__, self.reference_genome + '.fasta')):
                self.ref_genome_path = (os.path.join(snapperdb.__ref_genome_dir__, self.reference_genome + '.fasta'))
            elif os.path.exists(os.path.join(snapperdb.__ref_genome_dir__, self.reference_genome)):
                self.ref_genome_path = (os.path.join(snapperdb.__ref_genome_dir__, self.reference_genome))
            else:
               sys.stderr.write('Cant find reference genome %s\n' % self.reference_genome)
               sys.exit()

        except IOError:
            sys.stderr.write('Cant find reference genome %s' % self.reference_genome)
        
        # set vcf path
        self.make_tmp_dir(args)
        try:
            self.vcf_filehandle = args.vcf[0]
        except:
            self.vcf_filehandle = os.path.join(self.tmp_dir, os.path.pardir, '{0}.vcf'.format(self.sample_name))
# -------------------------------------------------------------------------------------------------

    def check_reference_bwa_indexed(self):
        indices = ['amb', 'ann', 'bwt', 'pac', 'sa']
        for i in indices:
            if os.path.exists(self.ref_genome_path + '.' + i):
                pass
            else:
                os.system('bwa index %s' % self.ref_genome_path)
        for i in indices:
            if os.path.exists(self.ref_genome_path + '.' + i):
                pass
            else:
                sys.stderr.write('Reference not indexed, you need to index with `bwa index <ref_genome.fa>`')
                sys.exit()
# ---------------------------------------DEPRECIATED----------------------------------------------------------

    def run_phoenix(self,args):
        self.check_reference_bwa_indexed()
        self.check_reference_gatk_indexed()
        self.mapper = 'bwa'
        self.variant_caller = 'gatk'        
        os.system('phenix.py run_snp_pipeline -r1 %s -r2 %s -r %s -o %s -m %s -v %s --sample-name %s --filters mq_score:%s,min_depth:%s,ad_ratio:%s --annotators coverage nlessness' % (args.fastqs[0], args.fastqs[1], self.ref_genome_path, self.tmp_dir, self.mapper, self.variant_caller, self.sample_name, self.mq_cutoff, self.depth_cutoff,self.ad_cutoff))

# -------------------------------------------------------------------------------------------------



    def run_bwa(self,args):
        self.check_reference_bwa_indexed()
        print ("bwa mem -Y -M -R '@RG\\tID:%s\\tSM:%s' -t %s %s %s %s | samclip --max 10 --ref %s.fai | samtools sort -n -l 0 -T /tmp --threads 1 -m 4000M | samtools fixmate -m --threads 1 - - | samtools sort -l 0 --threads 1 -m 4000M | samtools markdup --threads 1 -r -s - - > %s/%s.bam" % (self.sample_name, self.sample_name, self.mapper_threads, self.ref_genome_path, args.fastqs[0], args.fastqs[1],self.ref_genome_path, self.tmp_dir, self.sample_name))

        os.system("bwa mem -Y -M -R '@RG\\tID:%s\\tSM:%s' -t %s %s %s %s | samclip --max 10 --ref %s.fai | samtools sort -n -l 0 -T /tmp --threads 1 -m 4000M | samtools fixmate -m --threads 1 - - | samtools sort -l 0 --threads 1 -m 4000M | samtools markdup --threads 1 -r -s - - > %s/%s.bam" % (self.sample_name, self.sample_name, self.mapper_threads, self.ref_genome_path, args.fastqs[0], args.fastqs[1],self.ref_genome_path, self.tmp_dir, self.sample_name))
        os.system("samtools index %s/%s.bam" % (self.tmp_dir, self.sample_name))


# -------------------------------------------------------------------------------------------------

    def run_gatk(self,args):
        self.check_reference_gatk_indexed()
        print ("java -XX:+UseSerialGC -jar GenomeAnalysisTK.jar -T UnifiedGenotyper -R %s -I %s/%s.bam --sample_ploidy 2 --genotype_likelihoods_model SNP -rf BadCigar -out_mode EMIT_ALL_SITES -nt %s > %s/%s.vcf" % (self.ref_genome_path, self.tmp_dir, self.sample_name, self.variant_caller_threads, self.tmp_dir, self.sample_name))
        os.system("java -XX:+UseSerialGC -jar GenomeAnalysisTK.jar -T UnifiedGenotyper -R %s -I %s/%s.bam --sample_ploidy 2 --genotype_likelihoods_model SNP -rf BadCigar -out_mode EMIT_ALL_SITES -nt %s > %s/%s.vcf" % (self.ref_genome_path, self.tmp_dir, self.sample_name, self.variant_caller_threads , self.tmp_dir, self.sample_name))


# -------------------------------------------------------------------------------------------------


    def check_reference_gatk_indexed(self):
        try:
            picard_path = "java -jar "+os.environ['PICARDPATH']
        except KeyError:
            picard_path = 'picard CreateSequenceDictionary'

        indicies = ['dict', 'fa.fai']
        # print (os.path.splitext(self.ref_genome_path)[0] + '.' + i)
        if os.path.exists(os.path.splitext(self.ref_genome_path)[0] + '.' + 'dict'):
            pass
        else:
            picard_dict_path = os.path.splitext(self.ref_genome_path)[0]
            os.system('%s R= %s O= %s.dict'
                      % (picard_path, self.ref_genome_path, picard_dict_path))

        if os.path.exists(os.path.splitext(self.ref_genome_path)[0] + '.' + 'fa.fai'):
            pass
        else:
            os.system('samtools faidx %s' % self.ref_genome_path)
        ## double check that above has worked
        for i in indicies:
            if os.path.exists(os.path.splitext(self.ref_genome_path)[0] + '.' + i):
                pass
            else:
                sys.stderr.write('Reference not indexed for GATK, you need to index with according to ''https://www.broadinstitute.org/gatk/guide/article?id=1601/n')
                sys.exit()
# -------------------------------------------------------------------------------------------------


    def read_rec_file(self, rec_file):
        try:
            openfile = open(rec_file, 'r')
        except:
            print (rec_file + " not found ... ")
            sys.exit()
        rec_list = []
        for line in openfile:
            if line[0].isdigit():
                temp = (line.strip()).split('\t')
                rec_range = range((int(temp[0]) - 1), (int(temp[1]) - 1))
                rec_list = set(rec_list) | set(rec_range)
        self.rec_list = rec_list
