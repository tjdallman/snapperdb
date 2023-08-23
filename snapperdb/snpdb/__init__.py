__author__ = 'gidis'

import datetime
import inspect
import os
import re
import subprocess
import sys
import json
import logging
import psycopg2, psycopg2.extras
import snapperdb
from snapperdb.snpdb.snpdb import SNPdb
from snapperdb.gbru_vcf.vcf import Vcf
from snapperdb import parse_config
import pprint
import time
import math
import concurrent.futures as cf


# -------------------------------------------------------------------------------------------------

def snippy_to_db(args, config_dict):
    logger = logging.getLogger('snapperdb.snpdb.snippy_to_db')
    logger.info('Initialising SNPdb class')
    #create snpdb class
    snpdb = SNPdb(config_dict)

    #parse config into snpdb object
    logger.info('Parsing config dict')
    snpdb.parse_config_dict(config_dict)

    #connect to snpdb postgres
    snpdb._connect_to_snpdb()
    snpdb.snpdb_conn = psycopg2.connect(snpdb.conn_string)

    #read the reference fasta
    ref_seq_file = os.path.join(snapperdb.__ref_genome_dir__, snpdb.reference_genome + '.fa')
    ref_seq = read_multi_contig_fasta(ref_seq_file)

    #read_strain_file
    strain_seq = read_multi_contig_fasta(args.aln)

    #compare_aln_and_add_to_vcf_object
    vcf = Vcf()
    vcf.sample_name = args.name
    vcf.create_snippy_vcf(ref_seq, strain_seq, snpdb)

    logger.info('Uploading to SNPdb')
    #upload vcf
    snpdb.snpdb_upload(vcf,args)
    #annotate vars
    logger.info('Annotating new variants')

    snpdb.snpdb_annotate_vars(vcf)

def vcf_to_db(args, config_dict, vcf):
    #set up loggging
    logger = logging.getLogger('snapperdb.snpdb.vcf_to_db')
    logger.info('Initialising SNPdb class')

    #create snpdb class
    snpdb = SNPdb(config_dict)

    #parse config into snpdb object
    logger.info('Parsing config dict')
    snpdb.parse_config_dict(config_dict)

    #connect to snpdb postgres
    snpdb._connect_to_snpdb()
    snpdb.snpdb_conn = psycopg2.connect(snpdb.conn_string)

    #check stack?
    if inspect.stack()[0][3] == 'fastq_to_db':
        # fastq_to_db we will alread have a vcf object to work wih
        logger.info('You are running fastq_to_db.')

    elif inspect.stack()[0][3] == 'vcf_to_db':
        ## there is no existing vcf class here, but there will definitely be a vcf
        logger.info('You are running vcf_to_db. Initialising Vcf class.')
        vcf = Vcf()
        logger.info('Making SNPdb variables and output files')
        #set up variables
        snpdb.define_class_variables_and_make_output_files(args, vcf)

    #read the reference fasta
    ref_seq_file = os.path.join(snapperdb.__ref_genome_dir__, snpdb.reference_genome + '.fa')
    ref_seq = read_multi_contig_fasta(ref_seq_file)

    #read vcf
    #vcf.read_multi_contig_vcf(ref_seq)

    vcf.read_vcf_bcftools(ref_seq, snpdb)

    logger.info('Uploading to SNPdb')
    #upload vcf
    snpdb.snpdb_upload(vcf,args)
    #annotate vars
    logger.info('Annotating new variants')

    snpdb.snpdb_annotate_vars(vcf)

# -------------------------------------------------------------------------------------------------

def add_ref_cluster(args, config_dict):
    #set up loggging
    logger = logging.getLogger('snapperdb.snpdb.add_ref_cluster')
    logger.info('Initialising SNPdb class')

    #create snpdb class
    snpdb = SNPdb(config_dict)

    #parse config into snpdb object
    logger.info('Parsing config dict')
    snpdb.parse_config_dict(config_dict)

    #connect to snpdb postgres
    snpdb._connect_to_snpdb()
    snpdb.snpdb_conn = psycopg2.connect(snpdb.conn_string)

    snpdb.add_cluster()

# -------------------------------------------------------------------------------------------------

def make_snpdb(config_dict):
    snpdb = SNPdb(config_dict)
    snpdb._connect_to_snpdb()
    snpdb.make_snpdb()

# -------------------------------------------------------------------------------------------------

def read_file(file_name):
    #read list of strains for get the snps
    try:
        openfile = open(file_name, 'r')
    except:
        print ("### Strain list '" + file_name + "' not found ... ")
        print ("### Exiting "+str(datetime.datetime.now()))
        sys.exit()
    strain_list = []
    for line in openfile:
        strain_list.append(line.strip())
    return strain_list

# -------------------------------------------------------------------------------------------------

def read_multi_contig_fasta(ref):
    try:
        openfile = open(ref, 'r')
    except:
        print ("### Reference genome "+ ref + " not found ... ")
        print ("### Exiting "+str(datetime.datetime.now()))
        sys.exit()

    ref_seq = {}
    contig = ""
    for line in openfile:
        matchObj = re.search('>', line)
        if matchObj is None:
            for n in line.strip():
                ref_seq[contig[0]].append(n)
        else:
            contig = line[1:].strip().split()
            ref_seq[contig[0]] = []
    return ref_seq

# -------------------------------------------------------------------------------------------------

def read_rec_file_mc(rec_file):
    #read recombination file - tab delineated with contig
    try:
        openfile = open(rec_file, 'r')
    except:
        print ("### Recombination list "+ rec_file + " not found ... ")
        print ("### Exiting "+str(datetime.datetime.now()))
        sys.exit()

    rec_dict = {}
    for line in openfile:
        split_line = line.strip().split('\t')
        if split_line[0] in rec_dict:
            try:
                rec_dict[split_line[0]] += range(int(split_line[1]), (int(split_line[2])))
            except:
                print ('Error parsing ignored position list')
                sys.exit()
        else:
            rec_dict[split_line[0]] = []
            try:
                rec_dict[split_line[0]] += range(int(split_line[1]), (int(split_line[2])))
            except:
                print ('Error parsing ignored position list')
                sys.exit()
    #print rec_dict            
    return rec_dict

# -------------------------------------------------------------------------------------------------

def create_contig_index_for_consensus_genome(reference_genome):
    ## first, need to concatenate the reference genome in the order of
    ## alphanumerically sorted contig names
    concatenated_ref_genome = ''
    contig_names = sorted(reference_genome.keys())
    for c in contig_names:
        concatenated_ref_genome += ''.join(reference_genome[c])
    ## contig dict is going to be {(contig start in concat ref genome, contig stop in concat ref genome):contig_name}
    contig_index = {}
    ## start at 0
    i = 0
    ## need to iterate through the contig names in the sorted order, because this is the order that get_the_snps will go through them in
    for contig in contig_names:
        ## use i and j to move through ref genome
        ## we minus 1 here because of the python thing of not counting the last 'fence post' in the index, but counting it in len(). I think. Seems to work.
        j = i + len(reference_genome[contig]) - 1
        contig_index[(i, j)] = contig
        i += len(reference_genome[contig])

    pprint.pprint(contig_index)
    return contig_index

# -------------------------------------------------------------------------------------------------

def make_recomb_dict_from_gubbins(recombinant_sections, contig_index):
    rec_dict = {}
    for rs in recombinant_sections:
        for contig in contig_index:
            if contig[0] <= rs[0] <= contig[1]:
                if contig[0] <= rs[1] <= contig[1]:
                    recomb_start = (rs[0] - contig[0])
                    recomb_stop = (rs[1] - contig[0])
                    if contig_index[contig] in rec_dict:
                        rec_dict[contig_index[contig]] += range(recomb_start, recomb_stop)
                    else:
                        rec_dict[contig_index[contig]] = []
                        rec_dict[contig_index[contig]] += range(recomb_start, recomb_stop)
                else:
                    print ('The recombination section is not on one contig, this may not be a real recombination event, this is not currently being exlcuded.', rs, contig, contig_index[contig])
    for contig in rec_dict:
        rec_dict[contig] = set(rec_dict[contig])
    return rec_dict

# -------------------------------------------------------------------------------------------------

def read_rec_file_mc_gubbins(gubbins_rec_file, reference_genome):
    contig_index = create_contig_index_for_consensus_genome(reference_genome)
    try:
        openfile = open(gubbins_rec_file, 'r')
    except IOError:
        print ("### Gubbins GFF file "+ gubbins_rec_file + " not found ... ")
        print ("### Exiting "+str(datetime.datetime.now()))
        sys.exit()
    recombinant_sections = []
    for line in openfile.readlines():
        if not line.startswith('##'):
            split_line = line.strip().split('\t')
            recombinant_sections.append((int(split_line[3]), int(split_line[4])))

    rec_dict = make_recomb_dict_from_gubbins(recombinant_sections, contig_index)
    return rec_dict

 # -------------------------------------------------------------------------------------------------


def get_the_snps(args, config_dict):
    #set up logging
    logger = logging.getLogger('snapperdb.snpdb.get_the_snps')
    logger.info('Inititialising SnpDB Class')
    #initalise snpdb class
    snpdb = SNPdb(config_dict)
    #parse confif
    snpdb.parse_config_dict(config_dict)
    #read strainlist
    strain_list = read_file(args.strain_list)

    #connect to postgresdb
    snpdb._connect_to_snpdb()
    #get reference genome path
    ref_seq_file = os.path.join(snapperdb.__ref_genome_dir__, snpdb.reference_genome + '.fa')
    #read the reference fasta
    ref_seq = read_multi_contig_fasta(ref_seq_file)
    #add reference genome to strain_list
    strain_list.append(snpdb.reference_genome)

    #if recombination flag set
    if args.rec_file != 'N':
        logger.info('Reading recombination list')
        rec_dict = read_rec_file_mc(args.rec_file)
    elif args.gubbins_rec_file != None:
        logger.info('Reading gubbins recombination list')
        rec_dict = read_rec_file_mc_gubbins(args.gubbins_rec_file, ref_seq)
    else:
        #should we set this as none
        rec_dict = {}


    #query snadb
    snpdb.parse_args_for_get_the_snps_mc(args, strain_list, ref_seq, snpdb.reference_genome, rec_dict)

    snpdb.print_fasta_mc(args, rec_dict, ref_seq)

    #print matrix
    #if args.mat_flag == 'Y':
    #    snpdb.print_matrix(args.out)
    # print variant list
    if args.var_flag == 'Y':
        logger.info('Printing variants')
        snpdb.print_vars_mc(args,rec_dict)

# -------------------------------------------------------------------------------------------------

def export_json(args, config_dict):
    #set up logging
    logger = logging.getLogger('snapperdb.snpdb.export_json')
    logger.info('Inititialising SnpDB Class')
    #initalise snpdb class
    snpdb = SNPdb(config_dict)
    #parse confif
    snpdb.parse_config_dict(config_dict)
    #read strainlist
    strain_list = read_file(args.strain_list)

    #connect to postgresdb
    snpdb._connect_to_snpdb()
    #get reference genome path
    ref_seq_file = os.path.join(snapperdb.__ref_genome_dir__, snpdb.reference_genome + '.fa')
    #read the reference fasta
    ref_seq = read_multi_contig_fasta(ref_seq_file)


    snpdb.parse_args_for_export(args, strain_list, ref_seq)


  # -------------------------------------------------------------------------------------------------

def ignore_isolate(args,config_dict):

    #set up logging
    logger = logging.getLogger('snapperdb.snpdb.export_json')
    logger.info('Inititialising SnpDB Class')
    #initalise snpdb class
    snpdb = SNPdb(config_dict)
    #parse confif
    snpdb.parse_config_dict(config_dict)

    #connect to postgresdb
    snpdb._connect_to_snpdb()
    snpdb.snpdb_conn = psycopg2.connect(snpdb.conn_string)  

    #remove isolate
    snpdb.remove_isolate(args.ig_strain)


  # -------------------------------------------------------------------------------------------------

def accept_outlier(args,config_dict):

    #set up logging
    logger = logging.getLogger('snapperdb.snpdb.export_json')
    logger.info('Inititialising SnpDB Class')
    #initalise snpdb class
    snpdb = SNPdb(config_dict)
    #parse confif
    snpdb.parse_config_dict(config_dict)

    #connect to postgresdb
    snpdb._connect_to_snpdb()
    snpdb.snpdb_conn = psycopg2.connect(snpdb.conn_string)  

    #remove isolate
    snpdb.zscore_exception(args.out_strain)


  # -------------------------------------------------------------------------------------------------

def import_json(args):
    #set up logging
    logger = logging.getLogger('snapperdb.snpdb.import_json')
    json_path= untar_file(args.json_file)
    json_dict = {}

    #import json
    try:
        with open (json_path) as json_data:
            json_dict = json.load(json_data)
        json_data.close()
    except IOError:
        print ("Issue with JSON file. Exiting.")
        exit()

    #parse config
    args.config_file = json_dict['config_file']
    config_dict = parse_config(args)
    
    #initalise snpdb class
    snpdb = SNPdb(config_dict)
    #parse confif
    snpdb.parse_config_dict(config_dict)
    #get reference genome path
    ref_seq_file = os.path.join(snapperdb.__ref_genome_dir__, snpdb.reference_genome + '.fa')
    #read the reference fasta
    ref_seq = read_multi_contig_fasta(ref_seq_file)

    #create VCF class
    vcf = Vcf()
    vcf.parse_json_dict(json_dict, ref_seq)
    vcf.depth_average = json_dict['strain_stats']
    vcf.sample_name = json_dict['sample']

    logger.info('Uploading to SNPdb')
    #upload vcf
    #connect to snpdb postgres
    snpdb._connect_to_snpdb()
    snpdb.snpdb_conn = psycopg2.connect(snpdb.conn_string)    
    
    if args.write_flag == 'W':

        snpdb.snpdb_upload(vcf,args)
        #annotate vars
        logger.info('Annotating new variants')

        snpdb.snpdb_annotate_vars(vcf)
    elif args.write_flag == 'R':
        print ("under development")
        #snpdb.snpdb_query(vcf,args)

# -------------------------------------------------------------------------------------------------

def untar_file(file_path):
    if not file_path.endswith(".tar.gz"):
        print ("Not a Tar file. Exiting.")
        exit()
    bash_command = "tar -xzvf " + file_path
    flag = subprocess.call(bash_command, shell = True)
    if flag != 0:
        print ("Issue with Tar file. Exiting.")
        exit()
    new_path = file_path[0:-7]
    return new_path

# -------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------
def chunks(l, n):
    for i in range(0, len(l), n):
        yield l[i:i + n]

def update_distance_matrix(config_dict, args):
    logger = logging.getLogger('snapperdb.snpdb.update_distance_matrix')
    logger.info('Inititialising SnpDB Class')
    snpdb = SNPdb(config_dict)
    snpdb.parse_config_dict(config_dict)
    snpdb._connect_to_snpdb()
    logger.info('Getting strains')
    strain_list, update_strain = snpdb.get_strains()
    # # get_all_good_ids from snpdb2 takes a snp cutoff as well, here, we don't have a SNP cutoff so we set it arbitrarily high.
    snp_co = '1000000'
    if update_strain:
        print ("###  Populating distance matrix: " + str(datetime.datetime.now()))
        snpdb.parse_args_for_update_matrix(snp_co, strain_list)
        print (args.threads)
        if args.hpc == 'N' and args.threads == 'N':
            print ('### Launching serial update_distance_matrix ' + str(datetime.datetime.now()))
            snpdb.check_matrix(strain_list, update_strain)
            #snpdb.update_clusters(accept,place,full)
        elif args.hpc == 'N':
            print ('### Launching multi-processed update_distance_matrix ' + str(datetime.datetime.now()))
            
            present_strains = list(set(strain_list) - set(update_strain))
            chunksize = math.ceil(len(update_strain)/int(args.threads))
            
            #snpdb.shared_ig_pos()
            #snpdb.get_ignored_pos(strain_list)

            #remove connection before mp
            snpdb.snpdb_conn = None

            futures = []
            newrows = []
            #mp the comparison against the database
            with cf.ProcessPoolExecutor(max_workers=int(args.threads)) as executor:
                for idx, one_strain in enumerate(chunks(list(update_strain), chunksize)):
                    futures.append(executor.submit(snpdb.check_matrix_mp, present_strains, one_strain))
                for f in futures:
                    newrows = newrows + f.result()
            executor.shutdown()


            snpdb._connect_to_snpdb()
            snpdb.add_new_rows(newrows)

            futures = []
            newrows = []

            #remove connection before mp
            snpdb.snpdb_conn = None
            #mp the comparison against themselves
            with cf.ProcessPoolExecutor(max_workers=int(args.threads)) as executor:
                for idx, one_strain in enumerate(chunks(list(update_strain), chunksize)):
                    #print (one_strain, update_strain)
                    futures.append(executor.submit(snpdb.check_matrix_mp, update_strain, one_strain))
                    update_strain = set(update_strain) - set(one_strain)
                for f in futures:
                    newrows = newrows + f.result()

            snpdb._connect_to_snpdb()
            snpdb.add_new_rows(newrows)

        else:
            print ('### Launching parallel update_distance_matrix ' +str(datetime.datetime.now()))
            present_stains = list(set(strain_list) - set(update_strain))
            jobs_id_list = []
            for idx, one_strain in enumerate(chunks(list(update_strain), int(args.hpc))):
                jobs_id_list.extend(snpdb.write_qsubs_to_check_matrix(args, idx, one_strain, present_stains, config_dict['snpdb_name']))
            print ("start waiting: ",  datetime.datetime.now())
            _wait_for_jobs(jobs_id_list)
            print ("Finished waiting: ",  datetime.datetime.now())
            snpdb.check_matrix(update_strain, update_strain)
        uniq_dict = snpdb.get_table_unique()
        if any(uniq_dict.values()):
            print ('INCONSITENT TABLE DETECTED! ',  datetime.datetime.now())
            for k in uniq_dict:
                print ('\t', k, '\n','\t'*2, str(uniq_dict[k]))
    else:
        print ('### Nothing to update ' + str(datetime.datetime.now()))
# -------------------------------------------------------------------------------------------------

def qsub_to_check_matrix(config_dict, args):

    snpdb = SNPdb(config_dict)
    snpdb.parse_config_dict(config_dict)
    snpdb._connect_to_snpdb()
    snp_co = '1000000'
    added_list = []
    with open(args.added_list) as fi:
        for x in fi.readlines():
            added_list.append(x.strip())
    present_strains = []
    with open(args.present_strains) as fi:
        for x in fi.readlines():
            present_strains.append(x.strip())

    strain_list = list(set(present_strains) | set(added_list))
    snpdb.parse_args_for_update_matrix(snp_co, strain_list)
    snpdb.check_matrix(strain_list,added_list)


# -------------------------------------------------------------------------------------------------


def update_clusters(accept,place,full,config_dict):
    snpdb = SNPdb(config_dict)
    snpdb.parse_config_dict(config_dict)
    snpdb._connect_to_snpdb()
    snpdb.update_clusters(accept,place,full)



# -------------------------------------------------------------------------------------------------

def _wait_for_jobs(jobs_id_list):

    '''
        pre     :   jobs_id_list - list of job ID strings.
        post    :   Pauses execution of current function until all jobs in jobs_id_list have finished.
    '''
    
    do_exit = False
    while not do_exit:
        p = subprocess.Popen("qstat", shell=True, stdin=None, stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
        output_lines = [i.strip() for i in p.stdout.readlines()]
        current_jobs = [s.split(' ')[0] for s in output_lines]
        if any([i in current_jobs for i in jobs_id_list]):
            do_exit = False
        else:
            do_exit = True
        time.sleep(60) 
    return 0
# -------------------------------------------------------------------------------------------------


def get_nlessness(args, config_dict):

    '''
        pre     :   config_dict - path to config file for ebg of interest
                    args - subparser for nlessness determination 
        post    :   calls get_nlessness in snpdb.
    '''

    if args.strain_list:
        with open(args.strain_list, 'r') as strains_input:
            strain_list = strains_input.readlines()
            strain_list = [i.strip() for i in strain_list]
            strain_list = [i for i in strain_list if i != '']
    else: 
        strain_list = None
    snpdb = SNPdb(config_dict)
    snpdb._connect_to_snpdb()
    snpdb.get_nlessness(args.name, strain_list, args.all_isolates, args.print_to_csv)
    
# -------------------------------------------------------------------------------------------------





