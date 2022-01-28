from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import Entrez
import Bio
import textwrap
import csv
import copy
import time
import os
import argparse
import re

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description = 'Multigene blastn: pull all associated sequences from a blastn hit with the same specimen voucher',\
            epilog = textwrap.dedent('''
            Additional information:
            This script takes a fasta file of query sequences and blasts them against the nr database.
            It returns the csv file "collections_table.csv" containing all the all the blast matches and the associated hit data.
            If any additional sequences from the same specimen voucher are available then that accession number is included.
            A fasta file "blast_match_sequences.fasta" containing all of the blast matches and additional markers from the same
            specimen is written to the working directory. The first field in the fasta header is the same for each species collection,
            allowing different markers to be concatenated after sequence alignment.'''))

parser.add_argument('-i', '--input', type = str, metavar = '', required = True, help = 'enter input fasta file name')
parser.add_argument('-e', '--email', type = str, metavar = '', required = True, help = 'enter email address so NCBI can contact you if there is a problem')
parser.add_argument('-v', '--verbose', action = 'store_true', help = 'verbose, keep intermediate files')
args = parser.parse_args()


Entrez.email = args.email


def blastnncbi(fasta_file):
    '''
    Takes a fasta file of sequence queries and blasts them against ncbi nt database,
    excluding all environmental sequences. Writes the resulting file to Temp_File_blast_results.xml file
    '''
    print('Retrieving results from blastn')
    fasta_file = open(fasta_file).read()
    result = NCBIWWW.qblast('blastn', 'nt', fasta_file, entrez_query='NOT ("environmental samples"[All Fields] OR "metagenomes"[All Fields] OR "uncultured"[All Fields])')
    blast_file = open('Temp_File_blast_results.xml', 'w')
    blast_file.write(result.read())
    blast_file.close()
    result.close()
    print('blastn results retrieved')


def parse_xml_hits(open_xml):
    '''
    Takes the output from blastnncbi and returns a dictionary of accession numbers
    with the corresponding blast match percentage identity to the query sequence as
    a dictionary.
    '''

    open_xml = open(open_xml)
    blast_records = NCBIXML.parse(open_xml)
    number_records = 0
    percent_id = {}

    for blast_record in blast_records:
        number_records += 1
        accession = ''
        percent = 0

        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                #print(alignment.title)
                accession = alignment.accession
                percent = round(float(hsp.positives)/float(hsp.align_length)*100, 2)
                align_length = alignment.length

                #check if the accesion is already in percent_id and only update if the blast match is higher and it is longer

                if accession not in percent_id.keys():
                    percent_id[accession] = {'pairwise, align_length, query':'{}, {}, {}'.format(percent,align_length,number_records)}
                elif accession in percent_id.keys() and float(percent_id[accession]['pairwise, align_length, query'].split(',')[0]) < percent and float(percent_id[accession]['pairwise, align_length, query'].split(',')[1]) < align_length:
                    percent_id[accession] = {'pairwise, align_length, query':'{}, {}, {}'.format(percent,align_length,number_records)}

    print('Extracted blastn hit accessions')
    return percent_id


def txid_for_hits():
    '''
    From the keys in the dictionary produced from parse_xml_hits this script pulls
    all the txids and returns a list of a txids.
    '''
    list_of_ids = parse_xml_hits('Temp_File_blast_results.xml')
    list_of_ids = list(list_of_ids.keys())
    list_of_txids = []
    #get related sequeces from GB
    record_num = Entrez.read(Entrez.elink(db='taxonomy', dbfrom='nuccore', id=list_of_ids))

    #make list of txids
    for i in range(len(record_num)):
        txid = record_num[i]['LinkSetDb'][0]['Link'][0]['Id']
        if txid not in list_of_txids:
            list_of_txids.append(txid)

    print('Retrieved txids for blastn hit accessions')
    return list_of_txids


def fetch_seqs_from_txid():
    '''
    Takes a list of txids from txid_for_hits() and pulls down all associated info
    for each txid. From each txid record it adds the associated sequence id to the
    sequence_accession list. This list is then used to download all the sequence
    info from genbank and is written to the file 'Temp_File_taxa_hits_file.gb' for later parsing.
    I had to do the batch downloading as a gb file and not xml. For whatever reason, I
    cannot parse the resulting batch downloaded xml file with biopyhon. It throws a
    corrupt file error at the position of where the first batch is written to and the
    next batch starts.
    '''

    txid_list = txid_for_hits()
    sequence_search_handle = Entrez.elink(dbfrom='taxonomy', db='nuccore', id=txid_list)
    sequence_search_results = Entrez.read(sequence_search_handle)
    sequence_search_handle.close()

    sequence_accession = []

    # get sequence accession from the txid entry
    for i in range(len(sequence_search_results)):
        range_txid = range(len(sequence_search_results[i]['LinkSetDb'][0]['Link']))
        for j in range_txid:
            seq_GI = sequence_search_results[i]['LinkSetDb'][0]['Link'][j]['Id']
            if seq_GI not in sequence_accession:
                sequence_accession.append(seq_GI)


    # there is a pesky nuccore entry that doesn't belong in that database. If it is not removed
    #then the gb parser fails

    #if '1562123491' in sequence_accession:
    #    sequence_accession.remove('1562123491')
    #elif '1796269881' in sequence_accession:
    #    sequence_accession.remove('1796269881')


    #Use epost to get webenv and query_key to speed up retrieving records
    epost_var = Entrez.epost("nuccore", id=",".join(sequence_accession))
    epost_search_results = Entrez.read(epost_var)
    webenv = epost_search_results["WebEnv"]
    query_key = epost_search_results["QueryKey"]

    #get the number of records for batch downloading
    count = len(sequence_accession)

    batch_size = 200
    out_handle = open("Temp_File_taxa_hits_file.gb", "w")
    for start in range(0, count, batch_size):
        end = min(count, start+batch_size)
        fetch_handle = Entrez.efetch(db="nucleotide", id=sequence_accession, rettype="gbwithparts", retmode="text", retstart=start, retmax=batch_size, webenv=webenv, query_key=query_key) #rettype="gb"
        data = fetch_handle.read()
        fetch_handle.close()
        out_handle.write(data)
    print('Retrieved associated GB files from txids')
    out_handle.close()



def parse_gb_full_files():
    '''
    This function reads in Temp_File_taxa_hits_file.gb and parses out the species, voucher
    sequence etc. It returns a dictionary only containing the species_col and the
    associated region:accessions plus the blast match info.
    '''
    handle = open("Temp_File_taxa_hits_file.gb")
    record_iterator = SeqIO.parse(handle, "gb")
    count = 0
    no_voucher_num = 0
    dic_gb = {}

    for record in record_iterator:
        try:
            count += 1
            organism = ''
            accession = record.name
            isolate = ''
            specimen_voucher = ''
            sequence = record.seq
            gene = ''
            header_gene = ''
            gene_for_table = ''
            header = record.description.lower()
            strain = ''
            clone = ''
            species_col = ''

            #Check if the sequence is ITS, LSU or SSU
            if 'internal transcribed spacer' in header or 'its1' in header or 'its2' in header or 'its region' in header:
                header_gene = 'ITS'
            elif 'external transcribed spacer' in header:
                header_gene = 'ETS'
            elif '28s' in header or 'large subunit' in header or '25s' in header and 'internal transcribed spacer' not in header:
                header_gene = 'LSU'
            elif '18s' in header or 'small subunit' in header and 'internal transcribed spacer' not in header and '5.8s' not in header:
                header_gene = 'SSU'
            elif 'psbA-trnH' in header:
                header_gene = 'psbA-trnH'
            elif 'psbK-psbI' in header:
                header_gene = 'psbK-psbI'
            elif 'atpF-atpH' in header:
                header_gene = 'atpF-atpH'
            else:
                header_gene = 'unknown'

            #Check if the sequence is from a type specimen, often contained in the header
            if 'type' in header:
                header_gene = header_gene + '-TYPE'

            #Often the ITS and LSU regions are concatenated together. Make a length cut off of
            #700 bp and try and account for some length variation by making the min len difference of 100 (800 total length)
            elif header_gene == 'ITS' and len(sequence) > 700 and len(sequence) - 700 >= 100:
                header_gene = 'ITS and LSU?'
            else:
                header_gene = header_gene

            #Try and catch cases where the collection is undescribed and also has the voucher in the
            #organism line e.g. Cortinarius sp. OTA00001
            if record.annotations['organism'].split(' ')[1] == 'sp.':
                organism = record.annotations['organism'].split(' ')
                organism = organism[:2]
                organism = '_'.join(organism).replace('.', '')
            #replace white space with '_'
            elif record.annotations['organism'].split(' ')[1] != 'sp.':
                organism = record.annotations['organism'].replace(' ', '_')


            #Find some sort of voucher and gene
            for feature in record.features:
                for key, value in feature.qualifiers.items():
                    if 'isolate' in key:
                        isolate = ''.join(value)
                    elif 'specimen_voucher' in key:
                        #[David] specimen voucher are coded as "XXXXX/XX" by convention, but
                        #some people put it as "XXXXX_XX". So, change the last "_" to "\".
                        specimen_voucher_temp = ''.join(value)
                        specimen_voucher = re.sub(r"(\S+)_(\S+)$", r"\1/\2", specimen_voucher_temp)
                    elif 'strain' in key:
                        strain = ''.join(value)
                    elif 'clone' in key:
                        clone = ''.join(value)
                    #Check what the gene is. If the sequence is from a single gene with multiple exons
                    #then only the gene name will taken. If there are multiple genes present in the sequence
                    #then the different gene names will be concatenated together e.g. trnL-trnF
                    elif 'gene' in key and len(gene) == 0:
                        gene = ''.join(value)
                    elif 'gene' in key and ''.join(value) == gene:
                        gene = gene
                    elif 'gene' in key and str(gene).find(''.join(value)) == -1:#len(gene) != 0:
                        gene = gene + str('-' + ''.join(value))

            #Check what vouchers are present and concatenate organism name and voucher together
            #the specimen voucher has presidence over isolate, strain, clone
            if len(specimen_voucher) != 0:
                species_col = organism + '__' + specimen_voucher
            elif len(specimen_voucher) == 0 and len(isolate) != 0:
                species_col = organism + '__' + isolate
            elif len(specimen_voucher) == 0 and len(isolate) == 0 and len(strain) != 0:
                species_col = organism + '__' + strain
            elif len(specimen_voucher) == 0 and len(isolate) == 0 and len(strain) == 0 and len(clone) != 0:
                species_col = organism + '__' + clone
            elif len(specimen_voucher) == 0 and len(isolate) == 0 and len(strain) == 0 and len(clone) == 0:
                no_voucher_num += 1
                species_col = organism + '__no_voucher_{}'.format(str(no_voucher_num))

            #replace all spaces in collection names, otherwise it will cause issues when the fasta file is imported into other software
            species_col = species_col.replace(' ', '-')

            #check that ITS not in header
            if len(gene) != 0: #and header_gene == 'unknown'
                gene_for_table = gene
            else:
                gene_for_table = header_gene

            #add species_col to dic_gb with the associated regions and accessions
            if species_col not in dic_gb:
                dic_gb[species_col] = {gene_for_table:accession}
            elif species_col in dic_gb:
                dic_gb[species_col][gene_for_table] = accession

            #create a fasta file of all the sequences from the .gb file
            with open('Temp_File_temp_all_GB_sequences.fasta','a+') as fastafile:
                if len(sequence) > 0:
                    fastafile.write('>{} {} {}\n{}\n'.format(species_col, accession, gene_for_table, str(sequence)))
        except Bio.Seq.UndefinedSequenceError:
            break


    #check which blast matches are in dic_gb against blast_matches dictionary and store matches
    #with their regions:accessions and blast match data in new_dic dictionary

    blast_matches = parse_xml_hits('Temp_File_blast_results.xml')

    new_dic = {}

    for key, value in dic_gb.items():
        for k, v in value.items():
            for blast_key, blast_value in blast_matches.items():
                if blast_key in v:
                    new_dic[key] = value

                    value_tmp = copy.deepcopy(value)
                    value_tmp.update(blast_value)
                    new_dic[key] = value_tmp


    print('Parsed GB {} files'.format(count))
    return new_dic


#print(parse_gb_full_files())


def write_table_csv_fasta():
    '''
    Take the dictionary from parse_gb_full_files and writes a csv file with the blast matches
    Also writes a fasta file with all the sequeces from the blast results.
    '''

    ##check that temp_fasta_file.fasta or blast_match_sequences.fasta and remove if they do
    ##otherwise they will be appended to and create duplicate entries
    dir_name = os.getcwd()
    file_list = os.listdir(dir_name)
    for item in file_list:
        if item == 'Temp_File_temp_all_GB_sequences.fasta' or item == 'blast_match_sequences.fasta':
            os.remove(os.path.join(dir_name, item))

    #get dictionary from parse_gb_full_files and write to csv
    dict_data = parse_gb_full_files()
    csv_fields = ['Species and collection']
    for species_col, regions in dict_data.items():
        for region in regions.keys():
            if region not in csv_fields:
                csv_fields.append(region)

    file = open('collections_table.csv', 'w')
    writer = csv.DictWriter(file, fieldnames = csv_fields, extrasaction = 'ignore')
    writer.writeheader()
    for key,value in sorted(dict_data.items()):
        row = {'Species and collection':key}
        row.update(value)
        writer.writerow(row)
    file.close()

    #read temp_all_GB_sequences and write all the blast match sequences to blast_match_sequences.fasta
    temp_fasta_file = open('Temp_File_temp_all_GB_sequences.fasta', 'r')

    accession_list = []

    for species_col, regions in dict_data.items():
        for accession in regions.values():
            if accession[0].isalpha():
                accession_list.append(accession)

    for sequence in SeqIO.parse('Temp_File_temp_all_GB_sequences.fasta', 'fasta'):
        for accession in accession_list:
            if accession in sequence.description:
                with open('blast_match_sequences.fasta', 'a+') as fasta_file:
                    fasta_file.write('>'+str(sequence.description + '\n'))
                    fasta_file.write(str(sequence.seq + '\n'))


def remove_temp_files():
    dir_name = os.getcwd()
    test = os.listdir(dir_name)
    for item in test:
        if item.startswith("Temp_File_"):
            os.remove(os.path.join(dir_name, item))



if __name__ == '__main__':
    if args.verbose:
        start = time.time()
        blastnncbi(args.input)
        fetch_seqs_from_txid()
        write_table_csv_fasta()
        finish = time.time()
        sec_diff = round(finish - start, 1)
        min_dif = round(sec_diff/60, 1)
        print('Time in seconds {}, or {} minutes'.format(sec_diff, min_dif))
    else:
        start = time.time()
        blastnncbi(args.input)
        fetch_seqs_from_txid()
        write_table_csv_fasta()
        remove_temp_files()
        finish = time.time()
        sec_diff = round(finish - start, 1)
        min_dif = round(sec_diff/60, 1)
        print('Time in seconds {}, or {} minutes'.format(sec_diff, min_dif))
