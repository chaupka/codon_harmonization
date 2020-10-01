import re
from functions import get_codon_usage,visualize_codon_usage,optimize_codon_usage

print('\n')
print('Hi, I will optimize the codon usage of your gene.')
print('\n')
print('First, I need some input. Do NOT type any file extensions unless I ask you to. Also, no typos, please.')
print('\n')
gene = input('Enter name of gene to be optimized: ')
gene_OPT = gene + '_opt'
sequence = re.sub(r"[\n\t\s]*", "", input('Copy and paste the DNA sequence of your gene (only CDS): ').upper())

organism_ORIGIN = input('Enter name of original host organism: ')
CDS_source_nofile_ORIGIN = input("Enter name of its genbank file ('xxx'.gb) or csv file from Kazusa ('xxx'.csv): ")
CDS_type_ORIGIN = input("file extension: 'gb' or 'csv?'")
CDS_source_ORIGIN = CDS_source_nofile_ORIGIN + '.' + CDS_type_ORIGIN
save_table_ORIGIN = CDS_source_nofile_ORIGIN + '_codon_usage'
CDSdict_txt_ORIGIN = 'CDS_dict_' + save_table_ORIGIN + '.txt'
save_table_txt_ORIGIN = save_table_ORIGIN + '.txt'
save_table_csv_ORIGIN = save_table_ORIGIN + '.csv'
savefig_ORIGIN = 'codon_usage_' + gene + '_vs_' + organism_ORIGIN + '.pdf'

organism_DEST = input('Enter name of destination organism: ')
CDS_source_nofile_DEST = input("Enter name of its genbank file ('xxx'.gb) or csv file from Kazusa ('xxx'.csv): ")
CDS_type_DEST = input("file extension: 'gb' or 'csv?'")
CDS_source_DEST = CDS_source_nofile_DEST + '.' + CDS_type_DEST
save_table_DEST = CDS_source_nofile_DEST + '_codon_usage'
CDSdict_txt_DEST = 'CDS_dict_' + save_table_DEST + '.txt'
save_table_txt_DEST = save_table_DEST + '.txt'
save_table_csv_DEST = save_table_DEST + '.csv'
savefig_DEST = 'codon_usage_' + gene + '_vs_' + organism_DEST + '.pdf'

save_optimized_seq = gene_OPT + '_for_' + organism_DEST + '.txt'
savefig_OPT = 'codon_usage_' + gene_OPT + '_vs_' + organism_DEST + '.pdf'

print('I create codon usage tables for both organisms.')

get_codon_usage(CDS_source_nofile_ORIGIN,
                    CDS_type_ORIGIN,
                    CDS_source_ORIGIN,
                    CDSdict_txt_ORIGIN,
                    save_table_ORIGIN,
                    save_table_txt_ORIGIN,
                    save_table_csv_ORIGIN)
get_codon_usage(CDS_source_nofile_DEST,
                    CDS_type_DEST,
                    CDS_source_DEST,
                    CDSdict_txt_DEST,
                    save_table_DEST,
                    save_table_txt_DEST,
                    save_table_csv_DEST)

print("I create bar charts for you.")

visualize_codon_usage(organism_ORIGIN,
                      CDS_type_ORIGIN,
                      save_table_ORIGIN,
                      save_table_txt_ORIGIN,
                      gene,
                      sequence,
                      savefig_ORIGIN)
visualize_codon_usage(organism_DEST,
                      CDS_type_DEST,
                      save_table_DEST,
                      save_table_txt_DEST,
                      gene,
                      sequence,
                      savefig_DEST)
print('Now, you can take a look at the figures that I made for you.')
print("Do you still want to optimize " + gene + " for the usage of " + organism_DEST + "?")
optimization = input("Type 'Y' or 'N': ")

if optimization == 'Y' or optimization == 'y':
    print("OK, a '.txt' file containing the optimized sequence will be stored in your folder.")

    optimize_codon_usage(save_table_ORIGIN,
                         save_table_txt_ORIGIN,
                         save_table_DEST,
                         save_table_txt_DEST,
                         gene,
                         sequence,
                         organism_ORIGIN,
                         organism_DEST,
                         save_optimized_seq
                         )
    with open(save_optimized_seq, 'r') as opt_seq:
        sequence_OPT = opt_seq.read().replace('\n', '')
        opt_seq.close

    print('Codon usage optimized.')
    print('\n')
    print("Lastly, let's generate bar charts of the codon optimized gene.")

    visualize_codon_usage(organism_DEST,
                          CDS_type_DEST,
                          save_table_DEST,
                          save_table_txt_DEST,
                          gene_OPT,
                          sequence_OPT,
                          savefig_OPT)
        
    print('Bye.')
if optimization == 'N' or optimization == 'n':
    print('Ok, bye.')
