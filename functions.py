import json
from decimal import Decimal
from Bio import SeqIO
import csv
from decimal import Decimal
import matplotlib.pyplot as plt
import numpy as np

synonymous_dict = {
'C': ['TGT', 'TGC'],
'D': ['GAT', 'GAC'], 
'S': ['TCT', 'TCG', 'TCA', 'TCC', 'AGC', 'AGT'], 
'Q': ['CAA', 'CAG'], 
'M': ['ATG'], 
'N': ['AAC', 'AAT'], 
'P': ['CCT', 'CCG', 'CCA', 'CCC'], 
'K': ['AAG', 'AAA'], 
'*': ['TAG', 'TGA', 'TAA'], 
'T': ['ACC', 'ACA', 'ACG', 'ACT'], 
'F': ['TTT', 'TTC'], 
'A': ['GCA', 'GCC', 'GCG', 'GCT'], 
'G': ['GGT', 'GGG', 'GGA', 'GGC'], 
'I': ['ATC', 'ATA', 'ATT'], 
'L': ['TTA', 'TTG', 'CTC', 'CTT', 'CTG', 'CTA'], 
'H': ['CAT', 'CAC'], 
'R': ['CGA', 'CGC', 'CGG', 'CGT', 'AGG', 'AGA'], 
'W': ['TGG'], 
'V': ['GTA', 'GTC', 'GTG', 'GTT'], 
'E': ['GAG', 'GAA'], 
'Y': ['TAT', 'TAC'] 
}

aa_dict = {'Cys': 'C',
           'Asp': 'D',
           'Ser': 'S',
           'Gln': 'Q',
           'Met': 'M',
           'Asn': 'N',
           'Pro': 'P',
           'Lys': 'K',
           'End': '*',
           'Thr': 'T',
           'Phe': 'F',
           'Ala': 'A',
           'Gly': 'G',
           'Ile': 'I',
           'Leu': 'L',
           'His': 'H',
           'Arg': 'R',
           'Trp': 'W',
           'Val': 'V',
           'Glu': 'E',
           'Tyr': 'Y'
           }

def get_codon_usage(CDS_source_nofile,CDS_type,CDS_source,CDSdict_txt,save_table,save_table_txt,save_table_csv):
    CDSnum = []
    CDSs = []
    codons = []
    total_codons = []
    cutoff = 300
    codon_dict = {}
    aa_codon = 0
    reladap_codon = 0
    aa_list = []
    rowlist = []
    if CDS_type == 'gb':
        for rec in SeqIO.parse(CDS_source, CDS_type):
            if rec.features:
                for feature in rec.features:
                    if feature.type == "CDS":
                        if len(feature) >= cutoff:
                            CDSnum.append(len(CDSs) + 1)
                            CDSs.append(str(feature.location.extract(rec).seq))
                for CDSno in range(0, len(CDSs), 1):
                    for i in range(0, len(CDSs[CDSno]), 3):
                        codons.append(CDSs[CDSno][i:i+3])
                        total_codons.append(CDSs[CDSno][i:i+3])
                    CDSs[CDSno] = codons
                    codons = []

        num_total_codons = len(total_codons)
        CDSlist = zip(CDSnum,CDSs)
        CDSdict =dict(CDSlist)

        for aa,codons in synonymous_dict.items():
            for codon_index,codon in enumerate(codons):
                aa_codon += total_codons.count(codon)
                if total_codons.count(codon) >= reladap_codon:
                    reladap_codon = total_codons.count(codon)
                    
            for codon_index,codon in enumerate(codons):
                codon_dict[codon] = [total_codons.count(codon),
                                     round(total_codons.count(codon)/num_total_codons*1000,2),
                                     round(total_codons.count(codon)/aa_codon*100,2),
                                     round(total_codons.count(codon)/reladap_codon*100,0)]
            synonymous_dict[aa] = codon_dict
            codon_dict = {}
            aa_codon = 0
            reladap_codon = 0

    if CDS_type == 'csv':
        with open(CDS_source,'r', encoding = 'utf-8') as codon_usage:
            for rownum,row in enumerate(codon_usage):
                if row != ' \n' and rownum != 0:
                    rowa = row.split()
                    rowlist.append(rowa[1:4])
                if row[0:3] in aa_dict.keys():
                    aa_list.append(aa_dict[row[0:3]])
            codon_usage.close
        for row in rowlist:
            row = row.append('0')
        for row in rowlist:
            row = row.append('0')
        for aa,codons in synonymous_dict.items():
            for codonnum,codon in enumerate(codons):
                for rownum, row in enumerate(rowlist):
                    if row[0] == codon:
                        aa_codon += int(float(row[1]))
                        if int(float(row[1])) >= reladap_codon:
                            reladap_codon = int(float(row[1]))
            for codonnum,codon in enumerate(codons):
                for rownum, row in enumerate(rowlist):
                    if row[0] == codon:
                        row[3] = round(int(float(row[1]))/aa_codon*100,2)
                        row[4] = round(int(float(row[1]))/reladap_codon*100,0)
                        codon_dict[codon] = row[1:5]
            aa_codon = 0
            reladap_codon = 0
            synonymous_dict[aa] = codon_dict
            codon_dict = {}

    with open(save_table_txt, 'w') as codon_usage:
        codon_usage.write(json.dumps(synonymous_dict))
        codon_usage.close
    with open(save_table_csv, 'w') as codon_usage:
        fieldnames = ['amino acid', 'codon', 'number of codons', 'frequency per thousand', 'fraction','relative adaptiveness']
        w = csv.writer(codon_usage)
        w.writerow(fieldnames)
        for aa, codons in synonymous_dict.items():
            for codon,valuelist in codons.items():
                w.writerow([aa,codon,valuelist[0], valuelist[1],valuelist[2], valuelist[3]])
        codon_usage.close
    print(save_table_txt + ' saved.')
    print(save_table_csv + ' saved.')
    return CDS_source_nofile,CDS_type,CDS_source,CDSdict_txt,save_table,save_table_txt

def visualize_codon_usage(organism,CDS_type,save_table,save_table_txt,gene,sequence,savefig):
    space = ' '
    seqlist = []
    num_codon_list = []
    freq_adap_dict = {}
    freq_list = []
    adap_list = []
    freq_colorlist = []
    adap_colorlist = []
    codon_dict = {}

    for codon_index in range(0, len(sequence), 3):
        seqlist.append(sequence[codon_index:codon_index+3])

    with open(save_table_txt, 'r') as codon_usage:
        synonymous_dict = eval(codon_usage.read())
        for aa, codons in synonymous_dict.items():
            for codon, valuelist in codons.items():
                codon_dict[codon] = valuelist
        for codon_index,codon in enumerate(seqlist):
            for triplet in codon_dict.keys():
                if codon == triplet:
                    freq_adap_dict[codon_index+1] = [codon, codon_dict[codon][2], codon_dict[codon][3]]
                    freq_list.append(round(codon_dict[codon][2]))
                    if codon_dict[codon][2] >= 20:
                        freq_colorlist.append('black')
                    elif codon_dict[codon][2] >= 10:
                        freq_colorlist.append('grey')
                    else:
                        freq_colorlist.append('red')
            for triplet in codon_dict.keys():
                if codon == triplet:
                    adap_list.append(round(codon_dict[codon][3]))
                    if codon_dict[codon][3] >= 20:
                        adap_colorlist.append('black')
                    elif codon_dict[codon][3] >= 10:
                        adap_colorlist.append('grey')
                    else:
                        adap_colorlist.append('red')
            for aa, codons in synonymous_dict.items():
                for a, b in codons.items():
                    if codon == a:
                        seqlist[codon_index] = aa+'-'+a
        codon_usage.close

    N = len(seqlist)
    x = range(N)
    for number,codon in enumerate(seqlist):
        num_codon_list.append(str(number+1) + (4-len(str(number+1)))*2*space + codon)

    fig, (ax1,ax2) = plt.subplots(2,1,figsize=(len(seqlist)/5,15))
    plt.subplots_adjust(hspace=0.35)

    ax1.axis([-1, N, 0, 107])
    ax1.bar(x, freq_list, align='center', color=freq_colorlist)
    ax1.set_xticks(range(0,len(seqlist)))
    ax1.get_xaxis().set_tick_params(direction='out', pad=60)
    ax1.set_xticklabels(num_codon_list, rotation='vertical', size='9', verticalalignment = 'bottom')
    ax1.set_title('frequency', fontsize='40', color='blue', loc='left')
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)

    ax2.axis([-1, N, 0, 107])
    ax2.bar(x, adap_list,align='center', color=adap_colorlist)
    ax2.set_xticks(range(0,len(seqlist)))
    ax2.get_xaxis().set_tick_params(direction='out', pad=60)
    ax2.set_xticklabels(num_codon_list, rotation='vertical', size='9', verticalalignment = 'bottom')
    ax2.set_title('relative adaptiveness', fontsize='40', color='blue', loc='left')
    ax2.spines['right'].set_visible(False)
    ax2.spines['top'].set_visible(False)

    plt.text(0.75,1.04, '<10 in red', color='red', fontsize='15', transform=ax1.transAxes)
    plt.text(0.75,1, '<20 in grey', color='grey', fontsize='15', transform=ax1.transAxes)

    plt.text(0.25,1.04, str('gene: ') + gene, fontsize='15', transform=ax1.transAxes)
    plt.text(0.25,1, str('vs organism: ') + organism, fontsize='15', transform=ax1.transAxes)

    for a,b in zip(range(0,len(seqlist)),freq_list):
        ax1.text(a,b+4, str(b), horizontalalignment= 'center', verticalalignment = 'top', color=freq_colorlist[a], fontsize='7', rotation='vertical')

    for a,b in zip(range(0,len(seqlist)),adap_list):
        ax2.text(a,b+4, str(b), horizontalalignment= 'center', verticalalignment = 'top', color=adap_colorlist[a], fontsize='7', rotation='vertical')    

    plt.savefig(savefig, bbox_inches='tight')
    print('figure saved as ' + savefig)
    return organism,CDS_type,save_table,save_table_txt,gene,sequence,savefig

def optimize_codon_usage(save_table_ORIGIN,
                         save_table_txt_ORIGIN,
                         save_table_DEST,
                         save_table_txt_DEST,
                         gene,
                         sequence,
                         organism_ORIGIN,
                         organism_DEST,
                         save_optimized_seq):
    seqlist = []
    optimized_seqlist = []
    freq_list = []
    reladap_list = []
    short_freq_list = []
    long_freq_list = []
    best_freq_list = []

    for codon_index in range(0, len(sequence), 3):
        seqlist.append(sequence[codon_index:codon_index+3])

    with open(save_table_txt_ORIGIN, 'r') as origin_codon_usage:
        origin_synonymous_dict = eval(origin_codon_usage.read())
        for codon in seqlist:
            for origin_aa, origin_codons in origin_synonymous_dict.items():
                for origin_codon, origin_valuelist in origin_codons.items():
                    if codon == origin_codon:
                        freq_list.append(origin_valuelist[2])
                        reladap_list.append(origin_valuelist[3])
        origin_codon_usage.close()

    with open(save_table_txt_DEST, 'r') as dest_codon_usage:
        dest_synonymous_dict = eval(dest_codon_usage.read())
        for codonnum,codon in enumerate(seqlist):
            for dest_aa, dest_codons in dest_synonymous_dict.items():
                if codon in dest_codons.keys():
                    for codnum, cod in dest_codons.items():
                        short_freq_list.append(cod[2])
                    short_freq_array =np.array(short_freq_list)
                    long_freq_list.append(short_freq_list)
                    short_freq_list = []
        dest_codon_usage.close()

    for freqnum,freq in enumerate(freq_list):
        if freqnum == 0:
            best_freq_list.append('start')
            optimized_seqlist.append('ATG')
        elif reladap_list[freqnum] == 100.0:
            freq_array = np.array(long_freq_list[freqnum])
            max_freq = np.amax(freq_array)
            best_freq_list.append(max_freq)
        else:
            best_freqindex = (np.abs(np.array(long_freq_list[freqnum])-freq)).argmin()
            best_freq_list.append(long_freq_list[freqnum][best_freqindex])

    with open(save_table_txt_DEST, 'r') as dest_codon_usage:
        dest_synonymous_dict = eval(dest_codon_usage.read())
        for freqnum,freq in enumerate(best_freq_list):
            for dest_aa, dest_codons in dest_synonymous_dict.items():
                for codnum, cod in dest_codons.items():
                    if len(optimized_seqlist) < (freqnum+1):
                        if seqlist[freqnum] in dest_codons.keys() and freq in cod:
                            optimized_seqlist.append(codnum)
    optimized_seq = ''.join(map(str, optimized_seqlist))
    
    with open(save_optimized_seq, 'w') as opt_seq:
        opt_seq.write(optimized_seq)
        opt_seq.close
    print(save_optimized_seq + ' saved')
    return save_table_ORIGIN,save_table_txt_ORIGIN,save_table_DEST,save_table_txt_DEST,gene,sequence,save_optimized_seq
