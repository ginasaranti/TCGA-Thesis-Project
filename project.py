from os import listdir
from os.path import isfile, join
import scipy
from scipy.stats import ttest_ind
import seaborn as sns
import pandas as pd
import numpy as np
import itertools
import os
import statistics
import natsort
from matplotlib import pyplot as plt
from pprint import pprint

def prev():

    def correlations(df):
        subdf_list = separated_chr(df)


        total_corr_list = []
        for chromosome in subdf_list:
            corr_list = []
            for i in range(len(chromosome)):
                if i == 0:
                    corr1 = chromosome.iloc[i, 7:].to_list()
                    corr2 = chromosome.iloc[i + 1, 7:].to_list()
                    corr_list.append(scipy.stats.pearsonr(np.array(np.array(corr1)),np.array(np.array(corr2))))

                elif i == (len(chromosome) - 1) :
                    corr1 = chromosome.iloc[i, 7:].to_list()
                    corr2 = chromosome.iloc[i - 1, 7:].to_list()
                    corr_list.append(scipy.stats.pearsonr(np.array(np.array(corr1)),np.array(np.array((corr2)))))

                else:
                    corr = chromosome.iloc[i, 7:].to_list()
                    corr1 = chromosome.iloc[i + 1, 7:].to_list()
                    corr2 = chromosome.iloc[i - 1, 7:].to_list()
                    temp1 = scipy.stats.pearsonr(np.array(np.array(corr)),np.array(np.array(corr1)))
                    temp2 = scipy.stats.pearsonr(np.array(np.array(corr)),np.array(np.array(corr2)))
                    corr_list.append((temp1,temp2))
            total_corr_list.append(corr_list)

            for chromosome in subdf_list:
                chromosome.reset_index(drop= True, inplace = True)
        return total_corr_list

    def variance(df):
        subdf_list = separated_chr(df)

        total_variance_list = []
        for chromosome in subdf_list:
            variance_list = []
            for i in range(len(chromosome)):
                row = chromosome.iloc[i, 7:].to_list()
                variance_list.append(statistics.stdev(row))

            total_variance_list.append(variance_list)
        return total_variance_list



    #TCGA data is from 1222 patients for breast cancer only.

    onlyfiles = [f for f in listdir('/media/gina/9A53-2BCF/gina/all_files') if isfile(join('/media/gina/9A53-2BCF/gina/all_files', f))]
    os.chdir('/home/gina/GINA/all_files')

    names = ['']
    counts = ['']


    #Names list was created initially to check if the number of gene ids is the same in all files

    for file in onlyfiles:
        df = pd.read_csv(file, sep='\t', header = None, names = ["Gene_version", "Counts"])
        names.extend((df["Gene_version"]))

    names = names[1:]


    #All files contain 60483 genes so I proceed to the dataframe
    #Create TCGA dataframe

    for i in range(len(onlyfiles)):
        if i == 0 :
            tcga_df =  pd.read_csv(onlyfiles[i], sep='\t', header = None, names = ["Gene_version", "Counts"])
        else:
            counts = list(pd.read_csv(onlyfiles[i], sep='\t', header = None, names = ["Gene_version", "Counts"])["Counts"])
            tcga_df[i] = counts


    onlyfiles.insert(0, 'Gene_version')
    tcga_df.columns = onlyfiles

    #Bring only files to its previous state

    onlyfiles = [f for f in listdir('/media/gina/9A53-2BCF/gina/all_files') if isfile(join('/media/gina/9A53-2BCF/gina/all_files', f))]

    tcga_gene_version = tcga_df['Gene_version'].to_list()

    #FIND GENE COORDINATES

    #New dataframe from .gtf File

    gencode = pd.read_table("/media/gina/9A53-2BCF/gina/Homo_sapiens.GRCh38.104.gtf", comment="#", sep = "\t", names = ['seqname', 'source', 'feature', 'start' , 'end', 'score', 'strand', 'frame', 'attribute'])


    #Create a smaller dataframe with the columns of interest (chr, start-end sight, and attribute)


    gencode_genes = gencode[(gencode.feature == "gene")][['seqname', 'start', 'end', 'attribute']]

    #Resetting the index
    gencode_genes.reset_index(drop= True, inplace = True)

    gen_id_version_list = []
    gen_id_list = []
    gen_name_list = []
    gen_biotype_list = []

    #Taking certain infos from attribute column

    for i in range(len(gencode_genes)):
        g_id = gencode_genes['attribute'][i].split(';')[0].split(' ')[1].replace('"', '')
        version = gencode['attribute'][i].split(';')[1].split(' ')[2].replace('"', "")
        gen_id_list.append(g_id)
        gen_id_version_list.append(g_id+'.'+version)
        if 'gene_name' in gencode_genes['attribute'][i]:

            g_name = gencode_genes['attribute'][i].split(';')[2].split(' ')[2].replace('"', '')
            gen_name_list.append(g_name)

            gene_biotype = gencode_genes['attribute'][i].split(';')[4].split(' ')[2].replace('"', '')
            gen_biotype_list.append(gene_biotype)
        else:
            g_name = '-'
            gen_name_list.append(g_name)

            gen_biotype = gencode_genes['attribute'][i].split(';')[3].split(' ')[2].replace('"', '')
            gen_biotype_list.append(gene_biotype)

    #Adding everything to the GENCODE_GENES dataframe and removing the attribute column

    gencode_genes['gen_version'] = gen_id_version_list
    gencode_genes['gen_id'] = gen_id_list
    gencode_genes['gen_name'] = gen_name_list
    gencode_genes['gen_biotype'] = gen_biotype_list

    gencode_genes.drop(['attribute'], inplace = True, axis=1)

    #Merging the dataframes from TCGA and GENCODE, creating a new one that contains only the common genes

    tcga_id_list = []

    for i in range(len(tcga_gene_version)):
        tcga_gene_id = tcga_gene_version[i].split('.')[0]
        tcga_id_list.append(tcga_gene_id)

    tcga_df['Gene id'] = tcga_id_list

    common_df = gencode_genes.loc[gencode_genes['gen_id'].isin(tcga_id_list)]
    common_df = common_df.loc[common_df['gen_id'].isin(tcga_df['Gene id'])]
    common_df.reset_index(drop= True, inplace = True)

    #Sorting by chromosome
    common_df = pd.read_csv('/media/gina/9A53-2BCF/gina/aimiliosmegas.csv', sep ='\t', header = None, names = ["Chromosome", "Start", "End", "Gene Version", "Gene Id", "Gene Name", "Gene Biotype"])
    #Creating the final df by combining common_df and tcga_df (with the selected genes)

    final_gene_list = []
    final_gene_list = common_df['Gene Id'].to_list()
    final_df = tcga_df.loc[tcga_df['Gene id'].isin(final_gene_list)]
    final_df.reset_index(drop= True, inplace = True)

    final_df = pd.concat([common_df, final_df], axis=1)
    final_df.drop('Gene_version', axis=1, inplace=True)
    final_df = final_df.drop(final_df.columns[1229], axis=1)
    final_df.reset_index(drop= True, inplace = True)
    final_df = final_df.iloc[natsort.index_humansorted(final_df.Chromosome)]
    #final_df.to_csv('opos thes pes to.csv', sep = '\t')
    #Creating a list with separated chromosomes

    subdf_list = separated_chr(final_df)


    #Correlations between genes by chromosome

    total_corr_list = correlations(final_df)


    #List with the averages of the correlations separated by chromosome
    total_average_cor_list = []
    for chromosome_corr in total_corr_list:
        average_cor_list = []
        for i in range(len(chromosome_corr)):
            if i == 0:
                average = chromosome_corr[i][0]
                average_cor_list.append(average)

            elif i == (len(chromosome_corr) - 1):
                average = chromosome_corr[i][0]
                average_cor_list.append(average)

            else:
                average = (chromosome_corr[i][0][0] + chromosome_corr[i][1][0]) / 2
                average_cor_list.append(average)

        total_average_cor_list.append(average_cor_list)

    #Matching every patient with the disease type he diagnosed


    sample_df = pd.read_csv(r'/media/gina/9A53-2BCF/gina/final_sample_df.csv', sep = ',')
    sample_df.drop_duplicates(subset=['Sample ID'])
    sample_df.reset_index(drop= True, inplace = True)

    clinical_df = pd.read_excel(r'/media/gina/9A53-2BCF/gina/clinical.xlsx')
    clinical_df = clinical_df.drop_duplicates(subset=['case_submitter_id', 'primary_diagnosis'])
    clinical_df = clinical_df.rename(columns={'case_submitter_id': 'Case ID'})

    clinical_df.reset_index(drop= True, inplace = True)

    diagnosis_df = pd.concat([clinical_df['Case ID'], clinical_df['primary_diagnosis']], axis=1)
    diagnosis_df.reset_index(drop= True, inplace = True)


    #In sample_df['Sample Type'] could be Primary Tumor or Metastatic Cancer or Solid Tissue
    #We want patients with primary tumor as a primary diagnosis

    for i in range(len(sample_df['File Name'])):
       sample_df['File Name'][i] = sample_df['File Name'][i].replace('.gz', '')

    sample_df = sample_df.loc[sample_df['File Name'].isin(onlyfiles)]


    for i in range(len(sample_df['Sample Type'])):
        if sample_df['Sample Type'][i] != 'Primary Tumor':
            sample_df.drop([i],inplace = True)

    sample_df.reset_index(drop = True, inplace = True)

    temp_list = []
    for i in range(len(sample_df['Case ID'])):
        if sample_df['Case ID'][i] not in temp_list:
            temp_list.append(sample_df['Case ID'][i])
        else:
            sample_df.drop([i], inplace = True)

    sample_df.reset_index(drop = True, inplace = True)

    sample_df = sample_df.merge(diagnosis_df, on='Case ID')
    sample_df.reset_index(drop = True, inplace = True)


    #FIND ALL THE SUBTYPES

    #In clinical_df primary diagnosis is the subtype of the cancer (infiltrating ducy carcinoma, etc.)
    diagnosis_list = []

    for x in clinical_df['primary_diagnosis']:
        if x not in diagnosis_list:
            diagnosis_list.append(x)

    subtypes_counter = {subtypes: 0 for subtypes in diagnosis_list}

    for i in range(len(sample_df['primary_diagnosis'])):
        for subtype in diagnosis_list:
            if sample_df['primary_diagnosis'][i] == subtype:
                subtypes_counter[subtype] += 1


    # Match each patient with his primary diagnosis
    sample_df.to_csv('final_sample_df.csv', index=False)
    final_df = pd.read_csv(r'/media/gina/9A53-2BCF/gina/final_df.csv')
    infiltrating_duct_carcinoma_df =  pd.read_csv(r'/media/gina/9A53-2BCF/gina/final_sample_df.csv')
    lobular_carcinoma_df = pd.read_csv(r'/media/gina/9A53-2BCF/gina/final_sample_df.csv')

    for i in range(len(sample_df)):
        if infiltrating_duct_carcinoma_df['primary_diagnosis'][i] != 'Infiltrating duct carcinoma, NOS':
            infiltrating_duct_carcinoma_df.drop([i],inplace = True)


    for i in range(len(sample_df)):
        if lobular_carcinoma_df['primary_diagnosis'][i] != 'Lobular carcinoma, NOS':
            lobular_carcinoma_df.drop([i],inplace = True)

    infiltrating_list = infiltrating_duct_carcinoma_df['File Name'].tolist()
    lobular_list = lobular_carcinoma_df['File Name'].tolist()

    lobular_carcinoma_counts_df = final_df.iloc[:,:7 ].copy()
    for name in lobular_list:
            lobular_carcinoma_counts_df[name] = final_df[name]

    infiltrating_duct_carcinoma_counts_df = final_df.iloc[:,:7 ].copy()
    for name in infiltrating_list:
            infiltrating_duct_carcinoma_counts_df[name] = final_df[name]


    #Correlations by subtype

    lobular_corr = correlations(lobular_carcinoma_counts_df)
    infiltrating_corr = correlations(infiltrating_duct_carcinoma_counts_df)

    lob_average_cor_list = []
    for chromosome_corr in lobular_corr:
        average_cor_list = []
        for i in range(len(chromosome_corr)):
            if i == 0:
                average = chromosome_corr[i][0]
                average_cor_list.append(average)

            elif i == (len(chromosome_corr) - 1):
                average = chromosome_corr[i][0]
                average_cor_list.append(average)

            else:
                average = (chromosome_corr[i][0][0] + chromosome_corr[i][1][0]) / 2
                average_cor_list.append(average)

        lob_average_cor_list.append(average_cor_list)

    inf_average_cor_list = []
    for chromosome_corr in infiltrating_corr:
        average_cor_list = []
        for i in range(len(chromosome_corr)):
            if i == 0:
                average = chromosome_corr[i][0]
                average_cor_list.append(average)

            elif i == (len(chromosome_corr) - 1):
                average = chromosome_corr[i][0]
                average_cor_list.append(average)

            else:
                average = (chromosome_corr[i][0][0] + chromosome_corr[i][1][0]) / 2
                average_cor_list.append(average)

        inf_average_cor_list.append(average_cor_list)



    lob_subdf = separated_chr(lobular_carcinoma_counts_df)
    inf_subdf = separated_chr(infiltrating_duct_carcinoma_counts_df)

    lobular_carcinoma_counts_df.to_csv('lobular', sep = '\t')
    infiltrating_duct_carcinoma_counts_df.to_csv('infiltrating', sep='\t')
    #Variance
    lobular_variance = variance(lobular_carcinoma_counts_df)
    infiltrating_variance = variance(infiltrating_duct_carcinoma_counts_df)

    #Boxplots for lobular correlations
    data = pd.DataFrame((lob_average_cor_list[0]), columns = ['Chr1'])
    data['Chr1'] = data['Chr1'].fillna(0)
    for i in range(len(lob_average_cor_list)):
        if i == 0 :
            continue
        else:

            data['Chr' + str(i + 1)] = pd.Series(lob_average_cor_list[i])
            data['Chr' + str(i + 1)] = data['Chr' + str(i + 1)].fillna(0)

    data = data.rename(columns={'Chr23': 'ChrX', 'Chr24': 'ChrY', 'Chr25': 'MT'})

    palette = ['#FF2709', '#09FF10', '#0030D7', '#FA70B5']
    data_melted = pd.melt(data)
    data_melted = data_melted.rename(columns={'variable': 'Chromosomes', 'value': 'Correlation'})

    sns.set_style("whitegrid")
    sns.despine(left=True)
    sns.set(rc={'figure.figsize':(15.7,8.27)})
    sns.boxplot(x='Chromosomes', y='Correlation', data=data_melted,linewidth=1, showfliers = False).set_title('Lobular Carcinoma Correlations')

    #Boxplot for infiltrating correlations

    data = pd.DataFrame((inf_average_cor_list[0]), columns = ['Chr1'])
    data['Chr1'] = data['Chr1'].fillna(0)
    for i in range(len(inf_average_cor_list)):
        if i == 0 :
            continue
        else:

            data['Chr' + str(i + 1)] = pd.Series(inf_average_cor_list[i])
            data['Chr' + str(i + 1)] = data['Chr' + str(i + 1)].fillna(0)

    data = data.rename(columns={'Chr23': 'ChrX', 'Chr24': 'ChrY', 'Chr25': 'MT'})

    palette = ['#FF2709', '#09FF10', '#0030D7', '#FA70B5']
    data_melted = pd.melt(data)
    data_melted = data_melted.rename(columns={'variable': 'Chromosomes', 'value': 'Correlation'})

    sns.set_style("whitegrid")
    sns.despine(left=True)
    sns.set(rc={'figure.figsize':(15.7,8.27)})
    sns.boxplot(x='Chromosomes', y='Correlation', data=data_melted,linewidth=1, showfliers = False).set_title('Infiltrating Duct Carcinoma Correlations')

    # BOXPLOT FOR LOBULAR VARIANCE 
    data = pd.DataFrame((lobular_variance[0]), columns = ['Chr1'])
    data['Chr1'] = data['Chr1'].fillna(0)
    for i in range(len(lobular_variance)):
        if i == 0 :
            continue
        else:

            data['Chr' + str(i + 1)] = pd.Series(lobular_variance[i])
            data['Chr' + str(i + 1)] = data['Chr' + str(i + 1)].fillna(0)

    data = data.rename(columns={'Chr23': 'ChrX', 'Chr24': 'ChrY', 'Chr25': 'MT'})

    palette = ['#FF2709', '#09FF10', '#0030D7', '#FA70B5']
    data_melted = pd.melt(data)
    data_melted = data_melted.rename(columns={'variable': 'Chromosomes', 'value': 'Variance'})

    sns.set_style("whitegrid")
    sns.despine(left=True)
    sns.set(rc={'figure.figsize':(15.7,8.27)})
    sns.boxplot(x='Chromosomes', y='Variance', data=data_melted,linewidth=1, showfliers = False).set_title('Lobular Carcinoma Variance')
    
    #BOXPLOT FOR INFILTRATING VARIANCE
                                               
    data = pd.DataFrame((infiltrating_variance[0]), columns = ['Chr1'])
    data['Chr1'] = data['Chr1'].fillna(0)
    for i in range(len(infiltrating_variance)):
        if i == 0 :
            continue
        else:

            data['Chr' + str(i + 1)] = pd.Series(infiltrating_variance[i])
            data['Chr' + str(i + 1)] = data['Chr' + str(i + 1)].fillna(0)

    data = data.rename(columns={'Chr23': 'ChrX', 'Chr24': 'ChrY', 'Chr25': 'MT'})

    palette = ['#FF2709', '#09FF10', '#0030D7', '#FA70B5']
    data_melted = pd.melt(data)
    data_melted = data_melted.rename(columns={'variable': 'Chromosomes', 'value': 'Variance'})

    sns.set_style("whitegrid")
    sns.despine(left=True)
    sns.set(rc={'figure.figsize':(15.7,8.27)})
    sns.boxplot(x='Chromosomes', y='Variance', data=data_melted,linewidth=1, showfliers = False).set_title('Infiltrating Duct Carcinoma Variance')
                                                           
    
    #lineplots for lobular correlation


    for i in range(len(lob_average_cor_list)):

        data = pd.DataFrame(lob_subdf[i]["Start"])
        data['Correlation'] = lob_average_cor_list[i]
        figure = plt.figure()
        if i <= 21:
            sns.lineplot(x='Start', y='Correlation', data=data,linewidth=1).set(title=('Chromosome ' + str(i + 1)))
        elif i == 22:
            sns.lineplot(x='Start', y='Correlation', data=data,linewidth=1).set(title=('Chromosome X'))
        elif i == 23:
            sns.lineplot(x='Start', y='Correlation', data=data,linewidth=1).set(title=('Chromosome Y'))
        elif i == 24:
            sns.lineplot(x='Start', y='Correlation', data=data,linewidth=1).set(title=('MT'))

        plt.savefig('CHR' + str(i +1) + ".png")



    for i in range(len(inf_average_cor_list)):

        data = pd.DataFrame(inf_subdf[i]["Start"])
        data['Correlation'] = inf_average_cor_list[i]
        figure = plt.figure()
        if i <= 21:
            sns.lineplot(x='Start', y='Correlation', data=data,linewidth=1).set(title=('Chromosome ' + str(i + 1)))
        elif i == 22:
            sns.lineplot(x='Start', y='Correlation', data=data,linewidth=1).set(title=('Chromosome X'))
        elif i == 23:
            sns.lineplot(x='Start', y='Correlation', data=data,linewidth=1).set(title=('Chromosome Y'))
        elif i == 24:
            sns.lineplot(x='Start', y='Correlation', data=data,linewidth=1).set(title=('MT'))

        plt.savefig('CHR' + str(i +1) + ".png")



    def separated_chr(df):
        subdf_list = []
        for i in chrom_list:
            subdf_list.append(df[df['Chromosome'] == i])
        return subdf_list

#gene_names_per_chr = []
#for number in chrom_list:
    #gene_names = []
    #for i in range(len(lobular_counts_df)):
        #if lobular_counts_df['Chromosome'][i] == number:
            #gene_names.append(lobular_counts_df['Gene Name'][i])
    #gene_names_per_chr.append(gene_names)

########################################################################################################################################

def separate_chroms(df):
    chrom_list = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y','MT']

    chr_dict = {}
    for c in chrom_list:
        chr_dict[c] = df.loc[df.Chromosome == c]

    return chr_dict


def segment(df):
    window = 10_000
    curr_bound = 2 * window
    seg_num = 0

    chr_end = df.End.max()
    patient_count = len(df.columns) - 7

    print( "===========================================")
    print(f"[*] Analyzing chromosome: {df.Chromosome.iloc[0]}")
    print(f"[*] Chromosome end:       {chr_end}")
    print(f"[*] Patients:             {patient_count}")

    segments = []
    for i, row in df.iterrows():

        if row.Start >= curr_bound:
            segments.append([0 for _ in df.iloc[:, 7:]])
            seg_num += 1
            curr_bound += window

        seg_i = seg_num
        if (row.End > curr_bound) and (row.Start + row.End > 2 * curr_bound):
            seg_i = seg_num + 1

        diff = seg_i - len(segments)
        if diff >= 0:
            for _ in range(diff + 1):
                segments.append([0 for _ in df.iloc[:, 7:]])


        for j, val in enumerate(row.iloc[7:]):
            segments[seg_i][j] += val

    return segments

########################################################################################################################################

cancers = [None, None]
cancers[0] = pd.read_csv(os.getcwd() + '/infiltrating.csv', dtype={'Chromosome': 'str'}, sep = '\t')
cancers[1] = pd.read_csv(os.getcwd() + '/lobular.csv',      dtype={'Chromosome': 'str'}, sep = '\t')

gp_chrom = []
for df in cancers:
    #df.Chromosome = df.Chromosome.astype(str)
    gp_chrom.append(separate_chroms(df))

seg_counts_pchr_pcan = []
for d in gp_chrom:
    seg_counts_pchr_pcan.append({key: segment(value) for key, value in d.items()})
    print( "===========================================")

ttests = [ttest_ind(c1, c2, axis=1) for c1, c2 in zip(list(seg_counts_pchr_pcan[0].values()), list(seg_counts_pchr_pcan[1].values()))]


for i in enumerate(ttests):
    print(type(i))

sig_pv = []
unsig_pv = []
for ch_i, ch in enumerate(ttests):
    for pv_i, pv in enumerate(ch[1]):
        if pv <= 0.05:
            sig_pv.append((ch_i, pv_i, ch[0][pv_i]))
        else:
            unsig_pv.append((ch_i, pv_i, ch[0][pv_i]))

   
   
   data = pd.DataFrame((sig_pv[0][2]), columns = [sig_pv[0][1]])
    data['ch'] = data['ch'].fillna(0)
    data['pvalue'] = unsig_pv[0][2]
    
for i in enumerate(sig_pv[0]):
    
    data = data.rename(columns={'variable': 'Segment', 'value': 'P-value'})
    
    sns.set_style("whitegrid")
    sns.despine(left=True)

    sns.set(rc={'figure.figsize':(15.7,8.27)})
    palette = ['#FF2709', '#09FF10', '#0030D7', '#FA70B5']
    sns.scatterplot(data = data , x = "Segmnet", y = "P-value")

 plt.show()

data = pd.DataFrame(sig_pv[0][2])

for i in enumerate(sig_pv):
    if sig_pv[0] = 0:
        data = pd.DataFrame(sig_pv[i][])
        data['Correlation'] = lob_average_cor_list[i]





