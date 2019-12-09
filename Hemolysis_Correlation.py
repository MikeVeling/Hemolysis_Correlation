###############################################################################
#                            Import dependencies                              #
###############################################################################
import codecs, csv, os
from scipy.stats import pearsonr as R_val
from scipy.stats import spearmanr as S_R_val
cwd=os.getcwd()
if '\\' in cwd:
    slash='\\'
elif '/' in cwd:
    slash='/'
cwd+=slash
###############################################################################
#                                 user var                                    #
###############################################################################
data_len_cutoff=20
prot_data_location=cwd+'Proteins.csv'
heme_data_location=cwd+'Hemolysis.csv'
output_location=cwd+'Outputs'+slash
###############################################################################
#                                Functions                                    #
###############################################################################
protein_dic={}
hemolysis_dic={}
def the_opener(file_path):
    if '.csv' in file_path:
        return_list=[]
        with codecs.open(file_path, 'rU',encoding='utf-8-sig') as csvfile:
            reader=csv.reader(csvfile)
            for row in reader:
                return_list.append(row)
    elif '.tsv' in file_path or '.txt' in file_path:
        return_list=[]
        fixed_file='\n'.join(open(file_path).read().split('\r'))
        for row in fixed_file.split('\n'):
            if row != '':
                return_list.append(row.split('\t'))
    return return_list
def mkdir_safe(path):
    try:
        os.mkdir(path)
    except:
        pass
def the_saver(output_path, list_of_lists):
    with open(output_path,'w', newline='') as csvfile:
        for row in list_of_lists:
            csv.writer(csvfile, dialect='excel').writerow(row)
def make_data_float(data):
    return_data=[]
    for row in data:
        if row != '':
            return_data.append(float(row))
        else:
            return_data.append('')
    return return_data
def simp_2_sets(data_1,data_2):
    assert len(data_1)==len(data_2)
    return_data_1=[]
    return_data_2=[]
    for data_ID in range(len(data_1)):
        if data_1[data_ID] != '' and data_2[data_ID] != '':
            return_data_1.append(data_1[data_ID])
            return_data_2.append(data_2[data_ID])
    return return_data_1, return_data_2

def get_output(data_1,data_2,output):
    data_1_float=make_data_float(data_1)
    data_2_float=make_data_float(data_2)
    simp_data_1,simp_data_2=simp_2_sets(data_1_float,data_2_float)
    data_count=len(simp_data_1)
    if data_count>data_len_cutoff:
        pearson_coefic=str(R_val(simp_data_1,simp_data_2)[0])
        spearmans_coefic=str(S_R_val(simp_data_1,simp_data_2)[0])
    else:
        pearson_coefic='N/A'
        spearmans_coefic='N/A'
    if output == 'Spearman\'s':
        return spearmans_coefic
    elif output == 'Pearson\'s':
        return pearson_coefic
    elif output == 'Data Count':
        return data_count
###############################################################################
#                                Load data                                    #
###############################################################################
protein_data=the_opener(prot_data_location)
hemolysis_data=the_opener(heme_data_location)
###############################################################################
#                              Analyze data                                   #
###############################################################################
sider_infos=[]
for row in protein_data[1:]:
    header_info='|||'.join(row[0:5])
    sider_infos.append(header_info)
    data=row[5:]
    protein_dic[header_info]=data

header_infos=[]
for row in hemolysis_data[1:]:
    header_info='|||'.join(row[0:3])
    header_infos.append(header_info)
    data=row[3:]
    hemolysis_dic[header_info]=data

mkdir_safe(output_location)
for output in ['Spearman\'s','Pearson\'s','Data Count']: 
    output_file=[]
    for ting in range(0,3):
        if ting != 2:
            next_line=['','','','','']
        else:
            next_line=['Protein Name',
                       'Protein Discription',
                       'FASTA header',
                       'Item \ ID',
                       'Heritability']
        for header in header_infos:
            next_line.append(header.split('|||')[ting])
        output_file.append(next_line)
    
    done_sets=[]
    for sider_info in sider_infos:
        next_line=sider_info.split('|||')
        for header_info in header_infos:
            working_set=set([sider_info,header_info])
            if working_set not in done_sets:
                done_sets.append(working_set)
                data_1=protein_dic[sider_info]
                data_2=hemolysis_dic[header_info]
                next_line+=[get_output(data_1,data_2,output)]
            else:
                next_line+=['']
        output_file.append(next_line)
    the_saver(output_location+output+'.csv',output_file)
###############################################################################
#                                    Fin                                      #
###############################################################################