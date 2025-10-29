import Bio
import requests
import matplotlib.pyplot as plt
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import numpy as np
import pandas as pd
from pandas import errors as pd_errors
from pandas import ExcelWriter
from pandas import ExcelFile

import urllib.parse
import urllib.request
import pybiomart

import io
from io import StringIO
import base64

import seaborn as sns
from bioservices import UniProt
from Bio import SeqIO
from pybiomart import Dataset
from pybiomart import Server

import xlsxwriter
import csv

requests.__version__

BASE = 'http://www.uniprot.org'
KB_ENDPOINT = '/uniprot/'
TOOL_ENDPOINT = '/up loadlists/'
aminoAA=[]
aminoAA_count=[]
aminoAAperc=[]
aminoAA_countperc=[]

labels = []
values = []
colors = []

#lists for chromosome count
AA=[]
AA_count_list=[]
df_walk = []
count_walk = []
walk_result = []

def write_excel(df):
        output = io.BytesIO()
    # Use the BytesIO object as the filehandle
        writer = pd.ExcelWriter(output, engine='xlsxwriter')
    # Write the data frame to the BytesIO object and save it
        df.to_excel(writer, sheet_name='Sheet1')
        writer.save()
        excel_data = output.getvalue()
        b64 = base64.b64encode(excel_data)
        AAcountxls = b64.decode()
        return AAcountxls

def mergeDict(dict1, dict2):
    dict3 = dict([(k,dict2[k]) for k in dict1])         
    return dict3

@app.route('/bar', methods=['post'])
def bar():
 
    AA=[]
    aminoAA=[]
     
    countAA=[]
    aminoAA_count=[]
    
    excel_file='protein.xlsx'
    sheet_name='sheet1'
    #writer=pd.ExcelWriter(excel_file, engine='xlsxwriter')    
    user_input=request.form['input']  
    output_excel = io.BytesIO()
    query = 'accession:'+user_input 
    #columnlist = "id, entry name,length,mass, families, citationmapping, citation, comment(DISEASE), comment(TEMPERATURE DEPENDENCE), comment(PH DEPENDENCE), comment(FUNCTION), comment(TISSUE SPECIFICITY), genes(PREFERRED), comment(SUBUNIT), comment(SUBCELLULAR LOCATION)"
    columnlist = "accession, id, length, mass, protein_families, cc_disease, temp_dependence, ph_dependence, cc_function, cc_subunit, cc_subcellular_location, sequence" 
    service = UniProt()

    #query1 = 'proteomecomponent:"chromosome"'+chromosome_number+' AND organism:"Homo sapiens (Human) [9606]" AND proteome:up000005640 AND reviewed:yes'
    #columnlist_query = "id,entry name,length"
    #result_uniprot_query = service.search(user_input, frmt="tab", columns=columnlist_query)

    #df_queryAC = service.get_df(user_input)
    try:
        new_prot_data = service.search(query, frmt="tsv", columns=columnlist)

        df_queryAC = pd.read_table(io.StringIO(new_prot_data))

        sequence=df_queryAC.iloc[0]['Sequence']

        aminoAA, aminoAA_count, aminoAA_count_rel, AA, RAPP = prot_analysis_new(sequence)

        countAA=ProteinAnalysis(sequence)
        count_amino=countAA.count_amino_acids()

        #df_AAcount = pd.DataFrame(aminoAA)
        df_prot=pd.read_csv(StringIO(new_prot_data), sep='\t')
        #csv_amino_count=pd.read_csv(StringIO(df_count), sep='\t')
 
        prot_xls=df_prot.to_excel(output_excel, encoding='utf-8', index=False, header=True, sheet_name=sheet_name)
        #writer.save()
        #writer.close()
        output_excel.seek(0)
        #b64=base64.b64encode(output_excel.read()).decode()
        #linko= f'<a href="data:application/vnd.openxmlformats-officedocument.spreadsheetml.sheet;base64,{b64}" download="protein.xlsx">Download excel file</a>'

        return render_template('bar.html', csv_data=new_prot_data, df_prot=df_prot, df_queryAC=df_queryAC, title='Hystogram', prot_xls_get=prot_xls, RAPP=RAPP, max=17000, AC=user_input, barlabels=aminoAA, barvalues=aminoAA_count)
    except (KeyError, pd_errors.EmptyDataError):
        error='invalid AC'
        return render_template('protein_home.html', error=error)          

@app.route('/bar2', methods=['post'])
def bar2():
    AA=[]
    aminoAA=[]
    aminoAA1=[]
    aminoAA2=[]
    countAA=[]

    aminoAA_count1=[]
    aminoAA_count2=[]
    user_input1=request.form['input1']
    user_input2=request.form['input2']
    query1 = 'accession:'+user_input1 
    query2 = 'accession:'+user_input2 
    #aminoAA_c=[]    
    #columnlist = "id, entry name,length,mass, families, citationmapping, citation, comment(DISEASE), comment(TEMPERATURE DEPENDENCE), comment(PH DEPENDENCE), comment(FUNCTION), comment(TISSUE SPECIFICITY), genes(PREFERRED), comment(SUBUNIT), comment(SUBCELLULAR LOCATION)"
    columnlist = "accession, id, length, mass, protein_families, cc_disease, temp_dependence, ph_dependence, cc_function, cc_subunit, cc_subcellular_location, sequence" 
    service = UniProt()

    #df_queryAC1 = service.get_df(user_input1)
    #df_queryAC2 = service.get_df(user_input2)
    try:
        new_prot_data1 = service.search(query1, frmt="tsv", columns=columnlist)
        new_prot_data2 = service.search(query2, frmt="tsv", columns=columnlist)
        df_queryAC1 = pd.read_table(io.StringIO(new_prot_data1))
        df_queryAC2 = pd.read_table(io.StringIO(new_prot_data2))
        sequence1=df_queryAC1.iloc[0]['Sequence']
        sequence2=df_queryAC2.iloc[0]['Sequence']
        aminoAA1, aminoAA_count1, aminoAA_count_rel, AA1, RAPP = prot_analysis_new(sequence1)      
        aminoAA2, aminoAA_count2, aminoAA_count_rel2, AA2, RAPP2 = prot_analysis_new(sequence2)
        df_prot1=pd.read_csv(StringIO(new_prot_data1), sep='\t')   
        df_prot2=pd.read_csv(StringIO(new_prot_data2), sep='\t') 
        return render_template('bar2.html', df_prot1=df_prot1, df_prot2=df_prot2, df_queryAC1=df_queryAC1, df_queryAC2=df_queryAC2, title='Hystogram', RAPP=RAPP, RAPP2=RAPP2, max=17000, AC1=user_input1, AC2=user_input2, cnt1=AA1, cnt2=AA2, barlabels1=aminoAA1, barlabels2=aminoAA2, barvalues1=aminoAA_count1, barvalues2=aminoAA_count2)
    except (KeyError, pd_errors.EmptyDataError):
        error='invalid AC'
        return render_template('protein_home.html', error=error)  
         
@app.route('/bar3', methods=['post'])
def bar3():
    AA=[]
    aminoAA=[]
    aminoAA1=[]
    aminoAA2=[]
    aminoAA3=[]
    countAA=[]
    aminoAA_count1=[]
    aminoAA_count2=[]
    aminoAA_count3=[]
    user_input1=request.form['input1']
    user_input2=request.form['input2']
    user_input3=request.form['input3']
    query1 = 'accession:'+user_input1
    query2 = 'accession:'+user_input2
    query3 = 'accession:'+user_input3
    #columnlist = "id, entry name,length,mass, families, citationmapping, citation, comment(DISEASE), comment(TEMPERATURE DEPENDENCE), comment(PH DEPENDENCE), comment(FUNCTION), comment(TISSUE SPECIFICITY), genes(PREFERRED), comment(SUBUNIT), comment(SUBCELLULAR LOCATION)"
    columnlist = "accession, id, length, mass, protein_families, cc_disease, temp_dependence, ph_dependence, cc_function, cc_subunit, cc_subcellular_location, sequence"
    service = UniProt()
    try:

        new_prot_data1 = service.search(query1, frmt="tsv", columns=columnlist)
        new_prot_data2 = service.search(query2, frmt="tsv", columns=columnlist)
        new_prot_data3 = service.search(query3, frmt="tsv", columns=columnlist)
        df_queryAC1 = pd.read_table(io.StringIO(new_prot_data1))
        df_queryAC2 = pd.read_table(io.StringIO(new_prot_data2))
        df_queryAC3 = pd.read_table(io.StringIO(new_prot_data3))
        sequence1=df_queryAC1.iloc[0]['Sequence']
        sequence2=df_queryAC2.iloc[0]['Sequence']
        sequence3=df_queryAC3.iloc[0]['Sequence']
        aminoAA1, aminoAA_count1, aminoAA_count_rel, AA1, RAPP = prot_analysis_new(sequence1)
        aminoAA2, aminoAA_count2, aminoAA_count_rel2, AA2, RAPP2 = prot_analysis_new(sequence2)
        aminoAA3, aminoAA_count3, aminoAA_count_rel3, AA3, RAPP3 = prot_analysis_new(sequence3)
        df_prot1=pd.read_csv(StringIO(new_prot_data1), sep='\t')   
        df_prot2=pd.read_csv(StringIO(new_prot_data2), sep='\t') 
        df_prot3=pd.read_csv(StringIO(new_prot_data3), sep='\t') 
        return render_template('bar3.html', df_prot1=df_prot1, df_prot2=df_prot2, df_prot3=df_prot3, df_queryAC1=df_queryAC1, df_queryAC2=df_queryAC2, df_queryAC3=df_queryAC3, title='Hystogram', RAPP=RAPP, RAPP2=RAPP2, RAPP3=RAPP3, max=17000, AC1=user_input1, AC2=user_input2, AC3=user_input3, cnt1=AA1, cnt2=AA2, cnt3=AA3, barlabels1=aminoAA1, barlabels2=aminoAA2, barlabels3=aminoAA3, barvalues1=aminoAA_count1, barvalues2=aminoAA_count2, barvalues3=aminoAA_count3)
    except (KeyError, pd_errors.EmptyDataError):
        error='invalid AC'
        return render_template('protein_home.html', error=error)  

    

@app.route('/bar4', methods=['post'])
def bar4():
    AA=[]
    aminoAA=[]
    aminoAA1=[]
    aminoAA2=[]
    aminoAA3=[]
    aminoAA4=[]
    countAA=[]
    aminoAA_count1=[]
    aminoAA_count2=[]
    aminoAA_count3=[]
    aminoAA_count4=[]
    user_input1=request.form['input1']
    user_input2=request.form['input2']
    user_input3=request.form['input3']
    user_input4=request.form['input4']
    query1 = 'accession:'+user_input1
    query2 = 'accession:'+user_input2
    query3 = 'accession:'+user_input3
    query4 = 'accession:'+user_input4
    #columnlist = "id, entry name,length,mass, families, citationmapping, citation, comment(DISEASE), comment(TEMPERATURE DEPENDENCE), comment(PH DEPENDENCE), comment(FUNCTION), comment(TISSUE SPECIFICITY), genes(PREFERRED), comment(SUBUNIT), comment(SUBCELLULAR LOCATION)"
    columnlist = "accession, id, length, mass, protein_families, cc_disease, temp_dependence, ph_dependence, cc_function, cc_subunit, cc_subcellular_location, sequence"
    service = UniProt()
    
    try:

        new_prot_data1 = service.search(query1, frmt="tsv", columns=columnlist)
        new_prot_data2 = service.search(query2, frmt="tsv", columns=columnlist)
        new_prot_data3 = service.search(query3, frmt="tsv", columns=columnlist)
        new_prot_data4 = service.search(query4, frmt="tsv", columns=columnlist)
        df_queryAC1 = pd.read_table(io.StringIO(new_prot_data1))
        df_queryAC2 = pd.read_table(io.StringIO(new_prot_data2))
        df_queryAC3 = pd.read_table(io.StringIO(new_prot_data3))
        df_queryAC4 = pd.read_table(io.StringIO(new_prot_data4))
        sequence1=df_queryAC1.iloc[0]['Sequence']
        sequence2=df_queryAC2.iloc[0]['Sequence']
        sequence3=df_queryAC3.iloc[0]['Sequence']
        sequence4=df_queryAC4.iloc[0]['Sequence']
        aminoAA1, aminoAA_count1, aminoAA_count_rel, AA1, RAPP = prot_analysis_new(sequence1)
        aminoAA2, aminoAA_count2, aminoAA_count_rel2, AA2, RAPP2 = prot_analysis_new(sequence2)
        aminoAA3, aminoAA_count3, aminoAA_count_rel3, AA3, RAPP3 = prot_analysis_new(sequence3)
        aminoAA4, aminoAA_count4, aminoAA_count_rel4, AA4, RAPP4 = prot_analysis_new(sequence4)
        df_prot1=pd.read_csv(StringIO(new_prot_data1), sep='\t')
        df_prot2=pd.read_csv(StringIO(new_prot_data2), sep='\t')
        df_prot3=pd.read_csv(StringIO(new_prot_data3), sep='\t')
        df_prot4=pd.read_csv(StringIO(new_prot_data4), sep='\t')
        return render_template('bar4.html', df_prot1=df_prot1, df_prot2=df_prot2, df_prot3=df_prot3, df_prot4=df_prot4, df_queryAC1=df_queryAC1, df_queryAC2=df_queryAC2, df_queryAC3=df_queryAC3, df_queryAC4=df_queryAC4, title='Hystogram', max=17000, RAPP=RAPP, RAPP2=RAPP2, RAPP3=RAPP3, RAPP4=RAPP4, AC1=user_input1, AC2=user_input2, AC3=user_input3, AC4=user_input4, cnt1=AA1, cnt2=AA2, cnt3=AA3, cnt4=AA4, barlabels1=aminoAA1, barlabels2=aminoAA2, barlabels3=aminoAA3, barlabels4=aminoAA4, barvalues1=aminoAA_count1, barvalues2=aminoAA_count2, barvalues3=aminoAA_count3, barvalues4=aminoAA_count4)

    except (KeyError, pd_errors.EmptyDataError):
        error='invalid AC'
        return render_template('protein_home.html', error=error)  

@app.route('/protein_count', methods=['post'])
def protein_count():
    import numpy
    AA=[]
    aminoAA=[]
    countAA=[]
    aminoAA_count=[]
    aminoAA_count_rel=[]
    user_input=request.form['input']
    query='accession:'+user_input
    columnlist = "id, entry name,length,mass, families, citationmapping, citation, comment(DISEASE), comment(TEMPERATURE DEPENDENCE), comment(PH DEPENDENCE), comment(FUNCTION), comment(TISSUE SPECIFICITY), genes(PREFERRED), comment(SUBUNIT), comment(SUBCELLULAR LOCATION)"

    EQ=[]
    EP=[]
    YF=[]
    DN=[]
    RAPP=[]

    service = UniProt()

    try:
        new_prot_data = service.search(query, frmt="tab", columns=columnlist)
        df_queryAC = pd.read_table(io.StringIO(new_prot_data))
        Entry_name=df_queryAC.iloc[0]["Entry name"] 

        aminoAA, aminoAA_count, aminoAA_count_rel, AA, RAPP = prot_analysis_new(user_input) 
        df_prot=pd.read_csv(StringIO(new_prot_data), sep='\t')    
        df = pd.DataFrame(list(zip(aminoAA, aminoAA_count, aminoAA_count_rel)), columns =['AA', 'Abs freq '+Entry_name, 'Rel freq '+Entry_name])

        protAAcount=df.to_csv()

        #protAAcountxls=write_excel(df)
        
        return render_template("protein_count.html", df_prot=df_prot, df_queryAC=df_queryAC, RAPP=RAPP, tables=[df.to_html(classes='data')], titles=df.columns.values, protAACnt=protAAcount, ACprot=user_input, AAcount=AA, aminoAA_res=aminoAA, aminoAA_count_res=aminoAA_count, aminoAA_count_rel_res=aminoAA_count_rel)
         
    except (KeyError, pd_errors.EmptyDataError):
        error='invalid AC'
        return render_template('protein_home.html', error=error)          
    
    
@app.route('/protein_count2', methods=['post'])
def protein_count2():

    AA=[]
    aminoAA=[]
    countAA=[]
    aminoAA_count=[]
    aminoAA_count_rel=[]
    aminoAA_count_rel2=[]

    user_input1=request.form['input']
    user_input2=request.form['input2']
    
    query1='accession:'+user_input1
    query2='accession:'+user_input2

    columnlist = "id, entry name,length,mass, families, citationmapping, citation, comment(DISEASE), comment(TEMPERATURE DEPENDENCE), comment(PH DEPENDENCE), comment(FUNCTION), comment(TISSUE SPECIFICITY), genes(PREFERRED), comment(SUBUNIT), comment(SUBCELLULAR LOCATION)"
    
    service = UniProt()
    
    EQ=[]
    EP=[]
    YF=[]
    DN=[]
    RAPP=[]

    try:
        new_prot_data1 = service.search(query1, frmt="tab", columns=columnlist)
        new_prot_data2 = service.search(query2, frmt="tab", columns=columnlist)

        df_queryAC1 = pd.read_table(io.StringIO(new_prot_data1))
        df_queryAC2 = pd.read_table(io.StringIO(new_prot_data2))

        Entry_name1=df_queryAC1.iloc[0]["Entry name"] 
        Entry_name2=df_queryAC2.iloc[0]["Entry name"] 
       
        aminoAA, aminoAA_count, aminoAA_count_rel, AA, RAPP = prot_analysis_new(user_input1)     
        aminoAA2, aminoAA_count2, aminoAA_count_rel2, AA2, RAPP2 = prot_analysis_new(user_input2) 
        df_prot1=pd.read_csv(StringIO(new_prot_data1), sep='\t')    
        df_prot2=pd.read_csv(StringIO(new_prot_data2), sep='\t')  
        
        df = pd.DataFrame(list(zip(aminoAA, aminoAA_count, aminoAA_count2, aminoAA_count_rel, aminoAA_count_rel2)), columns =['AA', 'Abs freq '+Entry_name1, 'Abs freq '+Entry_name2, 'Rel freq '+Entry_name1, 'Rel freq '+Entry_name2])
        protAAcount=df.to_csv()
        #protAAcountxls=write_excel(df)
        #return render_template("protein_count.html", df_prot1=df_prot1,  df_prot2=df_prot2, df_queryAC=df_queryAC, RAPP=RAPP, tables=[df.to_html(classes='data')], titles=df.columns.values, protAACnt=protAAcount, ACprot=user_input, AAcount=AA, aminoAA_res=aminoAA, aminoAA_count_res=aminoAA_count, aminoAA_count_rel_res=aminoAA_count_rel)
        return render_template("protein_count2.html", df_prot1=df_prot1,  df_prot2=df_prot2, df_queryAC1=df_queryAC1, df_queryAC2=df_queryAC2, RAPP=RAPP, RAPP2=RAPP2, tables=[df.to_html(classes='data')], titles=df.columns.values, protAACnt=protAAcount, ACprot=user_input1, ACprot2=user_input2, AAcount=AA, aminoAA_res=aminoAA, aminoAA_res2=aminoAA2, aminoAA_count_res=aminoAA_count, aminoAA_count_res2=aminoAA_count2)
 
        #return render_template("protein_count.html", df_prot1=df_prot1,  df_prot2=df_prot2, df_queryAC=df_queryAC, RAPP=RAPP, tables=[df.to_html(classes='data')], titles=df.columns.values, protAACnt=protAAcount, ACprot=user_input, AAcount=AA, aminoAA_res=aminoAA, aminoAA_count_res=aminoAA_count, aminoAA_count_rel_res=aminoAA_count_rel)
        #return render_template("protein_count2.html", df_prot1=df_prot1,  df_prot2=df_prot2, df_queryAC1=df_queryAC1, df_queryAC2=df_queryAC2, RAPP=RAPP, RAPP2=RAPP2, tables=[df.to_html(classes='data')], titles=df.columns.values, protAACnt=protAAcount, ACprot=user_input1, ACprot2=user_input2, AAcount=AA, aminoAA_res=aminoAA, aminoAA_res2=aminoAA2, aminoAA_count_res=aminoAA_count, aminoAA_count_res2=aminoAA_count2)
    except (KeyError, pd_errors.EmptyDataError):
        error='invalid AC'
        return render_template('protein_home.html', error=error)          



@app.route('/protein_count3', methods=['post'])
def protein_count3():
    
    AA=[]
    aminoAA=[]
    countAA=[]
    aminoAA_count=[]
    aminoAA_count_rel=[]
    aminoAA_count_rel2=[]
    aminoAA_count_rel3=[]

    user_input1=request.form['input']
    user_input2=request.form['input2']
    user_input3=request.form['input3']
    
    query1='accession:'+user_input1
    query2='accession:'+user_input2
    query3='accession:'+user_input3

    columnlist = "id, entry name,length,mass, families, citationmapping, citation, comment(DISEASE), comment(TEMPERATURE DEPENDENCE), comment(PH DEPENDENCE), comment(FUNCTION), comment(TISSUE SPECIFICITY), genes(PREFERRED), comment(SUBUNIT), comment(SUBCELLULAR LOCATION)"
    
    service = UniProt()
   
   
    EQ=[]
    EP=[]
    YF=[]
    DN=[]
    RAPP=[]
    

    try:
        new_prot_data1 = service.search(query1, frmt="tab", columns=columnlist)
        new_prot_data2 = service.search(query2, frmt="tab", columns=columnlist)
        new_prot_data3 = service.search(query3, frmt="tab", columns=columnlist)

        df_queryAC1 = pd.read_table(io.StringIO(new_prot_data1))
        df_queryAC2 = pd.read_table(io.StringIO(new_prot_data2))
        df_queryAC3 = pd.read_table(io.StringIO(new_prot_data3))

        Entry_name1=df_queryAC1.iloc[0]["Entry name"] 
        Entry_name2=df_queryAC2.iloc[0]["Entry name"] 
        Entry_name3=df_queryAC3.iloc[0]["Entry name"] 
        
        aminoAA, aminoAA_count, aminoAA_count_rel, AA, RAPP = prot_analysis_new(user_input1)     
        aminoAA2, aminoAA_count2, aminoAA_count_rel2, AA2, RAPP2 = prot_analysis_new(user_input2) 
        aminoAA3, aminoAA_count3, aminoAA_count_rel3, AA3, RAPP3 = prot_analysis_new(user_input3)  

        df_prot1=pd.read_csv(StringIO(new_prot_data1), sep='\t')    
        df_prot2=pd.read_csv(StringIO(new_prot_data2), sep='\t')  
        df_prot3=pd.read_csv(StringIO(new_prot_data3), sep='\t')

   
        df = pd.DataFrame(list(zip(aminoAA, aminoAA_count, aminoAA_count2, aminoAA_count3, aminoAA_count_rel, aminoAA_count_rel2, aminoAA_count_rel3)), columns =['AA', 'Abs freq '+Entry_name1, 'Abs freq '+Entry_name2, 'Abs freq '+Entry_name3, 'Rel freq '+Entry_name1, 'Rel freq '+Entry_name2, 'Rel freq '+Entry_name3])
        protAAcount=df.to_csv()
        #protAAcountxls=write_excel(df)
        #return render_template("protein_count2.html", df_prot1=df_prot1,  df_prot2=df_prot2, df_queryAC1=df_queryAC1, df_queryAC2=df_queryAC2, RAPP=RAPP, RAPP2=RAPP2, tables=[df.to_html(classes='data')], titles=df.columns.values, protAACnt=protAAcount, ACprot=user_input1, ACprot2=user_input2, AAcount=AA, aminoAA_res=aminoAA, aminoAA_res2=aminoAA2, aminoAA_count_res=aminoAA_count, aminoAA_count_res2=aminoAA_count2)
        return render_template("protein_count3.html", df_prot1=df_prot1,  df_prot2=df_prot2, df_prot3=df_prot3, df_queryAC1=df_queryAC1, df_queryAC2=df_queryAC2, df_queryAC3=df_queryAC3, RAPP=RAPP, RAPP2=RAPP2, RAPP3=RAPP3, tables=[df.to_html(classes='data')], titles=df.columns.values, protAACnt=protAAcount, ACprot=user_input1, ACprot2=user_input2, ACprot3=user_input3, AAcount=AA, aminoAA_res=aminoAA, aminoAA_res2=aminoAA2, aminoAA_count_res=aminoAA_count, aminoAA_count_res2=aminoAA_count2)

    except (KeyError, pd_errors.EmptyDataError):
        error='invalid AC'
        return render_template('protein_home.html', error=error)          


@app.route('/protein_count4', methods=['post'])
def protein_count4():

    AA=[]
    aminoAA=[]
    countAA=[]
    aminoAA_count=[]
    aminoAA_count_rel=[]
    aminoAA_count_rel2=[]
    aminoAA_count_rel3=[]
    aminoAA_count_rel4=[]

    user_input1=request.form['input']
    user_input2=request.form['input2']
    user_input3=request.form['input3']
    user_input4=request.form['input4']

    query1='accession:'+user_input1
    query2='accession:'+user_input2
    query3='accession:'+user_input3
    query4='accession:'+user_input4

    EQ=[]
    EP=[]
    YF=[]
    DN=[]
    RAPP=[]    

    columnlist = "id, entry name,length,mass, families, citationmapping, citation, comment(DISEASE), comment(TEMPERATURE DEPENDENCE), comment(PH DEPENDENCE), comment(FUNCTION), comment(TISSUE SPECIFICITY), genes(PREFERRED), comment(SUBUNIT), comment(SUBCELLULAR LOCATION)"
    
    service = UniProt()
   
    try:
        new_prot_data1 = service.search(query1, frmt="tab", columns=columnlist)
        new_prot_data2 = service.search(query2, frmt="tab", columns=columnlist)
        new_prot_data3 = service.search(query3, frmt="tab", columns=columnlist)
        new_prot_data4 = service.search(query4, frmt="tab", columns=columnlist)

        df_queryAC1 = pd.read_table(io.StringIO(new_prot_data1))
        df_queryAC2 = pd.read_table(io.StringIO(new_prot_data2))
        df_queryAC3 = pd.read_table(io.StringIO(new_prot_data3))
        df_queryAC4 = pd.read_table(io.StringIO(new_prot_data4))

        Entry_name1=df_queryAC1.iloc[0]["Entry name"] 
        Entry_name2=df_queryAC2.iloc[0]["Entry name"] 
        Entry_name3=df_queryAC3.iloc[0]["Entry name"] 
        Entry_name4=df_queryAC4.iloc[0]["Entry name"] 

        aminoAA, aminoAA_count, aminoAA_count_rel, AA, RAPP = prot_analysis_new(user_input1)     
        aminoAA2, aminoAA_count2, aminoAA_count_rel2, AA2, RAPP2 = prot_analysis_new(user_input2) 
        aminoAA3, aminoAA_count3, aminoAA_count_rel3, AA3, RAPP3 = prot_analysis_new(user_input3)  
        aminoAA4, aminoAA_count4, aminoAA_count_rel4, AA4, RAPP4 = prot_analysis_new(user_input4)

        df_prot1=pd.read_csv(StringIO(new_prot_data1), sep='\t')    
        df_prot2=pd.read_csv(StringIO(new_prot_data2), sep='\t')  
        df_prot3=pd.read_csv(StringIO(new_prot_data3), sep='\t')
        df_prot4=pd.read_csv(StringIO(new_prot_data4), sep='\t')
        df = pd.DataFrame(list(zip(aminoAA, aminoAA_count, aminoAA_count2, aminoAA_count3, aminoAA_count4, aminoAA_count_rel, aminoAA_count_rel2, aminoAA_count_rel3, aminoAA_count_rel4)), columns =['AA', 'Abs freq '+Entry_name1, 'Abs freq '+Entry_name2, 'Abs freq '+Entry_name3, 'Abs freq '+Entry_name4, 'Rel freq '+Entry_name1, 'Rel freq '+Entry_name2, 'Rel freq '+Entry_name3, 'Rel freq '+Entry_name4])

        protAAcount=df.to_csv()
        #return render_template("protein_count2.html", df_prot1=df_prot1,  df_prot2=df_prot2, df_queryAC1=df_queryAC1, df_queryAC2=df_queryAC2, RAPP=RAPP, RAPP2=RAPP2, tables=[df.to_html(classes='data')], titles=df.columns.values, protAACnt=protAAcount, ACprot=user_input1, ACprot2=user_input2, AAcount=AA, aminoAA_res=aminoAA, aminoAA_res2=aminoAA2, aminoAA_count_res=aminoAA_count, aminoAA_count_res2=aminoAA_count2)
        return render_template("protein_count4.html", df_prot1=df_prot1,  df_prot2=df_prot2, df_prot3=df_prot3, df_prot4=df_prot4, df_queryAC1=df_queryAC1, df_queryAC2=df_queryAC2, df_queryAC3=df_queryAC3, df_queryAC4=df_queryAC4, RAPP=RAPP, RAPP2=RAPP2, RAPP3=RAPP3, RAPP4=RAPP4, tables=[df.to_html(classes='data')], titles=df.columns.values, protAACnt=protAAcount, ACprot=user_input1, ACprot2=user_input2, ACprot3=user_input3, ACprot4=user_input4, AAcount=AA, aminoAA_res=aminoAA, aminoAA_res2=aminoAA2, aminoAA_count_res=aminoAA_count, aminoAA_count_res2=aminoAA_count2)

    except (KeyError, pd_errors.EmptyDataError):
        error='invalid AC'
        return render_template('protein_home.html', error=error)          

@app.route('/chromosome_count', methods = ['POST', 'GET'])
def chromosome_count():

	aminoAA.clear()
	aminoAA_count.clear()
	chromosome_number = request.args.get('chromosome')
	server = Server(host='http://www.ensembl.org')
	dataset = Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')
	
	Ensembl_data = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name', 'start_position', 'end_position', 'chromosome_name', 'uniprotswissprot', 'peptide'], filters={'chromosome_name': [chromosome_number]})
	Ensembl_data['UniProtKB/Swiss-Prot ID'].replace('', np.nan, inplace=True)
	Ensembl_data = Ensembl_data.dropna()
	Ensembl_data['Peptide'] = Ensembl_data['Peptide'].str.rstrip('*')
	Ensembl_data['length'] = Ensembl_data['Peptide'].apply(len)
	Ensembl_data.rename(columns={'UniProtKB/Swiss-Prot ID':'id', 'Peptide':'sequence'}, inplace=True)
	service = UniProt()
	#sostituito ......="Homo sapiens (Human) [9606]" con organism_id="9606"
	query1 = 'proteomecomponent:"chromosome '+chromosome_number+'" AND organism_id:"9606" AND proteome:up000005640 AND reviewed:true'
	columnlist_chromosome = "accession,id,length,sequence"
	#il formato tabulare non è più indicato txt in Bioservices ma con tsv
	result_uniprot = service.search(query1, frmt="tsv", columns=columnlist_chromosome)
	df_uniprot = pd.read_table(io.StringIO(result_uniprot))
	df_uniprot.rename(columns={'Entry':'id', 'Length':'length'}, inplace=True)
	#Attenzione, sono cambiati i nomi delle colonne in Bioservices
	#"entry name" adesso è indicato con "id", "id" con "accession", mentre lenght e sequence  sono rimasti invariati
	df_uniprot.columns=['id','entry name','length','sequence']
	merge_col = pd.merge(Ensembl_data, df_uniprot, on=['id', 'length', 'sequence'], how='inner')
	AA_count_list.clear()
	for i, row in enumerate(merge_col.values):
		countAA = ProteinAnalysis(merge_col.loc[i]['sequence'])
		AA=countAA.count_amino_acids()
		AA_count_list.append(AA)
	df_count = pd.DataFrame(AA_count_list)
	Chromosome_count_AA = pd.concat([merge_col, df_count], axis=1)
	columns=list(Chromosome_count_AA.columns.values)
	Chromosome_count_AA=Chromosome_count_AA[['Gene stable ID', 'Gene name', 'Gene start (bp)', 'Gene end (bp)', 'id', 'length', 'entry name', 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'Chromosome/scaffold name', 'sequence']]
	Chromosome_count_AA=Chromosome_count_AA.sort_values(by=['Gene start (bp)'])
	df_html_chr=Chromosome_count_AA.to_html()
	chr_csv=Chromosome_count_AA.to_csv()
	chr_list=Chromosome_count_AA.values.tolist()

	def walk(Amino, x, T):
		count = []
		i=1
		j=x
		while j<=Amino.size-x:
			subset=Amino[i:j+1]
			countlessT=subset[subset < T].count()
			count.append(countlessT)
			i=i+1
			j=j+1
		return count
	e = walk(df_count['E'], 25, 49)
	e_column=e
	q = walk(df_count['Q'], 25, 49)
	q_column=q

	walk = {'E WALK': e_column,
		'Q WALK': q_column
		 }

	df2 = pd.DataFrame.from_dict(walk)
	Chromosome_count_AA_WALK = pd.concat([merge_col, df_count, df2], axis=1)
	Chromosome_count_AA_WALK=Chromosome_count_AA_WALK.sort_values(by=['Gene start (bp)'])
	

	Chromosome_count_AA_top = Chromosome_count_AA.head()
	return render_template("chromosome_count.html", chr_list_web=chr_list, df_html_chr=df_html_chr, df_chrom=Chromosome_count_AA, Chrom_top=Chromosome_count_AA_top, res_Uniprot=result_uniprot, chr_csv_get=chr_csv, res_Ensambl=Ensembl_data, chromosome_num=chromosome_number)


