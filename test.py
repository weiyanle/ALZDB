10
#coding=utf-8
import webbrowser
import os
from flask import Flask, render_template ,request, redirect, url_for, flash, json
from flask_wtf import Form
from wtforms import StringField, SubmitField, SelectField
from wtforms.validators import DataRequired
import pymysql
from decimal import *
import re
import csv


app = Flask(__name__)

app.secret_key = 'test'

class Form1(Form):
	rs_id=StringField(render_kw={'placeholder':'rsID'})
	submit1=SubmitField('Search')

class Form2(Form):
	gene_id=StringField(render_kw={'placeholder':'SYMBOL/ENTREZ'})
	stp=SelectField(choices=[('all','All'),('GWAS','Genome'),('trans','Transcriptome'),('pro','Proteome'),('methy','Epigenome')],coerce = str)
	submit2=SubmitField('Search')

class Form3(Form):
	position=StringField(render_kw={'placeholder':'position'})
	submit3=SubmitField('Search')

class Form4(Form):
	input_gse=StringField(render_kw={'placeholder':'GSE'})
	submit4=SubmitField('Search')

@app.route('/',methods=['GET','POST'])
def index():
	form1=Form1()
	form2=Form2()
	form3=Form3()
	form4=Form4()
	if form1.submit1.data and form1.validate():
		tp='rsID'
		rsid=form1.rs_id.data
		return redirect('/search/type='+tp+'/input='+rsid)
	if form2.submit2.data and form2.validate():
		geneid=form2.gene_id.data
		searchtp=form2.stp.data
		tp=searchtp+'_geneID'
		return redirect('/search/type='+tp+'/input='+geneid)
	if form3.submit3.data and form3.validate():
		pos=form3.position.data
		searchtp=form2.stp.data
		tp=searchtp+'_position'
		return redirect('/search/type='+tp+'/input='+pos)
	if form4.submit4.data and form4.validate():
		inputgse='GSE'+form4.input_gse.data
		tp='gse'
		return redirect('/search/type='+tp+'/input='+inputgse)
	return render_template('index.html',form1=form1,form2=form2,form3=form3,form4=form4)

@app.route('/GWAS/chr<chro>:<pos>',methods=['GET','POST'])
def choosegwas(chro,pos):
	return render_template('gwas.html',chro=chro,pos=pos)

@app.route('/GWAS/chr<chro>',methods=['GET','POST'])
def gwaschr(chro):
	fp={'1':'751343','2':'10587','3':'60363','4':'68786','5':'12114','6':'192106','7':'35131','8':'91936','9':'46587','10':'66397','11':'131128','12':'181266','13':'19198288','14':'20385550','15':'20150882','16':'79609','17':'12344','18':'11671','19':'252541','20':'	61270','21':'10971951','22':'17051775'}
	return redirect('/GWAS/chr'+chro+':'+fp[chro])


@app.route('/miRNA',methods=['GET','POST'])
def choosemirna():
	data = getsummary('mirna')
	return render_template('choosemirna.html',data=data)


@app.route('/mRNA',methods=['GET','POST'])
def choosemrna():
	data = getsummary('mrna')
	return render_template('choosemrna.html',data=data)

@app.route('/protein',methods=['GET','POST'])
def chooseprotein():
	data = getsummary('protein')
	return render_template('chooseprotein.html',data=data)

@app.route('/epigenome',methods=['GET','POST'])
def epigenome():
	return render_template('epigenome.html')

@app.route('/meta',methods=['GET','POST'])
def choosemeta():
	return render_template('choosemeta.html')

@app.route('/meta/<tissue>/index=<dir1>/<dir2>',methods=['GET','POST'])
def meta(tissue,dir1,dir2):
	return render_template('meta.html',tissue=tissue,dir1=dir1,dir2=dir2)

@app.route('/search/type=<tp>/input=<search_input>',methods=['GET','POST'])
def search(tp,search_input):
	poslist=[]
	if 'geneID' in tp:
		poslist=findpos(search_input)
	if '_' in tp:
		tp0=tp.split("_")[0]
		tp1=tp.split("_")[1]
	else:
		tp0=tp
		tp1=tp
	return render_template('search.html',type=tp,tp0=tp0,tp1=tp1,input=search_input,poslist=poslist)

@app.route('/json/search/type=<tp>/input=<search_input>',methods=['GET','POST'])
def jsonsearch(tp,search_input):
	db = pymysql.connect(
		host='101.132.188.227',
		port=3306,
		user='test',
		passwd='123456',
		db='ADDB')
	cur = db.cursor()
	if tp == 'rsID':
		search_rs_id = "select * from gwas where rsID='%s';" % (search_input)
		cur.execute(search_rs_id)
		table = cur.fetchall()

	if tp == 'chrID':
		search_chr = "select * from gwas where chrID='%s' LIMIT 1000;" % (search_input)
		cur.execute(search_chr)
		table = cur.fetchall()
		
	if tp == 'gse' or tp == 'gse_all':
		if search_input not in ['GSE16759_GPL8757','GSE46131','GSE46579','GSE48552']:
			if search_input=='GSE84422':
				if tp == 'gse':
					search_gse = "select * from mRNA where GSE='GSE84422_GPL96' or GSE='GSE84422_GPL97' or GSE='GSE84422_GPL570';" 
				else:
					search_gse = "select * from mRNA_all where GSE='GSE84422_GPL96' or GSE='GSE84422_GPL97' or GSE='GSE84422_GPL570';"				
			else:	
				if tp == 'gse':
					search_gse = "select * from mRNA where GSE='%s' union select * from protein where GSE='%s';" % (search_input,search_input)
				else:
					search_gse = "select * from mRNA_all where GSE='%s' union select * from protein_all where GSE='%s';" % (search_input,search_input)
			cur.execute(search_gse)
			table = cur.fetchall()
		else:
			if tp == 'gse':
				search_gse = "select * from miRNA where GSE='%s';" % (search_input)
				cur.execute(search_gse)
				table = cur.fetchall()
			else:
				search_gse = "select * from miRNA_all where GSE='%s';" % (search_input)
				cur.execute(search_gse)
				table = cur.fetchall()

	tps=tp.split('_')

	if tps[0] == 'GWAS' and tps[1]== 'position':
		tmp_list =search_input.split(":")
		chro = int((tmp_list[0].split('chr'))[1])
		input_location = tmp_list[1]
		if '-' in input_location:
			input_location1=input_location.split('-')[0]
			input_location2=input_location.split('-')[1]
			search_position= "select * from gwas where chrID ='%d' and position>='%d' and position<='%d';" % (chro,int(input_location1),int(input_location2))
		else:
			search_position= "select * from gwas where chrID ='%d' and position='%d';" % (chro,int(input_location))
		cur.execute(search_position)
		table = cur.fetchall()

	if tps[0] == 'GWAS' and tps[1]=='geneID':
		if search_input.isdigit():
			search_input=int(search_input)
			search_geneID = "select chr,start,end from gene_location where geneID='%d';" % (search_input)
		else:
			search_geneID = "select chr,start,end from gene_location where SYMBOL='%s';" % (search_input)
		cur.execute(search_geneID)
		data1 = cur.fetchall()
		chro=data1[0][0][3:]
		search_geneID = "select * from gwas where chrID='%s' and position >= '%d' and position <= '%d';" % (chro,data1[0][1],data1[0][2])
		cur.execute(search_geneID)
		table = cur.fetchall()

	if tps[0] == 'trans' or tps[0] =='pro':
		diff_dic={'trans':'mRNA','pro':'protein'}
		if_digit=search_input.isdigit()
		id_dic={True:'geneID',False:'SYMBOL'}
		if tps[1] =='geneID' :
			if len(tps)==2:
				search_geneID = "select * from "+diff_dic[tps[0]]+" where "+id_dic[if_digit]+"='%s';" % (search_input)
			else:
				search_geneID = "select * from "+diff_dic[tps[0]]+"_all where "+id_dic[if_digit]+"='%s';" % (search_input)
			cur.execute(search_geneID)
			table = cur.fetchall()
		elif tps[1] =='position':	
			data1=containgene(search_input)
			table=[]
			if len(tps)==2:
				for dt in data1:
					search_position = "select * from "+diff_dic[tps[0]]+" where geneID='%s' or SYMBOL='%s';" % (dt[0],dt[1])
					cur.execute(search_position)
					data2 = cur.fetchall()
					table=table+list(data2)
			else:
				for dt in data1:
					search_position = "select * from "+diff_dic[tps[0]]+"_all where geneID='%s' or SYMBOL='%s';" % (dt[0],dt[1])
					cur.execute(search_position)
					data2 = cur.fetchall()
					table=table+list(data2)

	if tps[0] == 'methy' and tps[1]=='geneID':
		if search_input.isdigit():
			search_input=int(search_input)
			search_geneID = "select chr,start,end from gene_location where geneID='%d';" % (search_input)
		else:
			search_geneID = "select chr,start,end from gene_location where SYMBOL='%s';" % (search_input)
		cur.execute(search_geneID)
		data1 = cur.fetchall()
		search_geneID = "select * from DMP_GSE66351 where chr='%s' and pos >= '%d' and pos <= '%d' union select * from DMP_GSE67419_new where chr='%s' and pos >= '%d' and pos <= '%d';" % (data1[0][0],data1[0][1],data1[0][2],data1[0][0],data1[0][1],data1[0][2])
		cur.execute(search_geneID)
		table = cur.fetchall()

	if tps[0] == 'methy' and tps[1]=='position':
		tmp_list =search_input.split(":")
		chro = tmp_list[0]
		input_location = tmp_list[1]
		if '-' in input_location:
			input_location1=input_location.split('-')[0]
			input_location2=input_location.split('-')[1]
			search_position = "select * from DMP_GSE66351 where chr='%s' and pos >= '%d' and pos <= '%d' union select * from DMP_GSE67419_new where chr='%s' and pos >= '%d' and pos <= '%d';" % (chro,int(input_location1),int(input_location2),chro,int(input_location1),int(input_location2))
		else:
			search_position = "select * from DMP_GSE66351 where chr='%s' and pos = '%d' union select * from DMP_GSE67419_new where chr='%s' and pos = '%d';" % (chro,int(input_location),chro,int(input_location))
		cur.execute(search_position)
		table = cur.fetchall()
	db.close()

	jtable=[]
	if tp=='rsID' or tp=='chrID' or tps[0]=='GWAS':
		for row in table:
			jtable.append({'chrID':row[1],'position':row[2],'rsID':row[3],'before':row[4],'after':row[5],'beta':row[6],'unknown':row[7],'P_Value':row[8],'link':'view in jBrowse'})
	elif (tps[0]=='gse' or tps[0]=='trans' or tps[0]=='pro') and table:
		if len(table[0])==13:
			for row in table:
				jtable.append({'ENTREZ':row[1],'SYMBOL':row[2],'logFC':rounding(row[3]),'AveExpr':rounding(row[4]),'t':rounding(row[5]),'P_Value':rounding(row[6]),'adj_P_Value':rounding(row[7]),'B':rounding(row[8]),'tissue':tissue(row[9]),'GSE':row[10],'diff':row[11],'type':row[12],'link':'view in jBrowse'})
		else:
			for row in table:
				jtable.append({'miRNA_ID':row[1],'logFC':rounding(row[2]),'AveExpr':rounding(row[3]),'t':rounding(row[4]),'P_Value':rounding(row[5]),'adj_P_Value':rounding(row[6]),'B':rounding(row[7]),'tissue':tissue(row[8]),'GSE':row[9],'diff':row[10],'type':row[11]})
	elif tps[0]=='methy':
		for row in table:
			jtable.append({'chr':row[1],'pos':row[2],'strand':row[3],'Name':row[4],'Islands_Name':row[11],'Relation_to_Island':row[12],'UCSC_RefGene_Name':row[13],'UCSC_RefGene_Accession':row[14],'UCSC_RefGene_Group':row[15],'Phantom':row[16],'DMR':row[17],'Enhancer':row[18],'HMM_Island':row[19],'Regulatory_Feature_Name':row[20],'Regulatory_Feature_Group':row[21],'DHS':row[22],'logFC':row[23],'AveExpr':row[24],'t':row[25],'P_Value':rounding(row[26]),'adj_P_Val':rounding(row[27]),'B':row[28],'GSE':row[29]})
	result = json.dumps(jtable)	

	if tps[0]=='diff':
		result=json.dumps(json.loads(jsonsearch('trans'+tp[4:],search_input))+json.loads(jsonsearch('pro'+tp[4:],search_input)))
	return result

@app.route('/json/search/methy_high/type=<tp>/input=<search_input>',methods=['GET','POST'])
def jsonmythyhighsearch(tp,search_input):
	tps=tp.split('_')
	if tps[1]=='geneID':
		db = pymysql.connect(
				host='101.132.188.227',
				port=3306,
				user='test',
				passwd='123456',
				db='ADDB')
		cur = db.cursor()
		if search_input.isdigit():
			search_geneID = "select chr,start,end from gene_location where geneID='%s';" % (search_input)
		else:
			search_geneID = "select chr,start,end from gene_location where SYMBOL='%s';" % (search_input)
		cur.execute(search_geneID)
		data1 = cur.fetchall()
		search_geneID = "select * from methyDiff_notable_hyper where chr='%s' and start >= '%d' and start <= '%d' union select * from methyDiff_notable_hypo where chr='%s' and start >= '%d' and start <= '%d';" % (data1[0][0],data1[0][1],data1[0][2],data1[0][0],data1[0][1],data1[0][2])
		cur.execute(search_geneID)
		table = cur.fetchall()
		db.close()

	if tps[1]=='position':
		tmp_list =search_input.split(":")
		chro = int((tmp_list[0].split('chr'))[1])
		input_location = tmp_list[1]
		db = pymysql.connect(
		host='101.132.188.227',
		port=3306,
		user='test',
		passwd='123456',
		db='ADDB')
		cur = db.cursor()
		if '-' in input_location:
			input_location1=input_location.split('-')[0]
			input_location2=input_location.split('-')[1]
			search_position = "select * from methyDiff_notable_hyper where chr='%s' and start >= '%s' and start <= '%s' union select * from methyDiff_notable_hypo where chr='%s' and start >= '%s' and start <= '%s';" % (tmp_list[0],input_location1,input_location2,tmp_list[0],input_location1,input_location2)
		else:
			search_position = "select * from methyDiff_notable_hyper where chr='%s' and start = '%d' union select * from methyDiff_notable_hypo where chr='%s' and start = '%d';" % (tmp_list[0],int(input_location),tmp_list[0],int(input_location))
		cur.execute(search_position)
		table = cur.fetchall()
		db.close()

	jtable=[]
	for row in table:
		jtable.append({'chr':row[1],'start':row[2],'end':row[3],'strand':row[4],'pvalue':rounding(row[5]),'qvalue':rounding(row[6]),'meth_diff':row[7],'GSE':'GSE46644','link':'view in jBrowse'})		
	result = json.dumps(jtable)	
	return result

@app.route('/json/search/rate/input=<search_input>',methods=['GET','POST'])
def jsonrate(search_input):
	db = pymysql.connect(
		host='101.132.188.227',
		port=3306,
		user='test',
		passwd='123456',
		db='ADDB')
	cur = db.cursor()
	if search_input.isdigit():
		search_geneID1 = "select * from mRNA where geneID='%s' union select * from protein where geneID='%s';" % (search_input,search_input)
		search_geneID2 = "select * from mRNA_all where geneID='%s' union select * from protein_all where geneID='%s';" % (search_input,search_input)
	else:
		search_geneID1 = "select * from mRNA where SYMBOL='%s' union select * from protein where SYMBOL='%s';" % (search_input,search_input)
		search_geneID2 = "select * from mRNA_all where SYMBOL='%s' union select * from protein_all where SYMBOL='%s';" % (search_input,search_input)
	cur.execute(search_geneID1)
	table1 = cur.fetchall()
	cur.execute(search_geneID2)
	table2 = cur.fetchall()
	db.close()	
	result = json.dumps([{'diff':len(table1),'all':len(table2)}])	
	return result

@app.route('/json/epi/GSE<gse>/<tp>',methods=['GET','POST'])
def jsonepi(gse,tp):
	if tp=='dmps':
		f=open(os.path.abspath('.')+'/static/series/methylation/GSE'+gse+'/json_DMPs_part.txt')
	elif tp=='dmrs':
		f=open(os.path.abspath('.')+'/static/series/methylation/GSE'+gse+'/json_DMRs.txt')
	elif tp=='methydiff':
		f=open(os.path.abspath('.')+'/static/series/methylation/GSE'+gse+'/json_methyDiff_notable_all.txt')
	st=f.read()
	return st

@app.route('/epigenome/GSE<gse>/index=<num>/<diff_num>',methods=['GET','POST'])
def epiindex(gse,num,diff_num):
	return render_template('epigse.html',gse=gse,num=num,diff_num=diff_num)

@app.route('/GSE<gse>/index=<num>/<diff_num>',methods=['GET','POST'])
def gseindex(gse,num,diff_num):
	gselist=gsecontent(gse,num,diff_num)
	tp=gselist[0]
	file_list=gselist[1]
	diff_list=gselist[2]
	all_file_list=gselist[3]
	all_list=gselist[4]
	if_file=gselist[5]
	if_diff=gselist[6]
	if_all=gselist[7]
	return render_template('gse.html',tp=tp,file_list=file_list,diff_list=diff_list,all_file_list=all_file_list,all_list=all_list,gse=gse,num=num,if_file=if_file,if_diff=if_diff,if_all=if_all,diff_num=diff_num)

@app.route('/GSE<gse>/index=<num>/<diff_num>/<entrez>',methods=['GET','POST'])
def gseentrez(gse,num,diff_num,entrez):
	gselist=gsecontent(gse,num,diff_num)
	tp=gselist[0]
	file_list=gselist[1]
	diff_list=gselist[2]
	all_file_list=gselist[3]
	all_list=gselist[4]
	if_file=gselist[5]
	if_diff=gselist[6]
	if_all=gselist[7]
	sym=''
	page='1'
	if '_all' not in diff_num and len(diff_list[file_list.index(diff_num)])==1:
		pos='0'
		sym='Not Found'
	else:
		if entrez=='0':
			if '_all' in diff_num:
				for en in all_list[all_file_list.index(diff_num[0:-4])]:
					if en['ENTREZ'].isdigit():
						entrez=en['ENTREZ']
						break;
			else:
				for en in diff_list[file_list.index(diff_num)]:
					if en['ENTREZ'].isdigit():
						entrez=en['ENTREZ']
						break;
		if '_all' in diff_num:
			for gene in all_list[all_file_list.index(diff_num[0:-4])]:
				if gene['ENTREZ']==entrez:
					page=(all_list[all_file_list.index(diff_num[0:-4])].index(gene)-1)//10+1
					break;
		else:
			for gene in diff_list[file_list.index(diff_num)]:
				if gene['ENTREZ']==entrez:
					page=(diff_list[file_list.index(diff_num)].index(gene)-1)//10+1
					break;
		if entrez.isdigit():
			db = pymysql.connect(
					host='101.132.188.227',
					port=3306,
					user='test',
					passwd='123456',
					db='ADDB')
			cur = db.cursor()
			search_geneID = "select chr,start,end,SYMBOL from gene_location where geneID='%s';" % (entrez)
			cur.execute(search_geneID)
			table = cur.fetchall()
			db.close()
			jbdic={'AD_neurons--fetal_brain.GSE34879':'1','AD_neurons--normal_neurons.GSE34879':'2','AD--control.GSE45596':'3','AD--control.GSE110226':'4','AD--control.GSE32645':'5','AD--control.GSE6613':'6','AD--control.GSE85426':'7','AD--control.GSE95587':'8','AD--control.GSE95810':'9','AD--MCI.GSE18309':'10','AD--NC.GSE18309':'11','AD--Normal.GSE16759':'12','AD--ps_AD.GSE95810':'13','AD_CA1--AD_CA3.GSE29378':'14','AD_CA1--Control_CA1.GSE29378':'15','AD_CA3--Control_CA3.GSE29378':'16','AD_hIPSC--control_hIPSCs.GSE34879':'17','AD--normal.GSE33000':'18','Braak_6--Braak_3.GSE61196':'19','Brain_AD--control.GSE30945':'20','AD_mutation--control.GSE39420':'21','AD_early_stage--control.GSE39420':'22','AD_mutation--AD_early_stage.GSE39420':'23','AD1--control1.GSE104141':'24','Brain_AD--Breast_NC.GSE30945':'25','Brain_AD--Colon_NC.GSE30945':'26','AD--control.GSE12685':'27','AD--control.GSE15222':'28','AD--control.GSE26927':'29','AD--control.GSE43326':'30','AD--control.GSE63060':'31','AD--control.GSE63061':'32','AD--control.GSE97760':'33','AD--control.GSE53697':'34','MCI--control.GSE63060':'35','MCI--control.GSE63061':'36','Control_CA1--Control_CA3.GSE29378':'37','AD2--control2.GSE104141':'38','C_AD--C_Control.GSE6834':'39','EC_AD--control.GSE5281':'40','EC_AD--EC_control.GSE48350':'41','FC_definitive--FC_normal.GSE13214':'42','FC_definitive--FC_pathologic.GSE13214':'43','FC_pathologic--FC_normal.GSE13214':'44','Braak_3--control.GSE61196':'45','Braak_6--control.GSE61196':'46','HIP_AD--control.GSE5281':'47','HI_AD--HI_control.GSE48350':'48','HI_definitive--FC_normal.GSE13214':'49','HI_definitive--FC_pathologic.GSE13214':'50','HI_pathologic--FC_normal.GSE13214':'51','V_VI--III_IV.GSE29652':'52','Incipient--Control.GSE1297':'53','iPS26B--H1.GSE42492':'54','iPS26B--H9.GSE42492':'55','iPS5--H1.GSE42492':'56','iPS5--H9.GSE42492':'57','III_IV--I_II.GSE29652':'58','AD--MCI.GSE63060':'59','AD--MCI.GSE63061':'60','MCI--NC.GSE18309':'61','Brain_AD--Mixture_NC.GSE30945':'62','Moderate--Control.GSE1297':'63','Moderate--Incipient.GSE1297':'64','MTG_AD--control.GSE5281':'65','AD_Female--NEC_Female.GSE4226':'66','AD_Female--NEC_Female.GSE4229':'67','AD_Male--NEC_Male.GSE4226':'68','AD_Male--NEC_Male.GSE4229':'69','NFH2--NFH46.GSE42492':'70','NFH2--OiPS3.GSE42492':'71','NFH2--OiPS6.GSE42492':'72','NFH46--iPS26B.GSE42492':'73','NFH46--iPS5.GSE42492':'74','Definite--control.GSE84422_GPL570':'75','Definite--control.GSE84422_GPL96':'76','Definite--control.GSE84422_GPL97':'77','Possible--control.GSE84422_GPL570':'78','Possible--control.GSE84422_GPL96':'79','Possible--control.GSE84422_GPL97':'80','Probable--control.GSE84422_GPL570':'81','Probable--control.GSE84422_GPL96':'82','Probable--control.GSE84422_GPL97':'83','Tangle--control.GSE4757':'84','OiPS3--H1.GSE42492':'85','OiPS3--H9.GSE42492':'86','OiPS3--iPS26B.GSE42492':'87','OiPS3--iPS5.GSE42492':'88','OiPS6--H1.GSE42492':'89','OiPS6--H9.GSE42492':'90','OiPS6--iPS26B.GSE42492':'91','OiPS6--iPS5.GSE42492':'92','AD--old_control.GSE104704':'93','old--young.GSE104704':'94','PC_AD--control.GSE5281':'95','Definite--Possible.GSE84422_GPL570':'96','Definite--Possible.GSE84422_GPL96':'97','Definite--Possible.GSE84422_GPL97':'98','Possible--Probable.GSE84422_GPL570':'99','Possible--Probable.GSE84422_GPL96':'100','Possible--Probable.GSE84422_GPL97':'101','PG_AD--PG_control.GSE48350':'102','ps_AD--control.GSE95810':'103','Definite--Probable.GSE84422_GPL570':'104','Definite--Probable.GSE84422_GPL96':'105','Definite--Probable.GSE84422_GPL97':'106','Severe--Control.GSE1297':'107','Severe--Incipient.GSE1297':'108','Severe--Moderate.GSE1297':'109','SFG_AD--control.GSE5281':'110','SFG_AD--SFG_control.GSE48350':'111','TC_AD--TC_Control.GSE6834':'112','VCX_AD--control.GSE5281':'113','V_VI--I_II.GSE29652':'114','AD--young_control.GSE104704':'115'}
			jbdic2={'AD--control.GSE39087':'1','AD--control.GSE62283':'2','AD--control.GSE74763':'3','AD--Older_control.GSE29654':'4','AD--Older_control.GSE29676':'5','AD--Younger_control.GSE29654':'6','AD--Younger_control.GSE29676':'7'}
			if '_all' in diff_num:
				tpnum=diff_num[0:-4]
			else:
				tpnum=diff_num
			if tpnum+'.GSE'+gse in jbdic:
				if table:
					if table[0][1]<10000:
						pos='http://106.14.221.173/jbrowse/?loc='+table[0][0]+'%3A0..'+str(table[0][2]+table[0][1])+'&tracks=DNA%2Cmicroarray%20data%20'+jbdic[tpnum+'.GSE'+gse]+'%2CGRCh37&highlight='
					else:
						pos='http://106.14.221.173/jbrowse/?loc='+table[0][0]+'%3A'+str(table[0][1]-10000)+'..'+str(table[0][2]+10000)+'&tracks=DNA%2Cmicroarray%20data%20'+jbdic[tpnum+'.GSE'+gse]+'%2CGRCh37&highlight='
					sym=table[0][3]
				else:
					pos='http://106.14.221.173/jbrowse/?loc=chr1%3A10000..30000&tracks=DNA%2Cmicroarray%20data%20'+jbdic[tpnum+'.GSE'+gse]+'%2CGRCh37&highlight='
					sym='Not Found'
			elif tpnum+'.GSE'+gse in jbdic2:
				if table:
					if table[0][1]<10000:
						pos='http://106.14.221.173/jbrowse/?loc='+table[0][0]+'%3A0..'+str(table[0][2]+table[0][1])+'&tracks=DNA%2Cmicroarray%20data%20'+jbdic2[tpnum+'.GSE'+gse]+'%2CGRCh37&highlight='
					else:
						pos='http://106.14.221.173/jbrowse/?loc='+table[0][0]+'%3A'+str(table[0][1]-10000)+'..'+str(table[0][2]+10000)+'&tracks=DNA%2Cmicroarray%20data%20'+jbdic2[tpnum+'.GSE'+gse]+'%2CGRCh37&highlight='
					sym=table[0][3]
				else:
					pos='http://106.14.221.173/jbrowse/?loc=chr1%3A10000..30000&tracks=DNA%2Cmicroarray%20data%20'+jbdic2[tpnum+'.GSE'+gse]+'%2CGRCh37&highlight='
					sym='Not Found'
			else:
				if table:
					if table[0][1]<10000:
						pos='http://106.14.221.173/jbrowse/?loc='+table[0][0]+'%3A0..'+str(table[0][2]+table[0][1])+'&tracks=GRCh37&highlight='
					else:
						pos='http://106.14.221.173/jbrowse/?loc='+table[0][0]+'%3A'+str(table[0][1]-10000)+'..'+str(table[0][2]+10000)+'&tracks=GRCh37&highlight='
					sym=table[0][3]
				else:
					pos='http://106.14.221.173/jbrowse/?loc=chr1%3A10000..30000&tracks=GRCh37&highlight='
					sym='Not Found'				
		else:
			pos='0'
			sym='Not Found'
	return render_template('gse.html',tp=tp,file_list=file_list,diff_list=diff_list,all_file_list=all_file_list,all_list=all_list,gse=gse,num=num,if_file=if_file,if_diff=if_diff,if_all=if_all,diff_num=diff_num,ent=entrez,pos=pos,sym=sym,page=page)


@app.route('/json/GSE<gse>/<diff_num>',methods=['GET','POST'])
def jsongse(gse,diff_num):
	all_gse={'104141': 'mRNA', '104704': 'mRNA', '110226': 'mRNA', '12685': 'mRNA', '1297': 'mRNA', '13214': 'mRNA', '15222': 'mRNA', '16759': 'mRNA', '18309': 'mRNA', '26927': 'mRNA', '26972': 'mRNA', '28146': 'mRNA', '29378': 'mRNA', '29652': 'mRNA', '30945': 'mRNA', '32645': 'mRNA', '33000': 'mRNA', '34879': 'mRNA', '36980': 'mRNA', '37263': 'mRNA', '37264': 'mRNA', '39420': 'mRNA', '4226': 'mRNA', '4229': 'mRNA', '42492': 'mRNA', '43326': 'mRNA', '44768': 'mRNA', '44770': 'mRNA', '44771': 'mRNA', '45596': 'mRNA', '4757': 'mRNA', '48350': 'mRNA', '5281': 'mRNA', '53697': 'mRNA', '61196': 'mRNA', '63060': 'mRNA', '63061': 'mRNA', '6613': 'mRNA', '6834': 'mRNA', '84422_GPL570': 'mRNA', '84422_GPL96': 'mRNA', '84422_GPL97': 'mRNA', '85426': 'mRNA', '93885': 'mRNA', '95587': 'mRNA', '95810': 'mRNA', '97760': 'mRNA', '29654': 'protein', '29676': 'protein', '39087': 'protein', '62283': 'protein', '70424': 'protein', '74763': 'protein','16759_GPL8757':'miRNA','46131':'miRNA','46579':'miRNA','48552':'miRNA'}
	if gse in all_gse:
		tp=all_gse[gse]
	else:
		tp='mRNA'
	abspath_copy=os.path.abspath('.')
	filepath=os.path.abspath('.')+'/static/series/'+tp+'/GSE'+gse+'/diff.'+diff_num+'.GSE'+gse+'/diff.'+diff_num+'.GSE'+gse+'.txt'
	all_filepath=os.path.abspath('.')+'/static/series/'+tp+'/ALL/GSE'+gse+'/all.'+diff_num[0:-4]+'.GSE'+gse+'.txt'
	genes=[]
	if '_all' in diff_num:
		file=open(all_filepath)
	else:
		file=open(filepath)
	diff=file.read()
	diffs=re.split('\t\n|\t|\n',diff)
	if tp == 'miRNA':
		if len(diffs)>10:
			for i in range(1,int((len(diffs))/7)):
				genes.append({
					'ENTREZ':'-',
					'SYMBOL':diffs[7*i-1],
					'logFC':diffs[7*i],
					'AveExpr':diffs[7*i+1],
					't':diffs[7*i+2],
					'P_Value':diffs[7*i+3],
					'adj_P_Val':diffs[7*i+4],
					'B':diffs[7*i+5],
					'link':'view in jBrowse',
					})
	else:
		if len(diffs)>10:
			for i in range(1,int((len(diffs))/8)):
				genes.append({
					'ENTREZ':diffs[8*i],
					'SYMBOL':diffs[8*i+1],
					'logFC':diffs[8*i+2],
					'AveExpr':diffs[8*i+3],
					't':diffs[8*i+4],
					'P_Value':diffs[8*i+5],
					'adj_P_Val':diffs[8*i+6],
					'B':diffs[8*i+7],
					'link':'view in jBrowse',
					})
	os.chdir(abspath_copy)
	result = json.dumps(genes)	
	return result

@app.route('/json/meta/<tissue>/<tp>',methods=['GET','POST'])
def jsonmeta(tissue,tp):
	filepath={'DE':'MetaDE/summary.csv','QC':'MetaQC/summaryTable.csv','Path':'outputsMetaPath/summary.csv'}
	head={'DE':'SYMBOL','QC':'Series','Path':'Pathway_Annotation'}
	dic=[]
	result =[]
	with open(os.path.abspath('.')+'/static/series/meta/'+tissue+'_result/'+filepath[tp],'r') as csv_f:
		reader=csv.reader(csv_f)
		fieldnames = next(reader)
		fieldnames[0]=head[tp]
		csv_reader = csv.DictReader(csv_f,fieldnames=fieldnames)
		for row in csv_reader:
			d = {}
			for k, v in row.items():
				d[k] = v
			dic.append(d)
	result = json.dumps(dic)
	return result

@app.template_filter('diff_txt')
def difftxt(gse):
	if os.path.exists(os.path.abspath('.')+'/static/series/diff_summary/GSE'+gse+'.txt'):
		file=open(os.path.abspath('.')+'/static/series/diff_summary/GSE'+gse+'.txt')
		diffstr=file.read()
		strs=diffstr.splitlines()
	else:
		strs=''
	return(strs)

@app.template_filter('ppi_link')
def ppilink(diff_gse):
	all_ppi={'diff.AD--Older_control.GSE29654':'https://version-11-0.string-db.org/cgi/network.pl?networkId=rBB97hqZHcNr',
			'diff.AD--Younger_control.GSE29654':'https://version-11-0.string-db.org/cgi/network.pl?networkId=Y5dqnC6itBwE',
			'diff.AD--Younger_control.GSE29676':'https://version-11-0.string-db.org/cgi/network.pl?networkId=VKIqLh8h9Jw4',
			'diff.AD--Older_control.GSE29676':'https://version-11-0.string-db.org/cgi/network.pl?networkId=4N3zL59Jqddq',
			'diff.AD--control.GSE39087':'https://version-11-0.string-db.org/cgi/network.pl?networkId=wpLzt15pkmoR',
			'diff.AD--control.GSE62283':'https://version-11-0.string-db.org/cgi/network.pl?networkId=x0bF0iOfBLcw',
			'diff.AD--control.GSE74763':'https://version-11-0.string-db.org/cgi/network.pl?networkId=jA7KEjIUL4Pc'
			}
	if diff_gse in all_ppi:
		strs=all_ppi[diff_gse]
	else:
		strs=''
	return(strs)

@app.template_filter('contain_gene')
def containgene(search_input):
	db = pymysql.connect(
		host='101.132.188.227',
		port=3306,
		user='test',
		passwd='123456',
		db='ADDB')
	cur = db.cursor()
	tmp_list =search_input.split(":")
	chro = tmp_list[0]
	input_location = tmp_list[1]
	if '-' in input_location:
		input_location1=input_location.split('-')[0]
		input_location2=input_location.split('-')[1]
		search_position = "select geneID,SYMBOL from gene_location where chr = '%s' and start <= '%s' and end >= '%s' ;"  % (chro,input_location2,input_location1)
	else:
		search_position = "select geneID,SYMBOL from gene_location where chr = '%s' and start <= '%s' and end >= '%s' ;"  % (chro,input_location,input_location)
	cur.execute(search_position)
	data1 = cur.fetchall()
	db.close()
	return data1

@app.template_filter('snp_pos')
def snppos(search_input):
	db = pymysql.connect(
		host='101.132.188.227',
		port=3306,
		user='test',
		passwd='123456',
		db='ADDB')
	cur = db.cursor()
	search_rs_id = "select * from gwas where rsID='%s';" % (search_input)
	cur.execute(search_rs_id)
	table = cur.fetchall()
	db.close()
	link="http://106.14.221.173/jbrowse/"
	if table:
		link="http://106.14.221.173/jbrowse/?loc=chr"+str(table[0][1])+"%3A"+str(int(table[0][2])-150)+".."+str(int(table[0][2])+150)+"&tracks=DNA%2Cmysnps%2CGRCh37&highlight="
	return link

@app.template_filter('all_num')
def getallnum(gse):
	all_gse={'104141': 'mRNA', '104704': 'mRNA', '110226': 'mRNA', '12685': 'mRNA', '1297': 'mRNA', '13214': 'mRNA', '15222': 'mRNA', '16759': 'mRNA', '18309': 'mRNA', '26927': 'mRNA', '26972': 'mRNA', '28146': 'mRNA', '29378': 'mRNA', '29652': 'mRNA', '30945': 'mRNA', '32645': 'mRNA', '33000': 'mRNA', '34879': 'mRNA', '36980': 'mRNA', '37263': 'mRNA', '37264': 'mRNA', '39420': 'mRNA', '4226': 'mRNA', '4229': 'mRNA', '42492': 'mRNA', '43326': 'mRNA', '44768': 'mRNA', '44770': 'mRNA', '44771': 'mRNA', '45596': 'mRNA', '4757': 'mRNA', '48350': 'mRNA', '5281': 'mRNA', '53697': 'mRNA', '61196': 'mRNA', '63060': 'mRNA', '63061': 'mRNA', '6613': 'mRNA', '6834': 'mRNA', '84422_GPL570': 'mRNA', '84422_GPL96': 'mRNA', '84422_GPL97': 'mRNA', '85426': 'mRNA', '93885': 'mRNA', '95587': 'mRNA', '95810': 'mRNA', '97760': 'mRNA', '29654': 'protein', '29676': 'protein', '39087': 'protein', '62283': 'protein', '70424': 'protein', '74763': 'protein','16759_GPL8757':'miRNA','46131':'miRNA','46579':'miRNA','48552':'miRNA'}
	data=getsummary(all_gse[gse].lower())
	for row in data:
		if row[0]==gse:
			return(row[5])
	return(0)

def findpos(geneid):
	db = pymysql.connect(
				host='101.132.188.227',
				port=3306,
				user='test',
				passwd='123456',
				db='ADDB')
	cur = db.cursor()
	if geneid.isdigit():
		search_geneID = "select * from gene_location where geneID='%s';" % (geneid)
	else:
		search_geneID = "select * from gene_location where SYMBOL='%s';" % (geneid)
	cur.execute(search_geneID)
	table = cur.fetchall()
	db.close()
	poslist=[]
	if table:
		poslist.append(table[0][1])
		poslist.append(table[0][2])
		poslist.append(table[0][3])
		poslist.append(table[0][4])
		poslist.append(table[0][5])
		poslist.append(table[0][6])
		poslist.append(table[0][7])
		poslist.append(table[0][8])
		poslist.append(table[0][9])
		poslist.append(table[0][10])
		if table[0][2]<10000:
			poslist.append('http://106.14.221.173/jbrowse/?loc='+table[0][1]+'%3A0..'+str(table[0][3]+table[0][2])+'&tracks=DNA%2CGRCh37&highlight=')
		else:
			poslist.append('http://106.14.221.173/jbrowse/?loc='+table[0][1]+'%3A'+str(table[0][2]-10000)+'..'+str(table[0][3]+10000)+'&tracks=DNA%2CGRCh37&highlight=')
	return(poslist)


def rounding(decimal_num):
	if decimal_num<0.0001 or decimal_num>10000:
		return(str('%.2e' % decimal_num))
	else:
		return(str(decimal_num.quantize(Decimal('0.0000'))))

def tissue(old):
	if ',' in old:
		return((old.replace(',',', ')).replace(', ',' (',1)+')')
	else:
		return(old)

def gsecontent(gse,num,diff_num):
	all_gse={'104141': 'mRNA', '104704': 'mRNA', '110226': 'mRNA', '12685': 'mRNA', '1297': 'mRNA', '13214': 'mRNA', '15222': 'mRNA', '16759': 'mRNA', '18309': 'mRNA', '26927': 'mRNA', '26972': 'mRNA', '28146': 'mRNA', '29378': 'mRNA', '29652': 'mRNA', '30945': 'mRNA', '32645': 'mRNA', '33000': 'mRNA', '34879': 'mRNA', '36980': 'mRNA', '37263': 'mRNA', '37264': 'mRNA', '39420': 'mRNA', '4226': 'mRNA', '4229': 'mRNA', '42492': 'mRNA', '43326': 'mRNA', '44768': 'mRNA', '44770': 'mRNA', '44771': 'mRNA', '45596': 'mRNA', '4757': 'mRNA', '48350': 'mRNA', '5281': 'mRNA', '53697': 'mRNA', '61196': 'mRNA', '63060': 'mRNA', '63061': 'mRNA', '6613': 'mRNA', '6834': 'mRNA', '84422_GPL570': 'mRNA', '84422_GPL96': 'mRNA', '84422_GPL97': 'mRNA', '85426': 'mRNA', '93885': 'mRNA', '95587': 'mRNA', '95810': 'mRNA', '97760': 'mRNA', '29654': 'protein', '29676': 'protein', '39087': 'protein', '62283': 'protein', '74763': 'protein','16759_GPL8757':'miRNA','46131':'miRNA','46579':'miRNA','48552':'miRNA'}
	if gse in all_gse:
		tp=all_gse[gse]
	else:
		tp='mRNA'
	if_file=False
	if_diff=False
	if_all=False
	file_list=[]
	diff_list=[]
	all_file_list=[]
	all_list=[]
	abspath_copy=os.path.abspath('.')
	filepath=os.path.abspath('.')+'/static/series/'+tp+'/GSE'+gse
	all_filepath=os.path.abspath('.')+'/static/series/'+tp+'/ALL/GSE'+gse
	if os.path.exists(filepath):
		if_file=True
		if os.path.exists(all_filepath):
			if_all=True
		pathdir = os.listdir(filepath) 
		os.chdir(filepath)
		for fi in pathdir:
			newdir = os.path.join(filepath,fi)
			if os.path.isdir(newdir):
				if os.path.basename(newdir)[0:5] == 'diff.': 
					file_list.append(newdir+'/'+os.path.basename(newdir)+'.txt')
		for path in file_list:
			file=open(path)
			file_list[file_list.index(path)]=path.split('.')[1]
			diff=file.read()
			diffs=re.split('\t\n|\t|\n',diff)
			genes=[
				{
				'ENTREZ':'ENTREZ',
				'SYMBOL':'SYMBOL',
				'logFC':'logFC',
				'AveExpr':'AveExpr',
				't':'t',
				'P.Value':'P.Value',
				'adj.P.Val':'adj.P.Val',
				'B':'B',
				}
				]
			if tp == 'miRNA':
				genes[0]['SYMBOL']='miRNA_ID'
				if len(diffs)>10:
					for i in range(1,int((len(diffs))/7)):
						genes.append({
							'ENTREZ':'-',
							'SYMBOL':diffs[7*i-1],
							'logFC':diffs[7*i],
							'AveExpr':diffs[7*i+1],
							't':diffs[7*i+2],
							'P.Value':diffs[7*i+3],
							'adj.P.Val':diffs[7*i+4],
							'B':diffs[7*i+5],
							})
				diff_list.append(genes)
			else:
				if len(diffs)>10:
					for i in range(1,int((len(diffs))/8)):
						genes.append({
							'ENTREZ':diffs[8*i],
							'SYMBOL':diffs[8*i+1],
							'logFC':diffs[8*i+2],
							'AveExpr':diffs[8*i+3],
							't':diffs[8*i+4],
							'P.Value':diffs[8*i+5],
							'adj.P.Val':diffs[8*i+6],
							'B':diffs[8*i+7],
							})
				diff_list.append(genes)
		genes_len=0
		for genes in diff_list:
			genes_len=genes_len+len(genes)-1
		if genes_len!=0:
			if_diff=True
	os.chdir(abspath_copy)
	if '_all' in diff_num:
		all_pathdir = os.listdir(all_filepath) 
		os.chdir(all_filepath)
		for fi in all_pathdir:
			all_newdir = os.path.join(all_filepath,fi)
			all_file_list.append(all_newdir)
		for path in all_file_list:
			file=open(path)
			all_file_list[all_file_list.index(path)]=path.split('.')[1]
			diff=file.read()
			diffs=re.split('\t\n|\t|\n',diff)
			genes=[
				{
				'ENTREZ':'ENTREZ',
				'SYMBOL':'SYMBOL',
				'logFC':'logFC',
				'AveExpr':'AveExpr',
				't':'t',
				'P.Value':'P.Value',
				'adj.P.Val':'adj.P.Val',
				'B':'B',
				}
				]
			if tp == 'miRNA':
				genes[0]['SYMBOL']='miRNA_ID'
				if len(diffs)>10:
					for i in range(1,int((len(diffs))/7)):
						genes.append({
							'ENTREZ':'-',
							'SYMBOL':diffs[7*i-1],
							'logFC':diffs[7*i],
							'AveExpr':diffs[7*i+1],
							't':diffs[7*i+2],
							'P.Value':diffs[7*i+3],
							'adj.P.Val':diffs[7*i+4],
							'B':diffs[7*i+5],
							})
				all_list.append(genes)
			else:
				if len(diffs)>10:
					for i in range(1,int((len(diffs))/8)):
						genes.append({
							'ENTREZ':diffs[8*i],
							'SYMBOL':diffs[8*i+1],
							'logFC':diffs[8*i+2],
							'AveExpr':diffs[8*i+3],
							't':diffs[8*i+4],
							'P.Value':diffs[8*i+5],
							'adj.P.Val':diffs[8*i+6],
							'B':diffs[8*i+7],
							})
				all_list.append(genes)
	os.chdir(abspath_copy)
	return tp,file_list,diff_list,all_file_list,all_list,if_file,if_diff,if_all

def getsummary(tp):
	data = []
	with open(os.path.abspath('.')+'/static/series/info_csv/'+tp+'.csv','r',encoding='utf-8') as f:
		reader = csv.reader(f,dialect='excel')
		for row in reader:
			row[1]=tissue(row[1])
			data.append(row)
		del data[0]
	return data

@app.errorhandler(404)
def page_not_found(error):
    return render_template('404.html'),404

if __name__ == '__main__':
	app.run(debug=True)
