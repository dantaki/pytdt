#!python
#cython: boundscheck=False, wraparound=False
import numpy as np
import pandas as pd
import scipy.stats
import statsmodels.api as sm
import math,sys
from tqdm import tqdm
cdef binom_calc(t):
	pat_p = scipy.stats.binom_test(t[1],t[0]+t[1])
	mat_p = scipy.stats.binom_test(t[3],t[2]+t[3])
	tot_p = scipy.stats.binom_test(t[1]+t[3],sum(t))
	return tuple((pat_p,mat_p,tot_p))
cdef binom_confid(t):
	pat_ci,mat_ci,tot_ci = (0,0),(0,0),(0,0)
	if t[0]+t[1]>0:
		pat_ci = tuple(sm.stats.proportion_confint(t[1],t[0]+t[1]))
	if t[2]+t[3]>0:
		mat_ci = tuple(sm.stats.proportion_confint(t[3],t[2]+t[3]))
	if sum(t)>0:
		tot_ci = tuple(sm.stats.proportion_confint(t[1]+t[3],sum(t)))
	return tuple((pat_ci,mat_ci,tot_ci))
cdef chisq_calc(t):
	pat_chisq,mat_chisq,tot_chisq= 'nan','nan','nan'
	if t[0]+t[1]!=0: pat_chisq = ((t[0]-t[1])**2)/float((t[0]+t[1]))
	if t[2]+t[3]!=0: mat_chisq = ((t[2]-t[3])**2)/float((t[2]+t[3]))
	if all(x==0 for x in t)==False: tot_chisq = (((t[1]+t[3])-(t[0]+t[2]))**2)/float((t[0]+t[1]+t[2]+t[3]))
	return pat_chisq,mat_chisq,tot_chisq
cdef chisq_contingency(t):
	term=0
	if 0 in t[0] or 0 in t[1]: term=1
	if term==1: return np.nan	
	else:  return scipy.stats.chi2_contingency(t)[1]
cdef confid_inter(t,aff,con):
	log_odds, std_err= ['nan','nan','nan'],['nan','nan','nan']
	for x in range(0,len(t)):
		i = t[x]
		if i!='inf' and i!=0: log_odds[x]=np.log(i)
	# paternal 
	if aff[0]!=0 and aff[1]!=0 and con[0]!=0 and con[1]!=0: std_err[0]=std_err_calc(aff[0],aff[1],con[0],con[1])
	# maternal 
	if aff[2]!=0 and aff[3]!=0 and con[2]!=0 and con[3]!=0: std_err[1]=std_err_calc(aff[2],aff[3],con[2],con[3])
	tot = (aff[0]+aff[2],aff[1]+aff[3],con[0]+con[2],con[1]+con[3])
	if tot[0]!=0 and tot[1]!=0 and tot[2]!=0 and tot[3]!=0: std_err[2]=std_err_calc(tot[0],tot[1],tot[2],tot[3])
	ci = [ ('nan','nan'),('nan','nan'),('nan','nan') ]
	for x in range(0,len(std_err)):
		i,j= log_odds[x],std_err[x]
		if i!='nan' and j!='nan': 
			ll = i - (1.96*j)
			hl = i + (1.96*j)
			ci[x]=(np.exp(ll),np.exp(hl))
	return tuple(ci)
def group_by(df,col,bins):
	"""
	group_by returns a binned group dict

	@df: pandas dataframe object
	columns, id, 'col'

	@col: string of a column in df to bin variants

	@bins: list of values to bin
	"""
	group={}
	for x in bins:
		group[x]=list(df.loc[df[col] >= x].id)
	return group
cdef odds_calc(t):
	# o = trans / non-trans
	pat_odds, mat_odds, tot_odds = 0,0,0
	if t[0] != 0: pat_odds = t[1]/float(t[0])
	if t[2] != 0: mat_odds = t[3]/float(t[2])
	if t[0]+t[2] !=0: tot_odds = (t[1]+t[3])/float((t[0]+t[2]))
	return pat_odds,mat_odds,tot_odds
cdef odds_diff(rats,con_ratio,t1,t2):
	difs=[]
	for x in range(0,len(rats)):
		i,j = rats[x],con_ratio[x]
		if i!='inf' and j!='inf': difs.append(i-j)
		else: difs.append('nan')
	return tuple(difs)
cdef odds_ratio(t1,t2):
	rats=[]
	for x in range(0,len(t1)):
		i,j = t1[x],t2[x]
		if j==0: rats.append('inf')
		else: rats.append(i/j)
	return tuple(rats)
cdef poo_table(t): return [[t[1],t[0]],[t[3],t[2]]]
cdef pval(t):
	p = ['nan','nan','nan']
	for x in range(0,len(t)):
		if t[x]!='nan': p[x]=scipy.stats.distributions.chi2.sf(t[x],1)
	return tuple(p)
cdef std_err_calc(a,b,c,d):
	return math.sqrt(((1/float(a)) + (1/float(b)) + (1/float(c)) + (1/float(d))))
cdef table(a,c):
	pat = [[a[1],a[0]],[c[1],c[0]]]
	mat = [[a[3],a[2]],[c[3],c[2]]]
	tot = [[a[1]+a[3],a[0]+a[2]],[c[1]+c[3],c[0]+c[2]]]
	return pat,mat,tot
cdef trans_rate(t):
	p_rate, m_rate, t_rate = 'nan','nan','nan'
	if t[0]+t[1] > 0: p_rate = t[1]/ float(t[0]+t[1])
	if t[2]+t[3] > 0: m_rate = t[3]/float(t[2]+t[3])
	if sum(t) > 0:  t_rate = float(t[3]+t[1])/ float(sum(t))
	return p_rate,m_rate,t_rate
"""
calculate stats for each group
"""
def manual_tdt(aff_tdt=None,con_tdt=None):
	stats = [ Tdt(aff_tdt,con_tdt) ]
	return stats
def tdt(group=None,df=None,tot=False,prog=False):
	"""
	tdt returns an object containing 
	group id, tdt stats
	@ Group
	dict of [group id ] = [list of snp ids]
	
	@df accepts a pandas data frame
	columns id,aff_pt, aff_pn, aff_mt, aff_mn, con_pt, con_pn, con_mt, con_mn

	@tot when True, only process total transmissions (affected+control)
	"""
	stats={}
	for group_id in tqdm(group,disable= not prog,file=sys.stdout):
		ids = group[group_id]
		grp_snps = df.loc[df['id'].isin(ids)]
		if len(grp_snps)==0:
			stats[group_id]=-1
		if tot == False:
			aff_tdt = tuple(map(sum,(grp_snps.aff_pn,grp_snps.aff_pt,grp_snps.aff_mn,grp_snps.aff_mt)))
			con_tdt = tuple(map(sum,(grp_snps.con_pn,grp_snps.con_pt,grp_snps.con_mn,grp_snps.con_mt)))
			stats[group_id] = Tdt(aff_tdt,con_tdt)
		if tot==True:
			tot_tdt = tuple(map(sum,((grp_snps.aff_pn+grp_snps.con_pn),(grp_snps.aff_pt+grp_snps.con_pt),(grp_snps.aff_mn+grp_snps.con_mn),(grp_snps.aff_mt+grp_snps.con_mt))))
			stats[group_id]= Tot_Tdt(tot_tdt)
	return stats
class Tdt():
	def __init__(self,aff_tdt=None,con_tdt=None):
		# calc odds ratio
		self.aff_tdt,self.con_tdt = aff_tdt,con_tdt
		self.aff_prate,self.aff_mrate,self.aff_trate = trans_rate(self.aff_tdt)
		self.con_prate,self.con_mrate,self.con_trate = trans_rate(self.con_tdt)
		self.aff_odds,self.con_odds = odds_calc(aff_tdt), odds_calc(con_tdt)
		self.odds_ratio = odds_ratio(self.aff_odds,self.con_odds)
		self.odds_ci = confid_inter(self.odds_ratio,aff_tdt,con_tdt)
		self.aff_pval,self.con_pval = binom_calc(self.aff_tdt),binom_calc(self.con_tdt)
		self.aff_ci, self.con_ci = binom_confid(self.aff_tdt),binom_confid(self.con_tdt)
		#self.con_odds_ratio = odds_ratio(self.con_odds,self.aff_odds)
		#self.odds_diff = odds_diff(self.odds_ratio,self.con_odds_ratio,self.aff_odds,self.con_odds)
		#self.con_odds_diff = tuple([x *-1 for x in self.odds_diff])
		#self.con_odds_ci = confid_inter(self.con_odds_ratio,aff_tdt,con_tdt)
		#self.aff_chisq,self.con_chisq = chisq_calc(aff_tdt),chisq_calc(con_tdt)
		#self.aff_pval , self.con_pval = pval(self.aff_chisq),pval(self.con_chisq)
		#self.aff_poo_odds, self.aff_poo_pval = scipy.stats.fisher_exact(poo_table(self.aff_tdt))
		#self.con_poo_odds, self.con_poo_pval = scipy.stats.fisher_exact(poo_table(self.con_tdt))
		pat_tab,mat_tab,tot_tab = table(self.aff_tdt,self.con_tdt)
		self.pat_chisq, self.mat_chisq, self.tot_chisq = chisq_contingency(pat_tab),chisq_contingency(mat_tab),chisq_contingency(tot_tab)
		self.pat_fisher_odds, self.pat_fisher_pval = scipy.stats.fisher_exact(pat_tab)
		self.mat_fisher_odds, self.mat_fisher_pval = scipy.stats.fisher_exact(mat_tab)
		self.tot_fisher_odds, self.tot_fisher_pval = scipy.stats.fisher_exact(tot_tab)
	def print_stats(self,group_id=None):
		print("{:^84}\nGroup ID: {}\n{:^84}".format('-'*84,group_id,'-'*84))
		print("{:^20} {:^20} {:^20} {:^20}".format('PHENOTYPE','PATERNAL T:N','MATERNAL T:N','TOTAL T:N'))
		print("{:^84}".format('-'*84))
		print("{:^20} {:^20} {:^20} {:^20}".format('affected','{}:{}'.format(self.aff_tdt[1],self.aff_tdt[0]),'{}:{}'.format(self.aff_tdt[3],self.aff_tdt[2]),'{}:{}'.format(self.aff_tdt[1]+self.aff_tdt[3],self.aff_tdt[0]+self.aff_tdt[2])))
		print("{:^20} {:^20} {:^20} {:^20}".format('control','{}:{}'.format(self.con_tdt[1],self.con_tdt[0]),'{}:{}'.format(self.con_tdt[3],self.con_tdt[2]),'{}:{}'.format(self.con_tdt[1]+self.con_tdt[3],self.con_tdt[0]+self.con_tdt[2])))
		print("{:^20} {:^20} {:^20} {:^20}".format('affected','{}'.format(self.aff_prate),'{}'.format(self.aff_mrate),'{}'.format(self.aff_trate)))
		print("{:^20} {:^20} {:^20} {:^20}".format('control','{}'.format(self.con_prate),'{}'.format(self.con_mrate),'{}'.format(self.con_trate)))
		print("{:^20} {:^20} {:^20} {:^20}".format('.','paternal','maternal','total'))
		print("{:^84}".format('-'*84))
		print("{:^20} {:^20.2f} {:^20.2f} {:^20.2f}".format('ODDS affected',*list(map(lambda y:0 if type(y) is str else y, self.aff_odds))))
		print("{:^20} {:^20.2f} {:^20.2f} {:^20.2f}".format('ODDS control',*list(map(lambda y:0 if type(y) is str else y, self.con_odds))))
		print("{:^20} {:^20.2f} {:^20.2f} {:^20.2f}".format('ODDS RATIO',*list(map(lambda y:0 if type(y) is str else y, self.odds_ratio))))
		print("{:^20} {:^20.2f} {:^20.2f} {:^20.2f}".format('ODDS DIFFERENCE',*list(map(lambda y:0 if type(y) is str else y, self.odds_diff))))
		aff_chisq,aff_pval,con_chisq,con_pval = list(map(lambda y:0 if type(y) is str else y, self.aff_chisq)),list(map(lambda y:0 if type(y) is str else y, self.aff_pval)),list(map(lambda y:0 if type(y) is str else y, self.con_chisq)),list(map(lambda y:0 if type(y) is str else y, self.con_pval))
		print("{:^20} {:^10.2f}:{:^10.2f} {:^10.2f}:{:^10.2f} {:^10.2f}:{:^10.2f}".format('CHISQ:pval affected',aff_chisq[0],aff_pval[0],aff_chisq[1],aff_pval[1],aff_chisq[2],aff_pval[2]))
		print("{:^20} {:^10.2f}:{:^10.2f} {:^10.2f}:{:^10.2f} {:^10.2f}:{:^10.2f}".format('CHISQ:pval control',con_chisq[0],con_pval[0],con_chisq[1],con_pval[1],con_chisq[2],con_pval[2]))
		print("{:^84}".format('-'*84))
class Tot_Tdt():
	def __init__(self,tot_tdt=None):
		self.tdt= tot_tdt
		self.prate,self.mrate,self.trate = trans_rate(self.tdt)
		#self.odds = odds_calc(self.tdt)
		#self.chisq = chisq_calc(self.tdt)
		#self.pval = pval(self.chisq)
		#self.fisher_odds, self.fisher_pval = scipy.stats.fisher_exact(poo_table(self.tdt))
