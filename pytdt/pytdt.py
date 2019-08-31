#!/usr/bin/env python3
import numpy as np
import pytdt.pytdtStats as pytdtStats
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
	return pytdtStats.tdt(group,df,tot,prog)
def manual_tdt(aff_tdt=None,con_tdt=None): return pytdtStats.manual_tdt(aff_tdt,con_tdt)