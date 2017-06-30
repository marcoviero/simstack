import pdb
import six
import numpy as np
from astropy.io import fits
import astropy.units as u
from utils import bin_ndarray as rebin
from utils import gauss_kern
from utils import clean_nans
from utils import clean_args
from astropy import cosmology
from astropy.cosmology import Planck15 as cosmo
from astropy.cosmology import Planck15, z_at_value

class Field_catalogs:
	def __init__(self, tbl, zkey='z_peak', mkey='lmass', rkey='ra', dkey='dec'):
		self.table = tbl
		self.nsrc = len(tbl)
		self.zkey = zkey#.lower()
		self.mkey = mkey#.lower()
		self.rkey = rkey#.lower()
		self.dkey = dkey#.lower()
		#pdb.set_trace()

	def separate_pops_by_index(self, cuts_dict, uvj = True):
		'''
		This is a generalized classifier of galaxies.
		cuts_dict is a dictionary with conditions,
		where the key is the index of the population (e.g., 0,1,2)
		E.g.;
		cuts_dict = {1:[criteria],3:[criteria],..,5:[criteria]}
		cuts_dict[2] = [criteria]
		E.g.;
		cuts_dict = {0:[],1:[],2:['agn',['F_ratio',50,False]],3:['sb',['lage',False,7.5]]}

		The critera contain ['key_string',greater-than,less-than] values with False
		when only one condition.  E.g.,
			cuts_dict[3] = ['dusty',[['lage', False, 7.5]]]
			cuts_dict[2] = ['agn',[['F_ratio', 40, False]]]
			cuts_dict[4] = ['sb',[['mips24',300,False],['lage', False, 7.5]]]
		Ncrit == len(cuts_dict) [-2 if UVJ used for sf/qt]
		Conditions should go in descending order... hard when a dictionary
		has no intrinsic order.  First operation is to determine the order and
		set conditions.
		'''
		sfg = np.ones(self.nsrc)
		npop = len(cuts_dict)
		if uvj == True:
			Ncrit = npop - 2
		else:
			Ncrit = npop
		for i in range(self.nsrc):
			# Go through conditions in descending order.
			# continue when one is satisfied
			for j in range(npop)[::-1][:Ncrit]:
				#name = cuts_dict[j][0]
				conditions = cuts_dict[j][1]
				ckey = conditions[0]
				if (conditions[1] == False) & (conditions[2] == False):
					if (self.table[ckey][i] == conditions[3]):
						sfg[i]= j
						continue
				elif conditions[1] == False:
					if (self.table[ckey][i] < conditions[2]):
						sfg[i] = j
						continue
				elif conditions[2] == False:
					if (self.table[ckey][i] > conditions[1]):
						sfg[i] = j
						continue
				else:
					if (self.table[ckey][i] > conditions[1]) & (self.table[ckey][i] < conditions[2]):
						sfg[i] = j
						continue
			# If no condition yet met then see if it's Quiescent
			if (sfg[i] == 1) & (uvj == True):
					if (self.table.rf_U_V.values[i] > 1.3) and (self.table.rf_V_J.values[i] < 1.5):
						if (self.table[self.zkey][i] < 1):
							if (self.table.rf_U_V.values[i] > (self.table.rf_V_J.values[i]*0.88+0.69) ): sfg[i]=0
						if (self.table[self.zkey][i] > 1):
							if (self.table.rf_U_V.values[i] > (self.table.rf_V_J.values[i]*0.88+0.59) ): sfg[i]=0

		self.table['sfg'] = sfg

	def separate_pops_by_name(self, cuts_dict):
		'''
		This is a generalized classifier of galaxies.
		cuts_dict is a dictionary with conditions,
		where the key is the population (e.g., 'sf', 'agn','sf0')
		E.g.;
		cuts_dict['agn'] = [pop_index,[criteria]]
		The critera contain ['key',greater-than,less-than] values with False
		when only one condition.  E.g.,
			cuts_dict['dusty'] = [3,[['lage', False, 7.5]]]
			cuts_dict['agn'] = [2,[['F_ratio', 40, False]]]
			cuts_dict['sb'] = [4,[['mips24',300,False],['lage', False, 7.5]]]

		Ncrit == len(cuts_dict) [-2 if UVJ used for sf/qt]
		Conditions should go in descending order... hard when a dictionary
		has no intrinsic order.  First operation is to determine the order and
		set conditions.
		'''
		sfg = np.ones(self.nsrc)
		npop = len(cuts_dict)
		if 'qt' in cuts_dict:
			Ncrit = npop - 2
			uvj = True
		else:
			Ncrit = npop
			uvj = False
		print Ncrit
		#Set (descending) order of cuts.
		#names      = [k for k in cuts_dict][::-1]
		#reverse_ind here is the arguments indices in reversed order
		reverse_ind = np.argsort([cuts_dict[k][0] for k in cuts_dict])[::-1]
		ind         = [cuts_dict[k][0] for k in cuts_dict]
		conditions  = [cuts_dict[k][1] for k in cuts_dict]
		for i in range(self.nsrc):
			# Go through conditions in descending order.
			# continue when one is satisfied
			for j in range(Ncrit):
				icut = reverse_ind[j]
				ckey = conditions[icut][0]
				if (conditions[icut][1] == False) & (conditions[icut][2] == False):
					if (self.table[ckey].values[i] == conditions[icut][3]):
						sfg[i]=ind[icut]
						continue
				elif conditions[icut][1] == False:
					if (self.table[ckey].values[i] < conditions[icut][2]):
						sfg[i]=ind[icut]
						continue
				elif conditions[icut][2] == False:
					if (self.table[ckey].values[i] > conditions[icut][1]):
						sfg[i]=ind[icut]
						continue
				else:
					if (self.table[ckey].values[i] > conditions[icut][1]) & (self.table[ckey].values[i] < conditions[icut][2]):
						sfg[i]=ind[icut]
						continue
			# If no condition yet met then see if it's Quiescent
			if (sfg[i] == 1) & (uvj == True):
					if (self.table.rf_U_V.values[i] > 1.3) and (self.table.rf_V_J.values[i] < 1.5):
						if (self.table[self.zkey][i] < 1):
							if (self.table.rf_U_V.values[i] > (self.table.rf_V_J.values[i]*0.88+0.69) ): sfg[i]=0
						if (self.table[self.zkey][i] > 1):
							if (self.table.rf_U_V.values[i] > (self.table.rf_V_J.values[i]*0.88+0.59) ): sfg[i]=0

		self.table['sfg'] = sfg

	def separate_sf_qt(self):
		sfg = np.ones(self.nsrc)
		#pdb.set_trace()
		for i in range(self.nsrc):
			if (self.table.rf_U_V.values[i] > 1.3) and (self.table.rf_V_J.values[i] < 1.5):
				if (self.table[self.zkey][i] < 1):
					if (self.table.rf_U_V.values[i] > (self.table.rf_V_J.values[i]*0.88+0.69) ): sfg[i]=0
				if (self.table[self.zkey][i] > 1):
					if (self.table.rf_U_V.values[i] > (self.table.rf_V_J.values[i]*0.88+0.59) ): sfg[i]=0
		#indsf = np.where(sfg == 1)
		#indqt = np.where(sfg == 0)
		self.table['sfg'] = sfg

	def separate_additional_feature(self, feature_name):
		return 'what am i doing?'

	def separate_sf_qt_agn(self,Fcut = 20):
		sfg = np.ones(self.nsrc)
		#pdb.set_trace()
		#AGN
		for i in range(self.nsrc):
			if (self.table.F_ratio[i] >= Fcut):
				sfg[i]=2
			else:
				if (self.table.rf_U_V.values[i] > 1.3) and (self.table.rf_V_J.values[i] < 1.5):
					if (self.table[self.zkey][i] < 1):
						if (self.table.rf_U_V.values[i] > (self.table.rf_V_J.values[i]*0.88+0.69) ): sfg[i]=0
					if (self.table[self.zkey][i] > 1):
						if (self.table.rf_U_V.values[i] > (self.table.rf_V_J.values[i]*0.88+0.59) ): sfg[i]=0
		#indsf = np.where(sfg == 1)
		#indqt = np.where(sfg == 0)
		self.table['sfg'] = sfg

	def separate_sf_qt_agn_Ahat(self,Fcut = 20, Ahat = 0.5):
		sfg = np.ones(self.nsrc)
		#pdb.set_trace()
		#AGN
		for i in range(self.nsrc):
			if (self.table.F_ratio[i] >= Fcut) & (self.table.a_hat_AGN[i] >= Ahat):
				sfg[i]=2
			else:
				if (self.table.rf_U_V.values[i] > 1.3) and (self.table.rf_V_J.values[i] < 1.5):
					if (self.table[self.zkey][i] < 1):
						if (self.table.rf_U_V.values[i] > (self.table.rf_V_J.values[i]*0.88+0.69) ): sfg[i]=0
					if (self.table[self.zkey][i] > 1):
						if (self.table.rf_U_V.values[i] > (self.table.rf_V_J.values[i]*0.88+0.59) ): sfg[i]=0
		#indsf = np.where(sfg == 1)
		#indqt = np.where(sfg == 0)
		self.table['sfg'] = sfg
	#def separate_nuv_r(self):

	def separate_uvj_pops(self, uvj_nodes = np.array([0.0,0.6,1.0,1.5,2.0,2.5])):
		sfg = np.ones(self.nsrc)
		ncolor= len(uvj_nodes)-1 + 1 #+1 is qt
		for i in range(self.nsrc):
			if (self.table.rf_U_V.values[i] > 1.3) and (self.table.rf_V_J.values[i] < 1.5):
				if (self.table[self.zkey][i] < 1):
					if (self.table.rf_U_V.values[i] > (self.table.rf_V_J.values[i]*0.88+0.69) ): sfg[i]=ncolor-1
				if (self.table[self.zkey][i] > 1):
					if (self.table.rf_U_V.values[i] > (self.table.rf_V_J.values[i]*0.88+0.59) ): sfg[i]=ncolor-1
			if sfg[i] < ncolor-1:
				for j in range(ncolor-1):
					if ((self.table.UVJ.values[i]) >= uvj_nodes[j]) & ((self.table.UVJ.values[i]) < uvj_nodes[j+1]): sfg[i]=j
		self.table['sfg'] = sfg

	def separate_6pops(self, Fcut = 40, tau_cut = 7.5, age_cut = 7.4):
		sfg = np.ones(self.nsrc)
		z0 = 9.3
		slope = (12.0 - 9.3 ) / (1.25)
		for i in range(self.nsrc):
			if (self.table.lage.values[i] < age_cut):
				sfg[i]=4 #SB
			elif (self.table.ltau.values[i] < tau_cut):
				sfg[i]=3 #Dusty
			elif (self.table[self.mkey][i] > self.table[self.zkey][i] * slope + z0) & (self.table.mips24[i] > 250):
				sfg[i]=5 #LOCAL
			elif (self.table.F_ratio[i] >= Fcut):
				sfg[i]=2 #AGN
			else:
				if (self.table.rf_U_V.values[i] > 1.3) and (self.table.rf_V_J.values[i] < 1.5):
					if (self.table[self.zkey][i] < 1):
						if (self.table.rf_U_V.values[i] > (self.table.rf_V_J.values[i]*0.88+0.69) ): sfg[i]=0
					if (self.table[self.zkey][i] > 1):
						if (self.table.rf_U_V.values[i] > (self.table.rf_V_J.values[i]*0.88+0.59) ): sfg[i]=0
		self.table['sfg'] = sfg

	def separate_5pops(self, Fcut = 25, MIPS24_cut = 200.0):
		sfg = np.ones(self.nsrc)
		z0 = 9.3
		slope = (12.0 - 9.3 ) / (1.25)
		for i in range(self.nsrc):
			#pdb.set_trace()
			if (self.table.mips24.values[i] >= (MIPS24_cut+100.)) & (self.table.F_ratio.values[i] < Fcut):
				sfg[i]=4
			elif (self.table[self.mkey][i] > self.table[self.zkey][i] * slope + z0) & (self.table.mips24.values[i] > MIPS24_cut):
				sfg[i]=3
			elif (self.table.F_ratio.values[i] >= Fcut):
				sfg[i]=2
			else:
				if (self.table.rf_U_V.values[i] > 1.3) and (self.table.rf_V_J.values[i] < 1.5):
					if (self.table[self.zkey][i] < 1):
						if (self.table.rf_U_V.values[i] > (self.table.rf_V_J.values[i]*0.88+0.69) ): sfg[i]=0
					if (self.table[self.zkey][i] > 1):
						if (self.table.rf_U_V.values[i] > (self.table.rf_V_J.values[i]*0.88+0.59) ): sfg[i]=0
		self.table['sfg'] = sfg

	def separate_5popsB(self, Fcut = 50, ssfr = 30):
		sfg = np.ones(self.nsrc)
		z0 = 9.3
		slope = (12.0 - 9.3 ) / (1.25)
		for i in range(self.nsrc):
			if (10**(self.table.lssfr.values[i])*1e9 >= ssfr ) & (self.table.F_ratio.values[i] < Fcut):
				sfg[i]=4
			elif (self.table[self.mkey][i] > self.table[self.zkey][i] * slope + z0) & (self.table.mips24.values[i] > 250):
				sfg[i]=3
			elif (self.table.F_ratio.values[i] >= Fcut):
				sfg[i]=2
			else:
				if (self.table.rf_U_V.values[i] > 1.3) and (self.table.rf_V_J.values[i] < 1.5):
					if (self.table[self.zkey][i] < 1):
						if (self.table.rf_U_V.values[i] > (self.table.rf_V_J.values[i]*0.88+0.69) ): sfg[i]=0
					if (self.table[self.zkey][i] > 1):
						if (self.table.rf_U_V.values[i] > (self.table.rf_V_J.values[i]*0.88+0.59) ): sfg[i]=0
		self.table['sfg'] = sfg

	def separate_4pops(self, Fcut = 40, age_cut = 7.5):
		sfg = np.ones(self.nsrc)
		for i in range(self.nsrc):
			if (self.table.lage.values[i] <= age_cut):
				sfg[i]=3
			elif (self.table.F_ratio.values[i] >= Fcut):
				sfg[i]=2
			else:
				if (self.table.rf_U_V.values[i] > 1.3) and (self.table.rf_V_J.values[i] < 1.5):
					if (self.table[self.zkey][i] < 1):
						if (self.table.rf_U_V.values[i] > (self.table.rf_V_J.values[i]*0.88+0.69) ): sfg[i]=0
					if (self.table[self.zkey][i] > 1):
						if (self.table.rf_U_V.values[i] > (self.table.rf_V_J.values[i]*0.88+0.59) ): sfg[i]=0
		self.table['sfg'] = sfg

	def separate_3pops(self,  Fcut = 40):
		sfg = np.ones(self.nsrc)
		for i in range(self.nsrc):
			if (self.table.F_ratio.values[i] >= Fcut):
				sfg[i]=2
			else:
				if (self.table.rf_U_V.values[i] > 1.3) and (self.table.rf_V_J.values[i] < 1.5):
					if (self.table[self.zkey][i] < 1):
						if (self.table.rf_U_V.values[i] > (self.table.rf_V_J.values[i]*0.88+0.69) ): sfg[i]=0
					if (self.table[self.zkey][i] > 1):
						if (self.table.rf_U_V.values[i] > (self.table.rf_V_J.values[i]*0.88+0.59) ): sfg[i]=0
		self.table['sfg'] = sfg

	def separate_3pops_sb(self,  age_cut = 7.4):
		sfg = np.ones(self.nsrc)
		for i in range(self.nsrc):
			if (self.table.lage.values[i] <= age_cut):
				sfg[i]=2
			else:
				if (self.table.rf_U_V.values[i] > 1.3) and (self.table.rf_V_J.values[i] < 1.5):
					if (self.table[self.zkey][i] < 1):
						if (self.table.rf_U_V.values[i] > (self.table.rf_V_J.values[i]*0.88+0.69) ): sfg[i]=0
					if (self.table[self.zkey][i] > 1):
						if (self.table.rf_U_V.values[i] > (self.table.rf_V_J.values[i]*0.88+0.59) ): sfg[i]=0
		self.table['sfg'] = sfg

	def get_subpop_ids(self, znodes, mnodes, pop_dict, linear_mass=1, lookback_time = False):
		self.subpop_ids = {}
		if lookback_time == True:
			age_universe = cosmo.age(0).value # 13.797617455819209 Gyr
			znodes = np.array([z_at_value(cosmo.age,(age_universe - i) * u.Gyr) for i in znodes])

		for iz in range(len(znodes[:-1])):
			for jm in range(len(mnodes[:-1])):
				for k in pop_dict:
					if linear_mass == 1:
						ind_mz =( (self.table.sfg.values == pop_dict[k][0]) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
							(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )
					else:
						ind_mz =( (self.table.sfg == pop_dict[k][0]) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
							(self.table[self.mkey] >= np.min(mnodes[jm:jm+2])) & (self.table[self.mkey] < np.max(mnodes[jm:jm+2])) )

					self.subpop_ids['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_'+k] = self.table.ID[ind_mz].values
					#self.subpop_ids['z_'+clean_args('{:.3f}'.format(znodes[iz]))+'_'+clean_args('{:.3f}'.format(znodes[iz+1]))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_'+k] = self.table.ID[ind_mz].values

	def get_sf_qt_mass_redshift_bins(self, znodes, mnodes):
		self.id_z_ms = {}
		for iz in range(len(znodes[:-1])):
			for jm in range(len(mnodes[:-1])):
				ind_mz_sf =( (self.table.sfg == 1) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
					(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )
				ind_mz_qt =( (self.table.sfg == 0) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
					(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )

				self.id_z_ms['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_sf'] = self.table.ID[ind_mz_sf].values
				self.id_z_ms['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_qt'] = self.table.ID[ind_mz_qt].values

	def get_sf_qt_mass_lookback_time_bins(self, tnodes, mnodes):
		self.id_lookt_mass = {}
		age_universe = cosmo.age(0).value # 13.797617455819209 Gyr
		znodes = np.array([z_at_value(cosmo.age,(age_universe - i) * u.Gyr) for i in tnodes])

		for iz in range(len(znodes[:-1])):
			for jm in range(len(mnodes[:-1])):
				ind_mt_sf =( (self.table.sfg == 1) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
					(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )
				ind_mt_qt =( (self.table.sfg == 0) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
					(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )

				self.id_lookt_mass['lookt_'+clean_args(str('{:.2f}'.format(tnodes[iz])))+'_'+clean_args(str('{:.2f}'.format(tnodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_sf'] = self.table.ID[ind_mt_sf].values
				self.id_lookt_mass['lookt_'+clean_args(str('{:.2f}'.format(tnodes[iz])))+'_'+clean_args(str('{:.2f}'.format(tnodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_qt'] = self.table.ID[ind_mt_qt].values

	def get_sf_qt_agn_mass_redshift_bins(self, znodes, mnodes):
		self.id_z_ms = {}
		for iz in range(len(znodes[:-1])):
			for jm in range(len(mnodes[:-1])):
				ind_mz_sf =( (self.table.sfg == 1) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
					(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )
				ind_mz_qt =( (self.table.sfg == 0) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
					(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )
				ind_mz_agn =( (self.table.sfg == 2) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
					(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )

				self.id_z_ms['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_sf'] = self.table.ID[ind_mz_sf].values
				self.id_z_ms['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_qt'] = self.table.ID[ind_mz_qt].values
				self.id_z_ms['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_agn'] = self.table.ID[ind_mz_agn].values

	def get_mass_redshift_uvj_bins(self, znodes, mnodes, cnodes, linear_mass=1):
		#pop_suf = ['qt','sf0','sf1','sf2','sf3']
		self.id_z_ms_pop = {}
		for iz in range(len(znodes[:-1])):
			for jm in range(len(mnodes[:-1])):
				for kc in range(len(cnodes)):
					if linear_mass == 1:
						ind_mz =( (self.table.sfg == kc) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
							(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )
					else:
						ind_mz =( (self.table.sfg == kc) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
							(self.table[self.mkey] >= np.min(mnodes[jm:jm+2])) & (self.table[self.mkey] < np.max(mnodes[jm:jm+2])) )

					if kc < len(cnodes)-1: csuf = 'sf'+str(kc)
					else: csuf = 'qt'
					self.id_z_ms_pop['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_'+csuf] = self.table.ID[ind_mz].values

	def get_3pops_mass_redshift_bins(self, znodes, mnodes, linear_mass=1):
		#pop_suf = ['sf','qt','agn']
		self.id_z_ms_3pop = {}
		for iz in range(len(znodes[:-1])):
			for jm in range(len(mnodes[:-1])):
				if linear_mass == 1:
					ind_mz_sf =( (self.table.sfg == 1) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )
					ind_mz_qt =( (self.table.sfg == 0) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )
					ind_mz_agn =( (self.table.sfg == 2) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )
				else:
					ind_mz_sf =( (self.table.sfg == 1) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(self.table[self.mkey] >= np.min(mnodes[jm:jm+2])) & (self.table[self.mkey] < np.max(mnodes[jm:jm+2])) )
					ind_mz_qt =( (self.table.sfg == 0) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(self.table[self.mkey] >= np.min(mnodes[jm:jm+2])) & (self.table[self.mkey] < np.max(mnodes[jm:jm+2])) )
					ind_mz_agn =( (self.table.sfg == 2) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(self.table[self.mkey] >= np.min(mnodes[jm:jm+2])) & (self.table[self.mkey] < np.max(mnodes[jm:jm+2])) )

				self.id_z_ms_3pop['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_sf'] = self.table.ID[ind_mz_sf].values
				self.id_z_ms_3pop['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_qt'] = self.table.ID[ind_mz_qt].values
				self.id_z_ms_3pop['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_agn'] = self.table.ID[ind_mz_agn].values

	def get_3pops_mass_redshift_bins_sb(self, znodes, mnodes, linear_mass=1):
		#pop_suf = ['sf','qt','sb']
		self.id_z_ms_3pop = {}
		for iz in range(len(znodes[:-1])):
			for jm in range(len(mnodes[:-1])):
				if linear_mass == 1:
					ind_mz_sf =( (self.table.sfg == 1) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )
					ind_mz_qt =( (self.table.sfg == 0) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )
					ind_mz_sb =( (self.table.sfg == 2) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )
				else:
					ind_mz_sf =( (self.table.sfg == 1) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(self.table[self.mkey] >= np.min(mnodes[jm:jm+2])) & (self.table[self.mkey] < np.max(mnodes[jm:jm+2])) )
					ind_mz_qt =( (self.table.sfg == 0) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(self.table[self.mkey] >= np.min(mnodes[jm:jm+2])) & (self.table[self.mkey] < np.max(mnodes[jm:jm+2])) )
					ind_mz_sb =( (self.table.sfg == 2) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(self.table[self.mkey] >= np.min(mnodes[jm:jm+2])) & (self.table[self.mkey] < np.max(mnodes[jm:jm+2])) )

				self.id_z_ms_3pop['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_sf'] = self.table.ID[ind_mz_sf].values
				self.id_z_ms_3pop['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_qt'] = self.table.ID[ind_mz_qt].values
				self.id_z_ms_3pop['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_sb'] = self.table.ID[ind_mz_sb].values

	def get_4pops_mass_redshift_bins(self, znodes, mnodes, linear_mass=1):
		#pop_suf = ['sf','qt','agn','sb']
		self.id_z_ms_4pop = {}
		for iz in range(len(znodes[:-1])):
			for jm in range(len(mnodes[:-1])):
				if linear_mass == 1:
					ind_mz_sf =( (self.table.sfg == 1) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )
					ind_mz_qt =( (self.table.sfg == 0) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )
					ind_mz_agn =( (self.table.sfg == 2) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )
					ind_mz_sb =( (self.table.sfg == 3) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )
				else:
					ind_mz_sf =( (self.table.sfg == 1) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(self.table[self.mkey] >= np.min(mnodes[jm:jm+2])) & (self.table[self.mkey] < np.max(mnodes[jm:jm+2])) )
					ind_mz_qt =( (self.table.sfg == 0) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(self.table[self.mkey] >= np.min(mnodes[jm:jm+2])) & (self.table[self.mkey] < np.max(mnodes[jm:jm+2])) )
					ind_mz_agn =( (self.table.sfg == 2) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(self.table[self.mkey] >= np.min(mnodes[jm:jm+2])) & (self.table[self.mkey] < np.max(mnodes[jm:jm+2])) )
					ind_mz_sb =( (self.table.sfg == 3) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(self.table[self.mkey] >= np.min(mnodes[jm:jm+2])) & (self.table[self.mkey] < np.max(mnodes[jm:jm+2])) )

				self.id_z_ms_4pop['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_sf'] = self.table.ID[ind_mz_sf].values
				self.id_z_ms_4pop['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_qt'] = self.table.ID[ind_mz_qt].values
				self.id_z_ms_4pop['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_agn'] = self.table.ID[ind_mz_agn].values
				self.id_z_ms_4pop['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_sb'] = self.table.ID[ind_mz_sb].values

	def get_5pops_mass_redshift_bins(self, znodes, mnodes, linear_mass=1):
		#pop_suf = ['sf','qt','agn','sb','loc']
		self.id_z_ms_5pop = {}
		for iz in range(len(znodes[:-1])):
			for jm in range(len(mnodes[:-1])):
				if linear_mass == 1:
					ind_mz_sf =( (self.table.sfg == 1) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )
					ind_mz_qt =( (self.table.sfg == 0) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )
					ind_mz_agn =( (self.table.sfg == 2) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )
					ind_mz_sb =( (self.table.sfg == 4) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )
					ind_mz_loc=( (self.table.sfg == 3) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )
				else:
					ind_mz_sf =( (self.table.sfg == 1) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(self.table[self.mkey] >= np.min(mnodes[jm:jm+2])) & (self.table[self.mkey] < np.max(mnodes[jm:jm+2])) )
					ind_mz_qt =( (self.table.sfg == 0) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(self.table[self.mkey] >= np.min(mnodes[jm:jm+2])) & (self.table[self.mkey] < np.max(mnodes[jm:jm+2])) )
					ind_mz_agn =( (self.table.sfg == 2) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(self.table[self.mkey] >= np.min(mnodes[jm:jm+2])) & (self.table[self.mkey] < np.max(mnodes[jm:jm+2])) )
					ind_mz_sb =( (self.table.sfg == 4) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(self.table[self.mkey] >= np.min(mnodes[jm:jm+2])) & (self.table[self.mkey] < np.max(mnodes[jm:jm+2])) )
					ind_mz_loc=( (self.table.sfg == 3) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(self.table[self.mkey] >= np.min(mnodes[jm:jm+2])) & (self.table[self.mkey] < np.max(mnodes[jm:jm+2])) )

				self.id_z_ms_5pop['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_sf'] = self.table.ID[ind_mz_sf].values
				self.id_z_ms_5pop['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_qt'] = self.table.ID[ind_mz_qt].values
				self.id_z_ms_5pop['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_agn'] = self.table.ID[ind_mz_agn].values
				self.id_z_ms_5pop['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_sb'] = self.table.ID[ind_mz_sb].values
				self.id_z_ms_5pop['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_loc'] = self.table.ID[ind_mz_loc].values

	def get_6pops_mass_redshift_bins(self, znodes, mnodes, linear_mass=1):
		#pop_suf = ['qt',sf','agn','dst',sb','loc']
		self.id_z_ms_6pop = {}
		for iz in range(len(znodes[:-1])):
			for jm in range(len(mnodes[:-1])):
				if linear_mass == 1:
					ind_mz_sf =( (self.table.sfg == 1) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )
					ind_mz_qt =( (self.table.sfg == 0) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )
					ind_mz_agn =( (self.table.sfg == 2) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )
					ind_mz_sb =( (self.table.sfg == 4) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )
					ind_mz_dst=( (self.table.sfg == 3) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )
					ind_mz_loc=( (self.table.sfg == 5) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )
				else:
					ind_mz_sf =( (self.table.sfg == 1) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(self.table[self.mkey] >= np.min(mnodes[jm:jm+2])) & (self.table[self.mkey] < np.max(mnodes[jm:jm+2])) )
					ind_mz_qt =( (self.table.sfg == 0) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(self.table[self.mkey] >= np.min(mnodes[jm:jm+2])) & (self.table[self.mkey] < np.max(mnodes[jm:jm+2])) )
					ind_mz_agn =( (self.table.sfg == 2) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(self.table[self.mkey] >= np.min(mnodes[jm:jm+2])) & (self.table[self.mkey] < np.max(mnodes[jm:jm+2])) )
					ind_mz_sb =( (self.table.sfg == 4) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(self.table[self.mkey] >= np.min(mnodes[jm:jm+2])) & (self.table[self.mkey] < np.max(mnodes[jm:jm+2])) )
					ind_mz_dst=( (self.table.sfg == 3) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(self.table[self.mkey] >= np.min(mnodes[jm:jm+2])) & (self.table[self.mkey] < np.max(mnodes[jm:jm+2])) )
					ind_mz_loc=( (self.table.sfg == 5) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(self.table[self.mkey] >= np.min(mnodes[jm:jm+2])) & (self.table[self.mkey] < np.max(mnodes[jm:jm+2])) )

				self.id_z_ms_6pop['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_sf'] = self.table.ID[ind_mz_sf].values
				self.id_z_ms_6pop['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_qt'] = self.table.ID[ind_mz_qt].values
				self.id_z_ms_6pop['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_agn'] = self.table.ID[ind_mz_agn].values
				self.id_z_ms_6pop['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_dst'] = self.table.ID[ind_mz_dst].values
				self.id_z_ms_6pop['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_sb'] = self.table.ID[ind_mz_sb].values
				self.id_z_ms_6pop['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+'_loc'] = self.table.ID[ind_mz_loc].values

	def get_general_redshift_bins(self, znodes, mnodes, sfg = 1, suffx = '', Fcut = 25, Ahat = 1.0, initialize_pop = False):
		if initialize_pop == True: self.id_z_ms = {}
		for iz in range(len(znodes[:-1])):
			for jm in range(len(mnodes[:-1])):
				ind_mz =( (self.table.sfg == 1) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
					(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )

				self.id_z_ms['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+suffx] = self.table.ID[ind_mz].values

	def get_mass_redshift_bins(self, znodes, mnodes, sfg = 1, pop_suffix = '', initialize_pop = False):
		if initialize_pop == True: self.id_z_ms_pop = {}
		for iz in range(len(znodes[:-1])):
			for jm in range(len(mnodes[:-1])):
				ind_mz =( (self.table.sfg == sfg) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
					(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) )

				self.id_z_ms_pop['z_'+clean_args(str('{:.2f}'.format(znodes[iz])))+'_'+clean_args(str('{:.2f}'.format(znodes[iz+1])))+'__m_'+clean_args(str('{:.2f}'.format(mnodes[jm])))+'_'+clean_args(str('{:.2f}'.format(mnodes[jm+1])))+pop_suffix] = self.table.ID[ind_mz].values

	def get_criteria_specific_redshift_bins(self, znodes, mnodes, sfg = 1, criteria = '', crange = [1.0], initialize_pop = False):
		pop = ['qt','sf']
		nc = len(crange)
		if initialize_pop == True: self.id_crit = {}

		for iz in range(len(znodes[:-1])):
			for jm in range(len(mnodes[:-1])):
				if nc > 1:
					for kc in range(len(crange[:-1])):
						ind_crit =( (self.table.sfg == sfg) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
							(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) &
							(clean_nans(self.table[criteria]) >= crange[kc]) & (clean_nans(self.table[criteria]) < crange[kc+1]) )

						arg = 'z_'+str('{:.2f}'.format(znodes[iz]))+'_'+str('{:.2f}'.format(znodes[iz+1]))+'__m_'+str('{:.2f}'.format(mnodes[jm]))+'_'+str('{:.2f}'.format(mnodes[jm+1]))+'__'+criteria+'_'+str('{:.2f}'.format(crange[kc],2))+'_'+str('{:.2f}'.format(crange[kc+1],2))+'_'+pop[sfg]
						self.id_crit[clean_args(arg)] = self.table.ID[ind_crit].values
				else:
					#above and below no?
					ind_above =( (self.table.sfg == sfg) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) &
						(clean_nans(self.table[criteria]) >= crange[0]) )
					ind_below =( (self.table.sfg == sfg) & (self.table[self.zkey] >= np.min(znodes[iz:iz+2])) & (self.table[self.zkey] < np.max(znodes[iz:iz+2])) &
						(10**self.table[self.mkey] >= 10**np.min(mnodes[jm:jm+2])) & (10**self.table[self.mkey] < 10**np.max(mnodes[jm:jm+2])) &
						(clean_nans(self.table[criteria]) < crange[0]) )

					arg = 'z_'+str('{:.2f}'.format(znodes[iz]))+'_'+str('{:.2f}'.format(znodes[iz+1]))+'__m_'+str('{:.2f}'.format(mnodes[jm]))+'_'+str('{:.2f}'.format(mnodes[jm+1]))+'__'+criteria+'_ge_'+str('{:.2f}'.format(crange[0],2))+'_'+pop[sfg]
					self.id_crit[clean_args(arg)] = self.table.ID[ind_above].values
					arg = 'z_'+str('{:.2f}'.format(znodes[iz]))+'_'+str('{:.2f}'.format(znodes[iz+1]))+'__m_'+str('{:.2f}'.format(mnodes[jm]))+'_'+str('{:.2f}'.format(mnodes[jm+1]))+'__'+criteria+'_lt_'+str('{:.2f}'.format(crange[0],2))+'_'+pop[sfg]
					self.id_crit[clean_args(arg)] = self.table.ID[ind_below].values


	def get_parent_child_redshift_bins(self,znodes):
		self.id_z_sed = {}
		for ch in self.table.parent.unique():
			for iz in range(len(znodes[:-1])):
				self.id_z_sed['z_'+clean_args(str(znodes[iz]))+'_'+clean_args(str(znodes[iz+1]))+'__sed'+str(ch)] = self.table.ID[ (self.table.parent == ch) & (self.table[self.zkey] >= znodes[iz]) & (self.table[self.zkey] < znodes[iz+1]) ].values

	def get_parent_child_bins(self):
		self.id_children = {}
		for ch in self.table.parent.unique():
			self.id_children['sed'+str(ch)] = self.table.ID[self.table.parent == ch].values

	def subset_positions(self,radec_ids):
		''' This positions function is very general.
			User supplies IDs dictionary, function returns RA/DEC dictionaries with the same keys'''
		#self.ra_dec = {}
		ra_dec = {}
		#ra = {}
		#dec = {}
		#pdb.set_trace()
		for k in radec_ids.keys():
			ra0  = self.table[self.rkey]
			dec0 = self.table[self.dkey]
			ra  = ra0[self.table.ID.isin(radec_ids[k])].values
			dec = dec0[self.table.ID.isin(radec_ids[k])].values
			#ra  = self.table.ra[self.table.ID.isin(radec_ids[k])].values
			#dec = self.table.dec[self.table.ID.isin(radec_ids[k])].values
			ra_dec[k] = [ra,dec]
		return ra_dec
