import numpy as np
import pandas as pd
import scanpy as sc
from scanpy import AnnData
import os
from scipy import stats
import anndata
import json
from itertools import combinations
import psutil
from random import sample

def intersection(lst1, lst2): 
	lst3 = [value for value in lst1 if value in lst2] 
	return lst3 

def unique(list1): 
	# insert the list to the set 
	list_set = set(list1) 
	# convert the set to the list 
	unique_list = (list(list_set)) 
	return unique_list

def my_reset(*varnames):
	#varnames are what you want to keep
	globals_ = globals()
	to_save = {v: globals_[v] for v in varnames}
	to_save['my_reset'] = my_reset  # lets keep this function by default
	del globals_
	get_ipython().run_line_magic("reset", '-f')
	globals().update(to_save)

	

#create gene x metacell matrix where each metacell is a region
def create_metacell(adata,metacell_axis,metacells,min_cells,min_expression):
    genes = list(adata.var_names)
    metacell_df = pd.DataFrame(genes)
    metacell_log_df = pd.DataFrame(genes)
    for mcell in metacells:
        adata_temp=adata[adata.obs[metacell_axis] == mcell].copy()
        adata_temp.X = adata_temp.layers['raw']
        
        ncells = adata_temp.n_obs

        #check that each region has at least min_cells
        if ncells >= min_cells:
            #collapse data into meta-cells by summing counts across all cells in region
            metacell = np.sum(adata_temp.X,axis=0).T
            totalcounts = np.sum(metacell)
            metacell = metacell/totalcounts * 1e5 #normalize to 100,000 counts total
            metacell_log = np.log10(metacell) #take the log

            #from Fenna's code - very low values may not be meaningful, so change them to zero
            metacell_log[np.isnan(metacell_log)] = 0
            metacell_log[metacell_log == -np.inf] = 0
            metacell_log[metacell_log == np.inf] = 0
            metacell_log[metacell_log < 0] = 0
            
            metacell_df[mcell] = metacell
            metacell_log_df[mcell] = metacell_log
            
            del metacell
        else:
            print('Insufficient number of cells in ' + mcell +'!')
            print(ncells)

        del adata_temp
        
    metacell_df.set_index(0,drop=True,inplace=True)
    metacell_log_df.set_index(0,drop=True,inplace=True)
    
    #impose minimum expression limit with OR logic - at least one region needs 10 per 100,000
    metacell_log_df = metacell_log_df[(metacell_df>=min_expression).any(axis=1)]
    metacell_df = metacell_df[(metacell_df>=min_expression).any(axis=1)]

    return metacell_df,metacell_log_df
    
    
def calc_logFC(metacell_log_df,mcell1,mcell2):
    col1=metacell_log_df[mcell1]
    col2=metacell_log_df[mcell2]

    FC = col1-col2 #calculate fold change

    FC = pd.DataFrame({'gene': FC.index, 'change':FC.values})
    
    FC['higher_group'] = 'tbd'
    FC['lower_group'] = 'tbd'

    for i in FC.index:
        if FC['change'][i]>0:
            FC['higher_group'][i]= mcell1
            FC['lower_group'][i] = mcell2
        else:
            FC['higher_group'][i]= mcell2
            FC['lower_group'][i] = mcell1   
    
    return FC



def find_DEGs_merge_reps_pw(adata,metacell_axis,metacells,DEG_thresh,min_cells,min_expression,n_random_downsample, 
                         replicate_key, min_fraction):
	reps = unique(list(adata.obs[replicate_key]))

	#loop through each replicate
	for idx,rep in enumerate(reps):

		if n_random_downsample > 0:
			rep_barcodes = list(adata[adata.obs[replicate_key]==rep].obs_names)
			random_barcodes = sample(rep_barcodes,n_random_downsample)
			adata_temp = adata[random_barcodes,:].copy()
		else:
			adata_temp = adata[adata.obs[replicate_key]==rep].copy()

		metacell_temp,log_metacell_temp=create_metacell(adata=adata_temp,
						              metacell_axis=metacell_axis,
						              metacells=metacells,
						              min_cells=min_cells,
						              min_expression=min_expression)

		mcells = []
		for i in unique(adata_temp.obs[metacell_axis]):
			if len(adata_temp[adata_temp.obs[metacell_axis] == i].obs.index)>50:
				mcells.append(i)

		#initialize lists
		DEGs_all_rep = []
		higher_all_rep = []
		lower_all_rep = []
		lfc_all_rep = []
		ngenes=[]
		comps=[]
		nDEGs=[]
		replist=[]

		replist_long=[]
		complist_long=[]

		combos = combinations(mcells,2)
		#for some reason this function gives a weird data structure so convert it to something useful
		region_combos = []
		for c in combos: region_combos.append(c)
		combo_num = list(range(0, len(region_combos)))

		#first make a dictionary of genes that meet the min inclusion threshold so we don't ask this over and over
		inclusion_dict = {}
		for group in mcells:
			adata_group = adata_temp[adata_temp.obs[metacell_axis]== group].copy()
			min_inclusion = len(adata_group.obs.index)*min_fraction
			sc.pp.filter_genes(adata_group, min_cells = min_inclusion)
			group_inclusion_genes = list(adata_group.var.index) 
			inclusion_dict[group] = group_inclusion_genes

		#loop through each pairwise metacell comparison
		for c in combo_num:

			mcell1 = region_combos[c][0]
			mcell2 = region_combos[c][1]

			print('Ranking pairing ' + str(c + 1) + ' of '+ str(len(combo_num))
			+', ' + str(mcell1) + ' and ' + str(mcell2))

			comp_name = mcell1 + '-' + mcell2
			#print(comp_name)

			FC = calc_logFC(log_metacell_temp,mcell1,mcell2)

			#take only DEGs over the threshold
			FC = FC[np.abs(FC['change'])>=DEG_thresh]

			#only take genes above the min fraction included
			group_inclusion_genes = inclusion_dict[mcell1]
			rest_inclusion_genes = inclusion_dict[mcell2]
			inclusion_genes = group_inclusion_genes + rest_inclusion_genes
			inclusion_genes = unique(inclusion_genes)

			FC = FC[FC['gene'].isin(inclusion_genes)]
			DEG_names_rep = list(FC['gene'])
			DEGs_all_rep.extend(DEG_names_rep)

			#record which region is higher
			higher_rep = FC['higher_group']
			higher_names_rep = list(higher_rep)
			higher_all_rep.extend(higher_names_rep)

			#record which region is lower
			lower_rep = FC['lower_group']
			lower_names_rep = list(lower_rep)
			lower_all_rep.extend(lower_names_rep)
			
			#record LFC
			this_lfc = FC['change'].abs()
			this_lfc.astype('float')
			#absolute value since we specify higher and lower group
			lfc_rep = list(this_lfc)
			lfc_all_rep.extend(lfc_rep)


			#longform replicate name list, for matching up with DEGs
			replist_long.extend([rep]*len(DEG_names_rep))

			#longform comparison name list, for matcching up with DEGs
			complist_long.extend([comp_name]*len(DEG_names_rep))

			#how many DEGs?
			nDEGs.append(len(DEG_names_rep))

			#how many genes were over minimum expression threshold?
			ngenes.append(len(list(FC.index))) 

			#collect the comparison name
			comps.append(comp_name)

			#collect the replicate name
			replist.append(rep)

			#zip up the outputs and create a dataframe
			rep_condensed_output = list(zip(replist,comps,ngenes,nDEGs))
			rep_long_output = list(zip(replist_long,complist_long,higher_all_rep,lower_all_rep,DEGs_all_rep,lfc_all_rep))

			#convert to dataframe format
			condensed_DEG_df_rep = pd.DataFrame(rep_condensed_output,
								        columns=['replicate','comparison','ngenes','nDEGs'])
			long_DEG_df_rep = pd.DataFrame(rep_long_output,
								        columns=['replicate','comparison','higher_group','lower_group','DEG_name', 'logFC'])       

		if idx == 0: #if in first replicate
			condensed_DEG_df = condensed_DEG_df_rep
			long_DEG_df = long_DEG_df_rep
		else:
			condensed_DEG_df = pd.concat([condensed_DEG_df,condensed_DEG_df_rep],axis=0)
			long_DEG_df = pd.concat([long_DEG_df,long_DEG_df_rep],axis=0)

	return condensed_DEG_df,long_DEG_df


def find_DEGs_merge_reps_1vr(adata,metacell_axis,metacells,DEG_thresh,min_cells,min_expression,n_random_downsample, 
                         replicate_key, min_fraction):
	reps = unique(list(adata.obs[replicate_key]))

	#loop through each replicate
	for idx,rep in enumerate(reps):

		if n_random_downsample > 0:
			rep_barcodes = list(adata[adata.obs[replicate_key]==rep].obs_names)
			random_barcodes = sample(rep_barcodes,n_random_downsample)
			adata_temp = adata[random_barcodes,:].copy()
		else:
			adata_temp = adata[adata.obs[replicate_key]==rep].copy()


		mcells = []
		for i in unique(adata_temp.obs[metacell_axis]):
			if len(adata_temp[adata_temp.obs[metacell_axis] == i].obs.index)>50:
				mcells.append(i)

		#initialize lists
		DEGs_all_rep = []
		higher_all_rep = []
		lower_all_rep = []
		ngenes=[]
		comps=[]
		nDEGs=[]
		replist=[]
		lfc_all_rep = []

		replist_long=[]
		complist_long=[]


		#loop through each metacell axis option
		for m in mcells:
			onevr_ID = []
			DEG_names_rep = []
			group_inclusion_genes = []
			rest_inclusion_genes = []

			for i in adata_temp.obs.index:
				if adata_temp.obs[metacell_axis][i] == m:
					onevr_ID.append(m)
				else:
					onevr_ID.append('rest')
			adata_temp.obs['onevr'] = onevr_ID


			metacell_temp,log_metacell_temp=create_metacell(adata=adata_temp,
							    metacell_axis= 'onevr',
							    metacells=[m, 'rest'],
							    min_cells=min_cells,
							    min_expression=min_expression)
			mcell1 = m
			mcell2 = 'rest'


			comp_name = mcell1 + '-' + mcell2
			#print(comp_name)

			FC = calc_logFC(log_metacell_temp,mcell1,mcell2)

			#take only DEGs over the threshold
			candidates = list(FC[np.abs(FC['change'])>=DEG_thresh]['gene'])

			adata_group = adata_temp[adata_temp.obs.onevr == m].copy()
			min_inclusion = len(adata_group.obs.index)*min_fraction
			sc.pp.filter_genes(adata_group, min_cells = min_inclusion)
			group_inclusion_genes = list(adata_group.var.index)

			adata_rest = adata_temp[adata_temp.obs.onevr == 'rest'].copy()
			min_inclusion = len(adata_rest.obs.index)*min_fraction
			sc.pp.filter_genes(adata_rest, min_cells = min_inclusion)
			rest_inclusion_genes = list(adata_rest.var.index)

			inclusion_genes = group_inclusion_genes + rest_inclusion_genes

			inclusion_genes = unique(inclusion_genes)

			DEG_names_rep = intersection(candidates, inclusion_genes)

			DEGs_all_rep.extend(DEG_names_rep)

			#record which region is higher
			higher_names_rep = []
			for g in DEG_names_rep:
					higher = str(FC[FC['gene'] == g]['higher_group'])
					if 'rest' in higher:
						  higher_names_rep.append('rest')
					else:
						  higher_names_rep.append(m)
			higher_all_rep.extend(higher_names_rep)

			#record which region is lower
			lower_names_rep = []
			for g in DEG_names_rep:
				lower = str(FC[FC['gene'] == g]['lower_group'])
				if 'rest' in lower:
					lower_names_rep.append('rest')
				else:
					lower_names_rep.append(m)
			lower_all_rep.extend(lower_names_rep)
			
			#record LFC
			lfc_rep = []
			for g in DEG_names_rep:
				this_lfc = float(FC[FC['gene'] == g]['change'])
				this_lfc = np.abs(this_lfc)
				lfc_rep.append(this_lfc)
			lfc_all_rep.extend(lfc_rep)


			#longform replicate name list, for matching up with DEGs
			replist_long.extend([rep]*len(DEG_names_rep))

			#longform comparison name list, for matcching up with DEGs
			complist_long.extend([comp_name]*len(DEG_names_rep))

			#how many DEGs?
			nDEGs.append(len(DEG_names_rep))

			#how many genes were over minimum expression threshold?
			ngenes.append(len(list(FC.index))) 

			#collect the comparison name
			comps.append(comp_name)

			#collect the replicate name
			replist.append(rep)

		#zip up the outputs and create a dataframe
		rep_condensed_output = list(zip(replist,comps,ngenes,nDEGs))
		rep_long_output = list(zip(replist_long,complist_long,higher_all_rep,lower_all_rep,DEGs_all_rep,lfc_all_rep))

		#convert to dataframe format
		condensed_DEG_df_rep = pd.DataFrame(rep_condensed_output,
					                      columns=['replicate','comparison','ngenes','nDEGs'])
		long_DEG_df_rep = pd.DataFrame(rep_long_output,
					                      columns=['replicate','comparison','higher_group','lower_group','DEG_name', 'logFC'])       

		if idx == 0: #if in first replicate
			condensed_DEG_df = condensed_DEG_df_rep
			long_DEG_df = long_DEG_df_rep
		else:
			condensed_DEG_df = pd.concat([condensed_DEG_df,condensed_DEG_df_rep],axis=0)
			long_DEG_df = pd.concat([long_DEG_df,long_DEG_df_rep],axis=0)

	return condensed_DEG_df,long_DEG_df

    
    
def pairwise_clustering_analysis_metacell(adata, cluster_res, lfc, inclusion, replicate_key, min_expression, min_cells):
	pairwise_DEGs = 0
	oneVrest_DEGs = 0
	cluster_num = 0
	
	print("Cluster resolution = " + str(cluster_res))
	sc.tl.leiden(adata,resolution=cluster_res)
	sc.pl.umap(adata,color='leiden', palette = 'tab20')
	
	#get names of clusters and get combos of clusters
	clusters = unique(adata.obs.leiden)
	cluster_num = len(clusters)

	   
	condensed_DEG_df_pw,long_DEG_df_pw = find_DEGs_merge_reps_pw(adata=adata, 
															metacell_axis='leiden',
															metacells=clusters,
															DEG_thresh=lfc,
															min_cells=min_cells, 
															min_expression=min_expression,
                                                    		n_random_downsample=0, 
                                                    		replicate_key = replicate_key, 
                                                            min_fraction = inclusion)
    
    #get comparison counts
	combo_DEGs = list(condensed_DEG_df_pw['nDEGs'])
	pairwise_DEGs = min(combo_DEGs)
	
	print("Pairwise counts = ")
	print(combo_DEGs)
	print("At " + str(cluster_res) + ' resolution, this minimum DEGs of any pairing:')
	print(pairwise_DEGs)
	
	    
	condensed_DEG_df_1vr,long_DEG_df_1vr = find_DEGs_merge_reps_1vr(adata=adata,
																	metacell_axis= 'leiden',
                                                    				metacells=clusters,
                                                    				DEG_thresh=lfc,
                                                    				min_cells=min_cells,
                                                    				min_expression=min_expression,
                                                    				n_random_downsample=0,
                                                   					replicate_key = replicate_key, 
                                                                    min_fraction = inclusion)
	
	filt_counts = list(condensed_DEG_df_1vr['nDEGs'])                      
	oneVrest_DEGs = min(filt_counts)
	
	print("One vs rest counts = ") 
	print(filt_counts)
	print("At " + str(cluster_res) + ' resolution, this minimum DEGs of any one vs rest cluster:')
	print(min(filt_counts))
		
	return cluster_num, pairwise_DEGs, oneVrest_DEGs

#DM 4/12/24
def fix_comparison_names(df, **kwargs): #kwargs: order_list
    #takes in pandas dataframe and outputs the comparison argument with uniform order
    #order is either custom from kwarg input or default to alphabetical
    
    df = df.reset_index(drop = True)
    order_list = kwargs.get('order_list', None)

    if order_list == None:
        for i in df.index:
            combo = df['comparison'][i]
            split = combo.split('-')
            split.sort()
            new_combo = '-'.join(split)
            df['comparison'][i] = new_combo
    
    else:    
        def custom_order(order_list):    
            keys = {x: i for i, x in enumerate(order_list)}
            def key_function(x):
                return keys[x]
            return key_function
        order = custom_order(order_list)
        for i in df.index:
            combo = df['comparison'][i]
            split = combo.split('-')
            split.sort(key = order)
            new_combo = '-'.join(split)
            df['comparison'][i] = new_combo
    return df
    
    
#DM 9/16/24
def all_gene_DE(adata,metacell_axis,metacells,min_cells,min_expression,n_random_downsample, 
                         replicate_key, min_fraction):
	reps = unique(list(adata.obs[replicate_key]))

	#loop through each replicate
	for idx,rep in enumerate(reps):

		if n_random_downsample > 0:
			rep_barcodes = list(adata[adata.obs[replicate_key]==rep].obs_names)
			random_barcodes = sample(rep_barcodes,n_random_downsample)
			adata_temp = adata[random_barcodes,:].copy()
		else:
			adata_temp = adata[adata.obs[replicate_key]==rep].copy()

		metacell_temp,log_metacell_temp=create_metacell(adata=adata_temp,
						              metacell_axis=metacell_axis,
						              metacells=metacells,
						              min_cells=min_cells,
						              min_expression=min_expression)

		mcells = []
		for i in unique(adata_temp.obs[metacell_axis]):
			if len(adata_temp[adata_temp.obs[metacell_axis] == i].obs.index)>50:
				mcells.append(i)
				
		#first make a dictionary of genes that meet the min inclusion threshold so we don't ask this over and over
		all_candidates = []
		for group in mcells:
			adata_group = adata_temp[adata_temp.obs[metacell_axis]== group].copy()
			min_inclusion = len(adata_group.obs.index)*min_fraction
			sc.pp.filter_genes(adata_group, min_cells = min_inclusion)
			group_inclusion_genes = list(adata_group.var.index) 
			all_candidates.extend(group_inclusion_genes)
		all_candidates = unique(all_candidates)

		#initialize lists
		DEGs_all_rep = []
		higher_all_rep = []
		lower_all_rep = []
		lfc_all_rep = []
		replist_long=[]
		complist_long=[]

		combos = combinations(mcells,2)
		#for some reason this function gives a weird data structure so convert it to something useful
		region_combos = []
		for c in combos: region_combos.append(c)
		combo_num = list(range(0, len(region_combos)))

		#loop through each pairwise metacell comparison
		for c in combo_num:

			mcell1 = region_combos[c][0]
			mcell2 = region_combos[c][1]

			print('Ranking pairing ' + str(c + 1) + ' of '+ str(len(combo_num))
			+', ' + str(mcell1) + ' and ' + str(mcell2))

			comp_name = mcell1 + '-' + mcell2

			FC = calc_logFC(log_metacell_temp,mcell1,mcell2)

			FC = FC[FC['gene'].isin(all_candidates)]
			DEG_names_rep = list(FC['gene'])
			DEGs_all_rep.extend(DEG_names_rep)

			#record which region is higher
			higher_rep = FC['higher_group']
			higher_names_rep = list(higher_rep)
			higher_all_rep.extend(higher_names_rep)

			#record which region is lower
			lower_rep = FC['lower_group']
			lower_names_rep = list(lower_rep)
			lower_all_rep.extend(lower_names_rep)
			
			#record LFC
			this_lfc = FC['change'].abs()
			this_lfc.astype('float')
			#absolute value since we specify higher and lower group
			lfc_rep = list(this_lfc)
			lfc_all_rep.extend(lfc_rep)


			#longform replicate name list, for matching up with DEGs
			replist_long.extend([rep]*len(DEG_names_rep))

			#longform comparison name list, for matcching up with DEGs
			complist_long.extend([comp_name]*len(DEG_names_rep))


			#zip up the outputs and create a dataframe
			rep_long_output = list(zip(replist_long,complist_long,higher_all_rep,lower_all_rep,DEGs_all_rep,lfc_all_rep))

			#convert to dataframe format
			long_DEG_df_rep = pd.DataFrame(rep_long_output,
								        columns=['replicate','comparison','higher_group','lower_group','DEG_name', 'logFC'])       

		if idx == 0: #if in first replicate
			long_DEG_df = long_DEG_df_rep
		else:
			long_DEG_df = pd.concat([long_DEG_df,long_DEG_df_rep],axis=0)

	return long_DEG_df
	
#DM 2/4/25
	
def find_DEGs_merge_reps_fraction_pw(adata,metacell_axis,metacells,DEG_thresh,min_cells,min_expression,n_random_downsample, 
                         replicate_key, min_fraction, balanced_reps):
    reps = unique(list(adata.obs[replicate_key]))

    if balanced_reps == True: 
        #loop through each replicate
        for rep in enumerate(reps):
        
            if n_random_downsample > 0:
                rep_barcodes = list(adata[adata.obs[replicate_key]==rep].obs_names)
                random_barcodes = sample(rep_barcodes,n_random_downsample)
                adata_temp = adata[random_barcodes,:].copy()
            else:
                adata_temp = adata[adata.obs[replicate_key]==rep].copy()
        
            metacell_temp,log_metacell_temp=create_metacell(adata=adata_temp,
                                          metacell_axis=metacell_axis,
                                          metacells=metacells,
                                          min_cells=min_cells,
                                          min_expression=min_expression)
        
            mcells = []
            for i in unique(adata_temp.obs[metacell_axis]):
                if len(adata_temp[adata_temp.obs[metacell_axis] == i].obs.index)>50:
                    mcells.append(i)
        
            #initialize lists
            DEGs_all_rep = []
            higher_all_rep = []
            lower_all_rep = []
            lfc_all_rep = []
            ngenes=[]
            comps=[]
            nDEGs=[]
            replist=[]
            high_fraction_list = []
            low_fraction_list = []
            
            replist_long=[]
            complist_long=[]
        
            combos = combinations(mcells,2)
            #for some reason this function gives a weird data structure so convert it to something useful
            region_combos = []
            for c in combos: region_combos.append(c)
            combo_num = list(range(0, len(region_combos)))
        
            #first make a dictionary of genes that meet the min inclusion threshold so we don't ask this over and over
            inclusion_dict = {}
            for group in mcells:
                adata_group = adata_temp[adata_temp.obs[metacell_axis]== group].copy()
                min_inclusion = len(adata_group.obs.index)*min_fraction
                sc.pp.filter_genes(adata_group, min_cells = min_inclusion)
                group_inclusion_genes = list(adata_group.var.index) 
                inclusion_dict[group] = group_inclusion_genes
        
            #loop through each pairwise metacell comparison
            for c in combo_num:
        
                mcell1 = region_combos[c][0]
                mcell2 = region_combos[c][1]
        
                print('Ranking pairing ' + str(c + 1) + ' of '+ str(len(combo_num))
                +', ' + str(mcell1) + ' and ' + str(mcell2))
        
                comp_name = mcell1 + '-' + mcell2
                #print(comp_name)
        
                FC = calc_logFC(log_metacell_temp,mcell1,mcell2)
        
                #take only DEGs over the threshold
                FC = FC[np.abs(FC['change'])>=DEG_thresh]
        
                #only take genes above the min fraction included
                group_inclusion_genes = inclusion_dict[mcell1]
                rest_inclusion_genes = inclusion_dict[mcell2]
                inclusion_genes = group_inclusion_genes + rest_inclusion_genes
                inclusion_genes = unique(inclusion_genes)
        
                FC = FC[FC['gene'].isin(inclusion_genes)]
                DEG_names_rep = list(FC['gene'])
                DEGs_all_rep.extend(DEG_names_rep)
        
                #record which region is higher
                higher_rep = FC['higher_group']
                higher_names_rep = list(higher_rep)
                higher_all_rep.extend(higher_names_rep)
        
                #record which region is lower
                lower_rep = FC['lower_group']
                lower_names_rep = list(lower_rep)
                lower_all_rep.extend(lower_names_rep)
                
                #record LFC
                this_lfc = FC['change'].abs()
                this_lfc.astype('float')
                #absolute value since we specify higher and lower group
                lfc_rep = list(this_lfc)
                lfc_all_rep.extend(lfc_rep)
        
                #longform replicate name list, for matching up with DEGs
                replist_long.extend([rep]*len(DEG_names_rep))
        
                #longform comparison name list, for matcching up with DEGs
                complist_long.extend([comp_name]*len(DEG_names_rep))
        
                #how many DEGs?
                nDEGs.append(len(DEG_names_rep))
        
                #how many genes were over minimum expression threshold?
                ngenes.append(len(list(FC.index))) 
        
                #collect the comparison name
                comps.append(comp_name)
        
                #collect the replicate name
                replist.append(rep)
        
                #zip up the outputs and create a dataframe
                rep_condensed_output = list(zip(replist,comps,ngenes,nDEGs))
                rep_long_output = list(zip(replist_long,complist_long,higher_all_rep,lower_all_rep,DEGs_all_rep,lfc_all_rep))

        
                #convert to dataframe format
                condensed_DEG_df_rep = pd.DataFrame(rep_condensed_output,
                                            columns=['replicate','comparison','ngenes','nDEGs'])
                long_DEG_df_rep = pd.DataFrame(rep_long_output,
                                            columns=['replicate','comparison','higher_group','lower_group','DEG_name', 'logFC'])  
                
        
            if idx == 0: #if in first replicate
                condensed_DEG_df = condensed_DEG_df_rep
                long_DEG_df = long_DEG_df_rep
            else:
                condensed_DEG_df = pd.concat([condensed_DEG_df,condensed_DEG_df_rep],axis=0)
                long_DEG_df = pd.concat([long_DEG_df,long_DEG_df_rep],axis=0)

        
        return condensed_DEG_df,long_DEG_df

    elif balanced_reps == False: 
        if n_random_downsample > 0:
                rep_barcodes = list(adata.obs_names)
                random_barcodes = sample(rep_barcodes,n_random_downsample)
                adata_temp = adata[random_barcodes,:].copy()
        else:
                adata_temp = adata.copy()
        
        metacell_temp,log_metacell_temp=create_metacell(adata=adata_temp,
                                      metacell_axis=metacell_axis,
                                      metacells=metacells,
                                      min_cells=min_cells,
                                      min_expression=min_expression)
    
        mcells = []
        for i in unique(adata_temp.obs[metacell_axis]):
            if len(adata_temp[adata_temp.obs[metacell_axis] == i].obs.index)>50:
                mcells.append(i)
    
        #initialize lists
        DEGs_all_rep = []
        higher_all_rep = []
        lower_all_rep = []
        lfc_all_rep = []
        ngenes=[]
        comps=[]
        nDEGs=[]
        replist=[]
        high_fraction_list = []
        low_fraction_list = []
    
        replist_long=[]
        complist_long=[]
    
        combos = combinations(mcells,2)
        #for some reason this function gives a weird data structure so convert it to something useful
        region_combos = []
        for c in combos: region_combos.append(c)
        combo_num = list(range(0, len(region_combos)))
    
        #first make a dictionary of genes that meet the min inclusion threshold so we don't ask this over and over
        inclusion_dict = {}
        for group in mcells:
            adata_group = adata_temp[adata_temp.obs[metacell_axis]== group].copy()
            min_inclusion = len(adata_group.obs.index)*min_fraction
            sc.pp.filter_genes(adata_group, min_cells = min_inclusion)
            group_inclusion_genes = list(adata_group.var.index) 
            inclusion_dict[group] = group_inclusion_genes
    
        #loop through each pairwise metacell comparison
        for c in combo_num:
    
            mcell1 = region_combos[c][0]
            mcell2 = region_combos[c][1]
    
            print('Ranking pairing ' + str(c + 1) + ' of '+ str(len(combo_num))
            +', ' + str(mcell1) + ' and ' + str(mcell2))
    
            comp_name = mcell1 + '-' + mcell2
            #print(comp_name)
    
            FC = calc_logFC(log_metacell_temp,mcell1,mcell2)
    
            #take only DEGs over the threshold
            FC = FC[np.abs(FC['change'])>=DEG_thresh]
    
            #only take genes above the min fraction included
            group_inclusion_genes = inclusion_dict[mcell1]
            rest_inclusion_genes = inclusion_dict[mcell2]
            inclusion_genes = group_inclusion_genes + rest_inclusion_genes
            inclusion_genes = unique(inclusion_genes)
    
            FC = FC[FC['gene'].isin(inclusion_genes)]
            DEG_names_rep = list(FC['gene'])
            DEGs_all_rep.extend(DEG_names_rep)
    
            #record which region is higher
            higher_rep = FC['higher_group']
            higher_names_rep = list(higher_rep)
            higher_all_rep.extend(higher_names_rep)
    
            #record which region is lower
            lower_rep = FC['lower_group']
            lower_names_rep = list(lower_rep)
            lower_all_rep.extend(lower_names_rep)
            
            #record LFC
            this_lfc = FC['change'].abs()
            this_lfc.astype('float')
            #absolute value since we specify higher and lower group
            lfc_rep = list(this_lfc)
            lfc_all_rep.extend(lfc_rep)
    
            #longform replicate name list, for matching up with DEGs
            replist_long.extend('m'*len(DEG_names_rep))
    
            #longform comparison name list, for matcching up with DEGs
            complist_long.extend([comp_name]*len(DEG_names_rep))
    
            #how many DEGs?
            nDEGs.append(len(DEG_names_rep))
    
            #how many genes were over minimum expression threshold?
            ngenes.append(len(list(FC.index))) 
    
            #collect the comparison name
            comps.append(comp_name)
    
            #collect the replicate name
            replist.append("merged")

            #how many reps is this gene in
            for gene in DEG_names_rep:
                higher_group = FC[FC['gene']==gene]['higher_group']
                higher_group = higher_group.iloc[0]
                high = adata_temp[adata_temp.obs[metacell_axis]== higher_group]
                high_total = len(unique(high.obs[replicate_key]))
                high_reps = unique(high.obs[replicate_key])
                high_33 = 0
                for i in high_reps:
                    high_rep = high[high.obs[replicate_key] == i]
                    fraction = len(high_rep[high_rep[:,gene].X.todense() > 0].obs.index)/len(high_rep.obs.index)
                    if fraction > 0.33:
                        high_33 += 1
                high_fraction = high_33/high_total
                high_fraction_list.append(high_fraction)

                lower_group = FC[FC['gene']==gene]['lower_group']
                lower_group = lower_group.iloc[0]
                low = adata_temp[adata_temp.obs[metacell_axis]== lower_group]
                low_total = len(unique(low.obs[replicate_key]))
                low_reps = unique(low.obs[replicate_key])
                low_33 = 0
                for i in low_reps:
                    low_rep = low[low.obs[replicate_key]==i]
                    fraction = len(low_rep[low_rep[:,gene].X.todense() > 0].obs.index)/len(low_rep.obs.index)
                    if fraction > min_fraction:
                        low_33 += 1
                low_fraction = low_33/low_total
                low_fraction_list.append(low_fraction)
    
            #zip up the outputs and create a dataframe
            rep_condensed_output = list(zip(replist,comps,ngenes,nDEGs))
            rep_long_output = list(zip(replist_long,complist_long,higher_all_rep,lower_all_rep,DEGs_all_rep,lfc_all_rep, 
                                       high_fraction_list,low_fraction_list))
    
            #convert to dataframe format
            condensed_DEG_df_rep = pd.DataFrame(rep_condensed_output,
                                        columns=['replicate','comparison','ngenes','nDEGs'])
            long_DEG_df_rep = pd.DataFrame(rep_long_output,
                                        columns=['replicate','comparison','higher_group','lower_group','DEG_name', 'logFC', 'higher_group_reps', 'lower_group_reps'])       
    

        condensed_DEG_df = condensed_DEG_df_rep
        long_DEG_df = long_DEG_df_rep

            
        return condensed_DEG_df,long_DEG_df
