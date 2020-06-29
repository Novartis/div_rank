#!/usr/bin/env python
# coding: utf-8


import pandas as pd  #we want pandas 0.25 here. The code has not been tested with pandas 1.0
import numpy as np
import os,sys
import json

def read_rank_data(config):
    path = config['path']
    rank_file = config['rank_file']
    sample_col = str(config['sample_name'])
    mol_rank_df = pd.read_csv(os.path.join(path,rank_file),low_memory=False)
    mol_rank_df.set_index(sample_col,inplace=True)
    if config['score_column'] != 'score':
        mol_rank_df['score'] = mol_rank_df[config['score_column']]
    if 'random_id_column' in config:
        mol_rank_df['random_id'] = mol_rank_df[config['random_id_column']]
    else:
        mol_ids = np.arange(len(mol_rank_df))
        np.random.shuffle(mol_ids)
        mol_rank_df['random_id'] = mol_ids    
    return(mol_rank_df)
    

def get_individual_class_df(path,file,sep,sample_name,class_id,min_size,score_data):
    df = pd.read_table(os.path.join(path,'cluster_data',file),sep=sep).rename(columns={sample_name:'sample_name',class_id:'class_id'})[['sample_name','class_id']]
    df = df.join(score_data,on='sample_name',how='inner')
    df['class_size'] = df.groupby('class_id')['class_id'].transform('size')
    #we now dertemine from which picking roun on class may be unlocked. The memebership quorum for unlocking at class at a given 
    #round is defined per class type in the configuration file     
    df['unlock_round'] = -1
    for rnd,class_size in reversed(list(enumerate(min_size))):
        df.loc[df['class_size'] >= class_size,'unlock_round'] = rnd
    #remove any records with out unlock round. This will happen if the membership quorum is always above 1 for all rounds
    df = df[df['unlock_round'] >= 0]
    print('class records: {0}, classified_mols: {1}, number of classes:{2} ({3})'.format(len(df),df['sample_name'].nunique(),df['class_id'].nunique(),file))
    return(df)
            

# assemble the overal class data from the individual files
def get_global_class_df(config, mol_rank_df):
    score_data = mol_rank_df[['score','random_id']]
    class_types_config = config['class_data']
    class_data = pd.concat({class_type : get_individual_class_df(path=config['path'],score_data=score_data,**cf) for class_type,cf in class_types_config.items() },names=['class_type'])
    class_data = class_data.reset_index()[['class_type','class_id','sample_name','score','random_id','unlock_round']]
    class_data['global_class_id'] = class_data.groupby(['class_type','class_id']).ngroup()

    missing= ~mol_rank_df.index.isin(class_data['sample_name'])
    print('Number of samples without any class record: {0}. These will picked last. Theres hould not be a significant numbr of samples here'.format(missing.sum()))

    return(class_data)

#execute the ranking. This may modify mol_rank_df and class_data in place, depending on the pandas 
#performing a copyoperation or not
    
def execute_div_rank(mol_rank_df,class_data,target_picking_quorum):
    #create columns to record the div rank results
    mol_rank_df['pick_seq'] = 0
    mol_rank_df['pick_round'] = 0
    mol_rank_df = mol_rank_df.reindex(columns=(list(mol_rank_df.columns)+list(class_data['class_type'].unique())),fill_value = 0)

    #add a picked flag column to class data
    class_data['picked'] = False

    #initialze counters
    current_pick_round = 1
    current_pick_seq = 0
    current_picking_quorum = 0.0

    #main picking loop
    while (mol_rank_df['pick_round']==0).any() and (current_picking_quorum < target_picking_quorum):
        print('Picking round:{0}'.format(current_pick_round))
        #get class data for molecules not yet picked and the classes already unlocked for this round
        pick_candidates = class_data[~class_data['picked'] & (class_data['unlock_round'] < current_pick_round)].copy()
        #group them by global class id
        cand_gp = pick_candidates.groupby('global_class_id')
        #rank them within group by their score
        pick_candidates['pick_rank'] = cand_gp['score'].rank(method='dense').astype(int)
        #retain only the top ranking class members
        #we make a copy as we will manipulate this ubframe
        pick_candidates = pick_candidates[pick_candidates['pick_rank'] == 1].copy()
        #if this is not the first picking round, we determine how many compounds with the the best picking score
        #or better there have already been picked for this class. If there more then current_pick_round, 
        #we skip the class
        if current_pick_round > 1:
            #find best pickable scores
            current_pick_score = cand_gp['score'].min().rename('current_pick_score')
            #join these to the already picked classes by global_class_id
            class_data_picked = class_data[class_data['picked']].join(current_pick_score,on=('global_class_id'),how='inner')
            #for each class count the already picked compounds with a score btter or equal to the current one
            #include the unloack round in the grouping (which is constant over a class id) to have it availabale over the grouping
            pick_count_data = class_data_picked[class_data_picked['score'] <= class_data_picked['current_pick_score']].groupby(['global_class_id','unlock_round']).size()
            pick_count_data.name = 'pick_count'
            pick_count_data = pick_count_data.reset_index(level='unlock_round',drop=False)
            #generate list of global-class_ids to be skipped as they are already above par
            #par measn the the picked number of members is larger than the current_pick_round minis the unlock_round for te class
            #if a cluass is unlocked at round 1, and at beginning of round 3 we have picked 2 members, we are aready avove par as 2 >= 3-1
            pick_exclude = pick_count_data[pick_count_data['pick_count']>=(current_pick_round-pick_count_data['unlock_round'])].index
            #remove these classes from the pick canddiates
            pick_candidates = pick_candidates[pick_candidates['global_class_id'].isin(pick_exclude)].copy()
        #in the remaining candidates count now the number of class types a molecules covers
        pick_candidates['class_type_counts'] = pick_candidates.groupby('random_id')['class_type'].transform('nunique')
        #sort by score, andn then by class_type_count, and then by the random id
        #the last sort criterion ensures that if several classes have a molecules in commnon which rank first according the first two criteria
        #we indeed pick the same, but otherwise randomly chose molecules
        pick_candidates.sort_values(by=['score','class_type_counts','random_id'],ascending=[True,False,True],inplace=True)
        #now we keep only te first molecules pe class
        pick_candidates.drop_duplicates(subset='global_class_id',keep='first',inplace=True)
        #the picked molecules are the unique list of cencepst for the top molecules per class
        #as the pick candidates are sorted by score, so are the molecules
        pick_mols = pick_candidates['sample_name'].drop_duplicates()
        print('Picked {0} molecules covering {1} classes'.format(len(pick_mols),len(pick_candidates)))
        #record the pick in the mol_rank_df
        mol_rank_df.loc[pick_mols,'pick_round'] = current_pick_round
        mol_rank_df.loc[pick_mols,'pick_seq'] = np.arange(current_pick_seq,current_pick_seq+len(pick_mols))
        #record also the class_type which lead to the pick
        pick_type_count = pick_candidates.groupby(['sample_name','class_type']).size().unstack('class_type',fill_value=0)
        mol_rank_df.loc[pick_type_count.index,pick_type_count.columns] = pick_type_count
        #update the 'picked' column in class_data 
        class_data.loc[class_data['sample_name'].isin(pick_mols),'picked'] = True
        #update the round counter
        current_pick_round += 1
        current_pick_seq += len(pick_mols)
        current_picking_quorum = float(current_pick_seq)/len(mol_rank_df)
        print('Current picking quorum: {0}'.format(current_picking_quorum))

    #wrap up: assign the a pick rank to the not yet picked compounds. These will be 
    mol_rank_df.loc[mol_rank_df['pick_round']==0,'pick_round'] = current_pick_round
    pick_seq = mol_rank_df.loc[mol_rank_df['pick_round']==current_pick_round,'score'].rank(method='first').astype(int)+current_pick_seq
    mol_rank_df.loc[pick_seq.index,'pick_seq']=pick_seq
    return(mol_rank_df)



def main():
    try:
        config_file = sys.argv[1]
    except: 
        raise ValueError('Usage: div_rank.py <config_file>')

    with open(config_file,'r') as jsonf:
        config = json.load(jsonf)

    print('read in config data')
    mol_rank_df = read_rank_data(config)
    print('finished reading in rank data')
    class_data = get_global_class_df(config, mol_rank_df)
    print('finished reading in class data')
    mol_rank_df = execute_div_rank(mol_rank_df,class_data,target_picking_quorum=config['target_picking_quorum'])
    print('finished picking')
    mol_rank_df.to_csv(os.path.join(config['path'],config['outfile']))
    print('written output')
    print('DONE')


if __name__ == "__main__":
    main()