import os
import json
import pandas as pd
import pathlib as pl


root_dir = pl.Path(os.getcwd())
inputs_dir = root_dir/'Inputs/DC'
outputs_dir = root_dir/'Outputs/DC'
pluvial_events_path = inputs_dir/'FFRD_DC_All_Pluvial_Events.xlsx'
pluvial_rep_events_path = outputs_dir/'FFRD_DC_Representative_Pluvial_Events.xlsx'


weight_files = []
for f in inputs_dir.glob('*.json'):
    if 'Weights_TRI' in f.stem: 
        weight_files.append(f)       
        print(f)

writer = pd.ExcelWriter(pluvial_rep_events_path)

for file in weight_files:
    sheet_name = file.stem.split('_Weights')[0]
    with open(file) as f:
        dic = json.load(f)
    mainbcn = list(dic['BCName'].keys())[0]
    weights_dic = dic['BCName'][mainbcn]
    df = pd.read_excel(pluvial_events_path, sheet_name = sheet_name, index_col = 0)
    rep_df = df.loc[list(weights_dic.keys())].copy()
    rep_df['Weight'] = [weights_dic[i] for i in rep_df.index]
    rep_df.sort_index(inplace = True)
    rep_df.to_excel(writer, sheet_name)

writer.save()
writer.close()