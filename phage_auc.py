from typing import List, Dict, Tuple
import numpy as np
import pandas as pd

import os
import sys
import argparse


class dataSchema():
    def __init__(self, excel_datafile: str, schema_sheetname="Schema", ignore_bacteria=['KP85']) -> None:
        self.excel_datafile = excel_datafile
        self.schema_sheetname = schema_sheetname

        self.platemap = self.get_platemap()
        self.replicate_list = list(self.platemap['Replicate'].unique())

        self.ignored_bacteria = ignore_bacteria
        ignored_set = set(["-"] + self.ignored_bacteria) 
        self.bacteria_list = list(set(self.platemap['sp1_ID'].unique()) - ignored_set)
    
    def get_platemap(self) -> pd.DataFrame:
        platemap = pd.read_excel(self.excel_datafile, sheet_name=self.schema_sheetname)
        platemap['cut-off_time'] = platemap['cut-off_time'].apply(lambda x: x.strip())
        return platemap

    def get_platemap_bacteria(self, bacteria_ID: str) -> pd.DataFrame:
        df = self.platemap.loc[self.platemap['sp1_ID']==bacteria_ID, :]
        return df
    
    def get_platemap_replicate(self, replicate_name: str) -> pd.DataFrame:
        df = self.platemap.loc[self.platemap['Replicate']==replicate_name,:]
        return df

class dataPerReplicate():
    def __init__(self, excel_datafile: str, replicate_name: str, 
                 t_start=0.0, t_end=24.0, t_interval=5.0/60.0) -> None:
        self.excel_datafile = excel_datafile 
        self.replicate_name = replicate_name

        self.t_start = t_start
        self.t_end = t_end
        self.t_interval = t_interval
        self.n_timepoints = int((t_end - t_start)/t_interval + 1)

        self.data_table, self.t_labels = self.get_data_table()
        self.t_floats = np.linspace(self.t_start, self.t_end, len(self.t_labels))
        self.t_dict = dict(zip(self.t_labels, self.t_floats))
    
    def get_data_table(self) -> (pd.DataFrame, List):
        df = pd.read_excel(self.excel_datafile,
                           sheet_name=self.replicate_name, skiprows=[0],
                           usecols=[0, *range(2, self.n_timepoints+2)])
        t_labels = [x.strip() for x in df.columns[1:]]
        df.columns = ['Well'] + t_labels
        return df, t_labels
    
    def AUC_per_bacteria(self, bacteria_ID: str, platemap: pd.DataFrame) -> pd.DataFrame:
        replicate_platemap = platemap.loc[platemap['Replicate']==self.replicate_name, :]
        schema = replicate_platemap.loc[replicate_platemap['sp1_ID']==bacteria_ID, :]

        blank_wells = schema.loc[schema['culture_type'].str.contains("Blank"), 'Well']
        sample_wells = schema.loc[~schema['Well'].isin(blank_wells), 'Well']

        cutoff_list = list(schema.loc[schema['culture_type'].str.contains("Positive control"), 'cut-off_time'].unique())
        if len(cutoff_list) > 1:
            raise Exception("cut-off time values must be similar for same positive controls within a replicate")
        else:
            cutoff = cutoff_list[0]

        auc_of_bact = schema.loc[schema['Well'].isin(sample_wells), 
                                 ['Replicate', 'Well', 'culture_type', 
                                  'sp1_ID', 'sp2_ID', 'cut-off_time']].reset_index(drop=True)
        
        auc_of_bact['AUC'] = self.compute_AUC(sample_wells=sample_wells, blank_wells=blank_wells,
                                              cutoff=cutoff, blank_corrected=False)

        auc_of_bact['AUC_blank_corrected'] = self.compute_AUC(sample_wells=sample_wells, blank_wells=blank_wells,
                                                        cutoff=cutoff, blank_corrected=True)        

        return auc_of_bact

    def compute_AUC(self, sample_wells: List, blank_wells: List, cutoff: str, 
                    blank_corrected=True) -> pd.Series:
        
        selected_time = self.t_labels[0:(self.t_labels.index(cutoff)+1)]
        ts = self.t_floats[0:(self.t_labels.index(cutoff)+1)]

        sample_data = self.data_table.loc[self.data_table['Well'].isin(sample_wells), selected_time]
        blank_data = self.data_table.loc[self.data_table['Well'].isin(blank_wells), selected_time]

        if blank_corrected:
            input_data = sample_data - blank_data.mean()
        else:
            input_data = sample_data
        
        auc = input_data.reset_index(drop=True).apply(lambda y: np.trapz(y, x=ts), axis=1)

        return auc

def main():
    args = parse_arguments()
    
    schema = dataSchema(excel_datafile=args.data_file, 
                        schema_sheetname=args.schema_sheetname,
                        ignore_bacteria=args.ignore_bacteria)
    
    rep_list = schema.replicate_list
    data_obj_list = []
    for rep in rep_list:
        data_obj = dataPerReplicate(excel_datafile=args.data_file,
                                    replicate_name=rep,
                                    t_start=args.t_start,
                                    t_end=args.t_end,
                                    t_interval=args.t_interval)
        data_obj_list.append(data_obj)
    
    bact_list = schema.bacteria_list

    output_table = pd.DataFrame()
    for bact in bact_list:
        for rep_data in data_obj_list:
            auc_table = rep_data.AUC_per_bacteria(bacteria_ID=bact, platemap=schema.platemap)
            output_table = pd.concat([output_table, auc_table], ignore_index=True)

    if args.average:
        write_table = output_table.groupby(['sp1_ID', 'sp2_ID', 'culture_type'])['AUC', 'AUC_blank_corrected'].mean().reset_index()
    else:
        write_table = output_table
    
    write_table.to_csv(args.out_file, index=False)

    print("Calculations ran successful! Output has been saved to output file: {}".format(args.out_file))

def parse_arguments():
    parser = argparse.ArgumentParser(description="Calculate AUC of growth curves with given cut-off timepoints.")
    
    parser.add_argument("-f", "--data_file", type=str, required=True,
                        help="Excel input file [required]")
    
    parser.add_argument("-o", "--out_file", type=str, 
                        default="Phage_AUC.csv", required=False,
                        help="CSV output file [required]")
    
    parser.add_argument("-a", "--average", action='store_true', required=False,
                        help="Enable averaging AUC values across replicates if specified.")

    parser.add_argument("-s", "--schema_sheetname", 
                        default="Schema", required=False,
                        help="Name of the Excel sheet for the schema of culture plates. [default: \'Schema\']")
    
    parser.add_argument("-ts", "--t_start", type=float, 
                        default=0.0, required=False,
                        help="Start time of experiments, in hour. [default: 0.0]")
    
    parser.add_argument("-te", "--t_end", type=float, 
                        default=24.0, required=False,
                        help="End time of experiments, in hour. [default: 24.0]")
    
    parser.add_argument("-ti", "--t_interval", type=float,
                        default=1/12, required=False,
                        help="Time interval between every 2 consecutive OD measurements, in hour. [default: 1/12 (5 min)]")
    
    parser.add_argument("-i", "--ignore_bacteria", action='append',
                        default=['KP85'], required=False,
                        help="ID of bacteria (sp1_ID) to ignore in calculations. Use multiple times when there're >1 bacteria to ignore.\n[default: 'KP85']")

    # If no arguments were provided
    if len(sys.argv) == 1:
        # parser.print_help(file=sys.stderr)
        try:
            with open("dump.txt") as lovely:
                for line in lovely:
                    print(line.strip('\n'))
        except:
            print("Hey ! You missed your arguments! AND MY F*CKING DUMP FILE!\n   ... It's dump.txt. Please find it!\n   ... Immediately!")
        sys.exit(1)

    args = parser.parse_args()

    return args

if __name__ == '__main__':
    main()