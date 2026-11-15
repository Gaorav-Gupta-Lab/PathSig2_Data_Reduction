"""
Script to normalize PathSig2 data.
"""

import csv
import os
from collections import defaultdict
import scipy.stats as stats
import numpy
from contextlib import suppress as Suppress
import pandas as pd
import seaborn
import Tool_Box as ToolBox

def main(file_path):
    # ToDo: This section need to be tied to a GUI.
    input_data_file = csv.reader(open(file_path), delimiter='\t')
    sample = ""
    fam_label = "NHEJ"
    hex_label = "Genomic"
    rox_label = "Germline"
    cy5_label = "LBR2.1"
    cy55_label = "Not Used"
    atto590_label = "TMEJ-AAG"
    replicates = 2
    output_header = "Sample\t"
    column_label_list = [fam_label, hex_label, cy5_label, cy55_label, rox_label, atto590_label]

    # Make header string for output file
    for i in range(len(column_label_list)):
        for j in range(replicates):
            output_header += ("{}\t".format(column_label_list[i]))
    output_header += "\n"

    fam_copies = ""
    hex_copies = ""
    cy5_copies = ""
    rox_copies = ""
    cy55_copies = ""
    atto590_copies = ""

    controlname = "Genomic"
    normalcontrol = ""

    raw_data_dict = {}
    row_count = 0

    for row in input_data_file:
        if row_count == 0:
            row_count += 1
            continue

        well = row[0]

        if row[1] != "No Template" and row[2] == "PathSig2":
            sample = (row[1][:-2])
        elif row[1] == "No Template" and row[2] == "PathSig2":
            sample = "No Template"

        copies = float(row[6])
        if row[5] == fam_label:
            fam_copies = copies
        elif row[5] == hex_label:
            hex_copies = copies
        elif row[5] == cy5_label:
            cy5_copies = copies
        elif row[5] == cy55_label:
            cy55_copies = copies
        elif row[5] == rox_label:
            rox_copies = copies
        elif row[5] == atto590_label:
            atto590_copies = copies

        if controlname == fam_label:
            normalcontrol = fam_copies
        elif controlname == hex_label:
            normalcontrol = hex_copies
        elif controlname == cy5_label:
            normalcontrol = cy5_copies
        elif controlname == cy55_label:
            normalcontrol = cy55_copies
        elif controlname == rox_label:
            normalcontrol = rox_copies
        elif controlname == atto590_label:
            normalcontrol = atto590_copies
        else:
            # ToDo: This should display in the GUI.
            print("Warning: No normal control defined.")
            exit("Warning: No normal control defined.")

        if not isinstance(normalcontrol, (int, float)):
            normalfam = 0
            normalhex = 0
            normalcy5 = 0
            normalcy55 = 0
            normalrox = 0
            normalatto590 = 0

        else:
            # It is possible that one or more of the values are not defined.
            if isinstance(fam_copies, (int, float)):
                normalfam = fam_copies / normalcontrol
            else:
                normalfam = ""

            if isinstance(hex_copies, (int, float)):
                normalhex = hex_copies / normalcontrol
            else:
                normalhex = ""

            if isinstance(cy5_copies, (int, float)):
                normalcy5 = cy5_copies / normalcontrol
            else:
                normalcy5 = ""

            if isinstance(cy55_copies, (int, float)):
                normalcy55 = cy55_copies / normalcontrol
            else:
                normalcy55 = ""

            if isinstance(rox_copies, (int, float)):
                normalrox = rox_copies / normalcontrol
            else:
                normalrox = ""

            if isinstance(atto590_copies, (int, float)):
                normalatto590 = atto590_copies / normalcontrol
            else:
                normalatto590 = ""

            normalized_values = [normalfam, normalhex, normalcy5, normalcy55, normalrox, normalatto590]
            raw_data_dict[well] = [[sample, fam_copies, hex_copies, cy5_copies, cy55_copies, rox_copies,
                                    atto590_copies], normalized_values]

    # Make sure the dictionary is sorted by wells.  The ddPCR data is grouped by wells.
    raw_data_dict = dict(sorted(raw_data_dict.items()))
    normalized_sample_data_dict = defaultdict(list)
    output_data = ""

    for well in raw_data_dict:
        sample_name = raw_data_dict[well][0][0]
        normalized_values = raw_data_dict[well][1]
        output_data += ("{}\t{}\t{}\t{}\t{}\t{}\t{}\n"
                        .format(sample_name, normalized_values[0], normalized_values[1], normalized_values[2],
                                normalized_values[3], normalized_values[4], normalized_values[5]))

        if sample_name not in normalized_sample_data_dict:
            normalized_sample_data_dict[sample_name].append([normalized_values[0]])
            normalized_sample_data_dict[sample_name].append([normalized_values[1]])
            normalized_sample_data_dict[sample_name].append([normalized_values[2]])
            normalized_sample_data_dict[sample_name].append([normalized_values[3]])
            normalized_sample_data_dict[sample_name].append([normalized_values[4]])
            normalized_sample_data_dict[sample_name].append([normalized_values[5]])
        else:
            normalized_sample_data_dict[sample_name][0].append(normalized_values[0])
            normalized_sample_data_dict[sample_name][1].append(normalized_values[1])
            normalized_sample_data_dict[sample_name][2].append(normalized_values[2])
            normalized_sample_data_dict[sample_name][3].append(normalized_values[3])
            normalized_sample_data_dict[sample_name][4].append(normalized_values[4])
            normalized_sample_data_dict[sample_name][5].append(normalized_values[5])

    df = pd.DataFrame.from_dict(normalized_sample_data_dict, orient='index', columns=column_label_list)

    outfile = open("D:{}normalized_data.tsv".format(os.sep), "w")
    outfile.write(output_header+output_data)
    outfile.close()

    fam_data = ""
    hex_data = ""
    cy5_data = ""
    cy55_data = ""
    rox_data = ""
    atto590_data = ""
    data_list = [fam_data, hex_data, cy5_data, cy55_data, rox_data, atto590_data]
    for sample_name in df.index:
        for i, column_label in enumerate(column_label_list):
            # ToolBox.debug_messenger(df.loc[sample_name:,column_label])
            data_list[i] = df.loc[sample_name,column_label]

        fam_data = normalized_sample_data_dict[sample_name][0]
        hex_data = normalized_sample_data_dict[sample_name][1]
        cy5_data = normalized_sample_data_dict[sample_name][2]
        cy55_data = normalized_sample_data_dict[sample_name][3]
        rox_data = normalized_sample_data_dict[sample_name][4]
        atto590_data = normalized_sample_data_dict[sample_name][5]
        # ToolBox.debug_messenger(data_list)

        # In case any data is None, this allows the program to continue without crashing.
        with Suppress(TypeError):
            print("FAM ", sample_name, numpy.mean(fam_data), stats.sem(fam_data))
            print("HEX", sample_name, numpy.mean(hex_data), stats.sem(hex_data))
            print("Cy5", sample_name, numpy.mean(cy5_data), stats.sem(cy5_data))
            print("Cy55", sample_name, numpy.mean(cy55_data), stats.sem(cy55_data))
            print("ROX", sample_name, numpy.mean(rox_data), stats.sem(rox_data))
            print("Atto590", sample_name, numpy.mean(atto590_data), stats.sem(atto590_data))
            continue
            
if __name__ == '__main__':
    main(file_path = "D:{}ddPCR_Test_Data.tsv".format(os.sep))
