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
import matplotlib.pyplot as plt
import libraries.Tool_Box as ToolBox

def plot_grouped_histogram(output_folder, plotting_df):
    grouped = plotting_df.groupby(["Sample", "Target"], sort=False, as_index=False).agg(
        average_data=('Data', 'mean'), sem_data=('Data', 'sem')
    )
    targets = grouped["Target"].unique()
    if len(targets) == 0:
        raise ValueError("No data available to plot")

    for target in targets:
        fig, ax = plt.subplots(figsize=(10, 5))

        target_grouped = grouped[grouped["Target"] == target]
        samples = target_grouped["Sample"].tolist()

        # seaborn.barplot(data=target_grouped, x="Sample", y="average_data", order=samples, ax=ax)

        for position, row in enumerate(target_grouped.itertuples(index=False)):
            if pd.notna(row.sem_data):
                ax.errorbar(
                    position,
                    row.average_data,
                    yerr=row.sem_data,
                    fmt="none",
                    elinewidth=0.5,
                    color="black",
                    capsize=5,
                    capthick=0.5
                )

        seaborn.barplot(
            data=target_grouped,
            x="Sample",
            y="average_data",
            hue="Sample",
            palette="rocket",
            order=samples, ax=ax
        )

        seaborn.stripplot(
            data=plotting_df[plotting_df["Target"] == target],
            x="Sample",
            y="Data",
            order=samples,
            color="black",
            size=2.5,
            marker="o",
            jitter=False,
            ax=ax
        )

        ax.set_title(str(target))
        ax.set_xticks(range(len(samples)))
        ax.set_xlabel("")
        ax.set_xticklabels(samples, rotation=30, ha="right")
        ax.set_ylabel("Mean +/- SEM normalized data")

        fig.tight_layout()

        # Sanitize target name in case it contains characters invalid for Windows filenames (like / or :)
        safe_target = str(target).replace("/", "_").replace(":", "_")
        output_filepath = "{}{}_plot.png".format(output_folder, safe_target)
        # filepath = os.path.join("D:\\", f"{safe_target}_Test_plot.png")

        fig.savefig(output_filepath, dpi=600, bbox_inches="tight")
        plt.show()
        plt.close(fig)

def original_plot_grouped_histogram(plotting_df):
    grouped = plotting_df.groupby(["Sample", "Target"], as_index=False).agg(
        average_data=('Data', 'mean'), sem_data=('Data', 'sem')
    )

    # Not sure this is really needed.  Since grouped is already a DataFrame, it should have a unique "Target" column.
    targets = grouped["Target"].unique()
    if len(targets) == 0:
        raise ValueError("No data available to plot")

    fig, axes = plt.subplots(len(targets), 1, figsize=(10, 5 * len(targets)), squeeze=False)

    for ax, target in zip(axes.flat, targets):
        target_grouped = grouped[grouped["Target"] == target]
        samples = target_grouped["Sample"].tolist()
        seaborn.barplot(data=target_grouped, x="Sample", y="average_data", order=samples, ax=ax)
        for position, row in enumerate(target_grouped.itertuples(index=False)):
            if pd.notna(row.sem_data):
                ax.errorbar(position, row.average_data, yerr=row.sem_data, fmt="none", color="black", capsize=4)
        seaborn.stripplot(data=plotting_df[plotting_df["Target"] == target], x="Sample", y="Data",
                          order=samples, color="black", marker="o", jitter=False, ax=ax)
        ax.set_title(str(target))
        ax.set_xticks(samples)
        ax.set_xticklabels(samples, rotation=30, ha="right")
        ax.set_ylabel("Mean +/- SEM normalized data")
        fig.tight_layout()
        plt.savefig("D:{}{}_Test_plot.png".format(os.sep, target), dpi=600, bbox_inches="tight")
        plt.show()
    plt.close()
    return fig, axes.flat

class DataProcessing:
    def __init__(self, tsv_file, ouput_folder, channel_label_list):
        self.tsv_file = tsv_file
        self.ouput_folder = ouput_folder
        self.column_labels = channel_label_list
        self.fam_label = channel_label_list[0]
        self.hex_label = channel_label_list[1]
        self.cy5_label = channel_label_list[2]
        self.cy55_label = channel_label_list[3]
        self.rox_label = channel_label_list[4]
        self.atto590_label = channel_label_list[5]
        self.main()

    def main(self):

        ToolBox.debug_messenger("We here?")
        input_data_file = csv.reader(open(self.tsv_file), delimiter='\t')

        # ToDo: This section need to be tied to a GUI.
        sample = ""
        """
        fam_label = "NHEJ"
        hex_label = "Genomic"
        rox_label = "Germline"
        cy5_label = "LBR2.1"
        cy55_label = "Not Used"
        atto590_label = "TMEJ-AAG"
        """
        replicates = 2
        output_header = "Sample\t"
        #column_label_list = [fam_label, hex_label, cy5_label, cy55_label, rox_label, atto590_label]

        # Make header string for output file
        for i in range(len(self.column_labels)):
            for j in range(replicates):
                output_header += ("{}\t".format(self.column_labels[i]))
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

        # Parse data file
        for row in input_data_file:
            # skip header row
            if row_count == 0:
                row_count += 1
                continue

            well = row[0]

            if row[1] != "No Template" and row[2] == "PathSig2":
                # Strip off the last two characters of the sample name.  This removes the sample duplicate number.
                sample = (row[1][:-2])
            elif row[1] == "No Template" and row[2] == "PathSig2":
                sample = "No Template"

            # Get copies per uL from data.
            copies = float(row[6])
            if row[5] == self.fam_label:
                fam_copies = copies
            elif row[5] == self.hex_label:
                hex_copies = copies
            elif row[5] == self.cy5_label:
                cy5_copies = copies
            elif row[5] == self.cy55_label:
                cy55_copies = copies
            elif row[5] == self.rox_label:
                rox_copies = copies
            elif row[5] == self.atto590_label:
                atto590_copies = copies

            # Assign normal control value based on control name
            if controlname == self.fam_label:
                normalcontrol = fam_copies
            elif controlname == self.hex_label:
                normalcontrol = hex_copies
            elif controlname == self.cy5_label:
                normalcontrol = cy5_copies
            elif controlname == self.cy55_label:
                normalcontrol = cy55_copies
            elif controlname == self.rox_label:
                normalcontrol = rox_copies
            elif controlname == self.atto590_label:
                normalcontrol = atto590_copies
            else:
                # ToDo: This should display in the GUI.
                print("Warning: No normal control defined.")
                exit("Warning: No normal control defined.")

            # Normalize data to the normal control.  It is possible that one or more of the values are not defined.
            # If that is the case, set the normalized value to an empty string.
            if normalcontrol and isinstance(fam_copies, (int, float)):
                normalfam = fam_copies / normalcontrol
            else:
                normalfam = ""

            if normalcontrol and isinstance(hex_copies, (int, float)):
                normalhex = hex_copies / normalcontrol
            else:
                normalhex = ""

            if normalcontrol and isinstance(cy5_copies, (int, float)):
                normalcy5 = cy5_copies / normalcontrol
            else:
                normalcy5 = ""

            if normalcontrol and isinstance(cy55_copies, (int, float)):
                normalcy55 = cy55_copies / normalcontrol
            else:
                normalcy55 = ""

            if normalcontrol and isinstance(rox_copies, (int, float)):
                normalrox = rox_copies / normalcontrol
            else:
                normalrox = ""

            if normalcontrol and isinstance(atto590_copies, (int, float)):
                normalatto590 = atto590_copies / normalcontrol
            else:
                normalatto590 = ""

            # Build a dictionary entry for the well of raw data and normalized values.
            normalized_values = [normalfam, normalhex, normalcy5, normalcy55, normalrox, normalatto590]
            raw_data_dict[well] = \
                [[sample, fam_copies, hex_copies, cy5_copies, cy55_copies, rox_copies, atto590_copies],
                 normalized_values]

        # Make sure the dictionary is sorted by wells.  The ddPCR data is grouped by wells.
        # raw_data_dict = dict(sorted(raw_data_dict.items()))
        normalized_sample_data_dict = defaultdict(list)
        output_data = ""
        data_counter = 0
        ploting_data_dict = {}

        sorted_raw_data_dict = dict(sorted(raw_data_dict.items(), key=lambda item: (int(item[0][1:]), item[0][0])))
        for well in sorted_raw_data_dict:
            sample_name = raw_data_dict[well][0][0]
            normalized_values = raw_data_dict[well][1]
            output_data += ("{}\t{}\t{}\t{}\t{}\t{}\t{}\n"
                            .format(sample_name, normalized_values[0], normalized_values[1], normalized_values[2],
                                    normalized_values[3], normalized_values[4], normalized_values[5]))

            # This will be converted to pandas later
            for label, normal_value in zip(self.column_labels, normalized_values):
                if isinstance(normal_value, (int, float)):
                    ploting_data_dict[data_counter] = sample_name, label, normal_value
                    data_counter += 1

            # ToDo: this dictionary is for tsv file output.  Refactor.  Should be able to get data directly into string.
            """
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
            """
        # df = pd.DataFrame.from_dict(normalized_sample_data_dict, orient='index', columns=column_label_list)

        # sorted_data = dict(sorted(ploting_data_dict.items(), key=lambda item: item[1][0]))
        plotting_df_column_labels = ["Sample", "Target", "Data"]

        plotting_df = pd.DataFrame.from_dict(ploting_data_dict, orient='index', columns=plotting_df_column_labels)

        # Remove rows that contain "Genomic", "No RNP", or "No Template" in the "Target" column.  These are not useful
        # for plotting
        plotting_df = plotting_df[~plotting_df["Target"].str.contains("Genomic", case=False, na=False)]
        plotting_df = plotting_df[~plotting_df["Sample"].str.contains(" No RNP", case=False, na=False)]
        plotting_df = plotting_df[~plotting_df["Sample"].str.contains("No Template", case=False, na=False)]

        plot_grouped_histogram(self.ouput_folder, plotting_df)

        outfile = open("{}Normalized_Data.tsv".format(self.ouput_folder), "w")
        outfile.write(output_header + output_data)
        outfile.close()

        exit("Last Dragons")
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
                data_list[i] = df.loc[sample_name, column_label]

            fam_data = normalized_sample_data_dict[sample_name][0]
            hex_data = normalized_sample_data_dict[sample_name][1]
            cy5_data = normalized_sample_data_dict[sample_name][2]
            cy55_data = normalized_sample_data_dict[sample_name][3]
            rox_data = normalized_sample_data_dict[sample_name][4]
            atto590_data = normalized_sample_data_dict[sample_name][5]

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
    DataProcessing(tsv_file = "D:{}ddPCR_Test_Data.tsv".format(os.sep))
