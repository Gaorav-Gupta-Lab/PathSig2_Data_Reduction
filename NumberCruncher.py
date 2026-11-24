"""
Script to normalize PathSig2 data.
"""

import csv
import os
import pandas as pd
import seaborn
import matplotlib.pyplot as plt
import libraries.Tool_Box as ToolBox

def plot_grouped_histogram(output_folder, plotting_df):
    """
    Main plotting function for grouped histograms.
    """
    grouped = plotting_df.groupby(["Sample", "Target"], sort=False, as_index=False).agg(
        average_data=('Data', 'mean'), sem_data=('Data', 'sem')
    )
    targets = grouped["Target"].unique()
    if len(targets) == 0:
        raise ValueError("No data available to plot")

    # Create one plot for each target
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

        fig.savefig(output_filepath, dpi=600, bbox_inches="tight")
        plt.show()
        plt.close(fig)

def original_plot_grouped_histogram(plotting_df):
    """
    Not used
    """
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

        sample = ""
        # ToDo: replicates should be from GUI or count number of samples and generate the header below.
        replicates = 2
        output_header = "Sample\t"

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
        normalcontrol_data = ""

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
                normalcontrol_data = fam_copies
            elif controlname == self.hex_label:
                normalcontrol_data = hex_copies
            elif controlname == self.cy5_label:
                normalcontrol_data = cy5_copies
            elif controlname == self.cy55_label:
                normalcontrol_data = cy55_copies
            elif controlname == self.rox_label:
                normalcontrol_data = rox_copies
            elif controlname == self.atto590_label:
                normalcontrol_data = atto590_copies
            else:
                # ToDo: This should display in the GUI.
                print("Warning: No normal control defined.")
                exit("Warning: No normal control defined.")

            # Normalize data to the normal control.  It is possible that one or more of the values are not defined.
            # If that is the case, set the normalized value to an empty string.
            if normalcontrol_data and isinstance(fam_copies, (int, float)):
                normalfam = fam_copies / float(normalcontrol_data)
            else:
                normalfam = ""

            if normalcontrol_data and isinstance(hex_copies, (int, float)):
                normalhex = hex_copies / float(normalcontrol_data)
            else:
                normalhex = ""

            if normalcontrol_data and isinstance(cy5_copies, (int, float)):
                normalcy5 = cy5_copies / float(normalcontrol_data)
            else:
                normalcy5 = ""

            if normalcontrol_data and isinstance(cy55_copies, (int, float)):
                normalcy55 = cy55_copies / float(normalcontrol_data)
            else:
                normalcy55 = ""

            if normalcontrol_data and isinstance(rox_copies, (int, float)):
                normalrox = rox_copies / float(normalcontrol_data)
            else:
                normalrox = ""

            if normalcontrol_data and isinstance(atto590_copies, (int, float)):
                normalatto590 = atto590_copies /  float(normalcontrol_data)
            else:
                normalatto590 = ""

            # Build a dictionary entry for the well of raw data and normalized values.
            normalized_values = [normalfam, normalhex, normalcy5, normalcy55, normalrox, normalatto590]
            raw_data_dict[well] = \
                [[sample, fam_copies, hex_copies, cy5_copies, cy55_copies, rox_copies, atto590_copies],
                 normalized_values]

        # Make sure the dictionary is sorted by wells.  The ddPCR data is grouped by wells.
        # raw_data_dict = dict(sorted(raw_data_dict.items()))

        output_data = ""
        data_counter = 0
        plotting_data_dict = {}

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
                    plotting_data_dict[data_counter] = sample_name, label, normal_value
                    data_counter += 1

        plotting_df_column_labels = ["Sample", "Target", "Data"]

        plotting_df = pd.DataFrame.from_dict(plotting_data_dict, orient='index', columns=plotting_df_column_labels)

        # Remove rows that contain "Genomic", "No RNP", or "No Template" in the "Target" column.  These are not useful
        # for plotting
        plotting_df = plotting_df[~plotting_df["Target"].str.contains("Genomic", case=False, na=False)]
        plotting_df = plotting_df[~plotting_df["Sample"].str.contains(" No RNP", case=False, na=False)]
        plotting_df = plotting_df[~plotting_df["Sample"].str.contains("No Template", case=False, na=False)]

        plot_grouped_histogram(self.ouput_folder, plotting_df)

        outfile = open("{}Normalized_Data.tsv".format(self.ouput_folder), "w")
        outfile.write(output_header + output_data)
        outfile.close()


if __name__ == '__main__':
    DataProcessing(tsv_file = "D:{}ddPCR_Test_Data.tsv".format(os.sep))
