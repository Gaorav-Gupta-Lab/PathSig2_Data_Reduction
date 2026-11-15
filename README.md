# PathSig2_Data_Reduction
### Normalization and Plotting of PathSig2 Data

## Sample Naming Convention

The sample naming convention used in PathSig2 data reduction tool is as follows:

- The sample name is extracted from the second column of the input data.
- Duplicate samples, either biological or technical, must have a space then a single letter at the end.  The last two characters of the sample name are removed.
- The no template or water control must be named "No Template"If the sample name is "No Template", it is skipped and not processed further.