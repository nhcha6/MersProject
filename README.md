# Mers Splicer: User Manual

## Input Files

The program receives two input files:
* An MGF file containing tandem mass spectrometry (MS/MS) data.
* A FASTA file containing protein sequences.

These files must be successfully uploaded before the remaining inputs on the first tab are enabled to the user.

![Input Mgf File](https://raw.githubusercontent.com/arpitbajaj98/MersProject/master/docs/UserManualimgs/inputmgffile.png)   ![Input Fasta File](https://raw.githubusercontent.com/arpitbajaj98/MersProject/master/docs/UserManualimgs/inputfasta.png)

**Figure 1: an example input MGF File (left) and Fasta File (Right).**

### Uploading a file

The files are uploaded by clicking the “Select MGF File” and “Upload FASTA File” buttons. Additional FASTA files can be uploaded by selecting the “Add Another FASTA” button. On click of these buttons, a window will open to select the relevant file. Once selected, the program will either confirm the file upload or inform you of an issues.

If you click the “Select MGF File” or “Upload FASTA File” buttons after having successfully uploaded a file, the program will behave exactly as before. If you successfully upload a second file, it will simply replace the first one. You can add as many additional FASTA files as you please, but currently there is not method of deleting them. If you upload an incorrect FASTA file via the “Add Another FASTA” button you will need to restart close and re-open the program.

The program will check that the files uploaded have the correct extensions, and will ensure that if trans output is selected the input FASTA contains at least two proteins. However, it is unlikely that other miscellaneous errors in the files will be picked up and such inconsistencies are likely to produce bugs. Figure 2 details the process for uploading a Fasta or MGF visually.

![Uploading Files](https://raw.githubusercontent.com/arpitbajaj98/MersProject/master/docs/UserManualimgs/uploadfiles.png)

**Figure 2: infographic showing how to upload a Fasta or MGF file to the program.**

### Protein ID Syntax

The naming style for records in the FASTA file is important. The program returns the name of the protein that a peptide was found in, which it extracts from the protein ID in the input file. If you desire this data, please ensure your protein is named using the format:
				>sp|P61981|1433G_HUMAN 14-3-3 protein gamma……..
The program treats the text between the two | symbols as the name of the protein, and this is what will be included in the output files. For example, the above protein name will be simplified to P61981.

### Protein Sequences and Miscellaneous Characters

Proteins sequences which include characters that do not represent an amino-acid will not stop the program from running. Peptides are created as normal, and if they happen to contain one of these unwanted characters they are simply ignored.

### MGF Spectra Without Charge State

If some of the spectra in the input MGF file are missing the charge-state parameter, they are ignored. Spectra with charge states greater than +6 are similarly ignored.

Are there any other conditions in the input file which are ignored??

### Produce Intensity Plot

Next to the “Upload MGF File” button, there is a checkbox labelled “Produce Intensity Plot”. If this checkbox is ticked when the “Upload MGF File” button is clicked, the program will present the “Intensity Plot” to the user if an MGF file is successfully uploaded. This plot shows the distribution of maximum ion intensities of each spectra. That is, it plots ion intensity on the x-axis and the percentage of spectra in the MGF file which possess at least one ion with a greater intensity value on the y-axis. This plot can be used to choose the desired “Intensity Threshold” (discussed later in the section “First Tab User Inputs”). Figure 3 details how to create an intensity plot visually.

![Intensity Plot](https://raw.githubusercontent.com/arpitbajaj98/MersProject/master/docs/UserManualimgs/mgfplot.png)