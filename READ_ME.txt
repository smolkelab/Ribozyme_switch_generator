Developer: Calvin Schmidt
Questions to: calvinschmdt@gmail.com.

Setup:
If on windows:
	Download VirtualBox and install a linux virtual machine.
	Step-by-step instructions: http://www.psychocats.net/ubuntu/virtualbox

Programs needed:
Python 3
Modules:
    Biopython - http://biopython.org/ for installation instructions, for manipulating sequences.
    Tensorflow, Keras - keras.io

For RNAstructure:
1. Download RNAstructure Linux Text Interface from this site: http://rna.urmc.rochester.edu/
2. Use the command line to navigate into the folder containing the downloaded .tgz file. <cd Foldername> changes the current folder.
3. Unpack the file using the command <tar -zxvf Filename.tgz>.
4. Navigate to the RNAstructure folder.
5. Build the package using the command <make all>.

Ribozyme generation:
    1. Run Generate_candidate_list.py. Make sure to enter integers for prompted loops sizes and that the aptamer sequence is all caps and contains only A, U, C, and G.
        - This generates a .pkl file with the list of all ribozyme sequences using this aptamer.
    2. Enter the command <export PATH=$PATH:~/Desktop/RNAstructure/exe/> where the path is the path to the exe folder in the built RNAstructure program.
    3. Enter the command <export DATAPATH=~/Desktop/RNAstructure/data_tables/> with the path replaced here as well.
    4. Run Fold_candidate_list.py
        - This generates a .pkl file with all the folded and analyzed ribozyme sequences.
    5. Make sure the ribozyme structures and aptamer structures are accurate. Getting rid of the ribozyme loops enables more flexible tracking of ribozyme formation.
    6. Run Predict_activities.py. Make sure all the models are being loaded in and used.
        - This generates a .csv file with the loop sequences and predicted basal gene-regulatory activity for each sequence.

    Tips:
    Each N added increases processing time by 5x. 6-7 Ns can be finished overnight depending on the complexity of the aptamer, context, and programs desired.
        - The parameter finder attempts to predict how long it will take to run the library. It is accurate within an order of magnitude.
    Context may disrupt predictions, try both with and without these.


