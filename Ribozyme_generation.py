import os
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import Util_functions
import os.path

def RNAStructure_minimal_generator(sequence, structures, text, temp = 310):
    '''
    Creates dot bracket files from RNA sequences that represent minimum free energy structures. Saves them as .txt
    files.
    :param sequence: Input RNA sequence
    :param structures: Number of minimum free energy structures desired
    :param text: Name that you wish to save the .txt files under
    :return: Returns .fasta files based on sequence, .ct and dot bracket fies for structures. Saves all in same folder
    '''

    fasta = SeqRecord(Seq(sequence),"temp","temp")
    fasta_list = [fasta]
    SeqIO.write(fasta_list,text+".fasta","fasta")
    os.system("Fold '"+text+".fasta' '"+text+".ct' -w 3 -p 100 -t " + str(temp))

    # Generates dotbracket file for each structure in desired number of structures.
    for i in range(1, structures+1):
        os.system("ct2dot '"+text+".ct' "+str(i)+" '"+text+"."+str(i)+".txt'")

def RNAStructure_get_reference_structures(sequence, type, left_ribozyme = '', temp = 310):
    '''
    Gets the reference structure for the RNAStructure program. Can be switched to get ribozyme or aptamer. Checks to see
    if hanging ends should be cut from aptamer or loops cut from ribozyme.
    :param sequence: String denoting the sequence being evaluated.
    :param type: String denoting whether evaluating an aptamer or a ribozyme.
    :param left_ribozyme: String denoting the sequence of the left side of the ribozyme. This side goes all the way up
        to the tip of loop 2, where the aptamer is normally added. Does not have to be accurate if the loops are going
        to be cut. Set to empty by default.
    :return: Returns list containing sequence and structure of aptamer, or list containing list of lists of ribozyme
        parts and list containing loop sequences an structure.
    '''
    if not os.path.exists("Reference"):
        os.makedirs("Reference")
    os.chdir("Reference")

    RNAStructure_minimal_generator(sequence, 1, "reference", temp)
    # Gets ribozyme parts structures (five_ribostruct, three_ribostruct)
    fileopen = open("reference.1.txt")
    structure = fileopen.readlines()[2][:-1]

    os.chdir('..')
    shutil.rmtree("Reference")

    print('RNAStructure reference ' + type + ' structure: ')

    print(sequence)
    print(structure)

    # Makes sure the reference structure is the desired one.
    if input("Is this correct? (y/n) ") == 'n':
        structure = input("Please enter correct dotbracket structure: ")

        while len(structure) != len(sequence):
            structure = input("Length does not match structural length, please enter correct structure: ")

        while len(structure.replace("(","").replace(")","").replace("|","").replace(".","")) > 0:
            structure = input("Non dotbracket notation entered, please enter correct structure: ")

    # Gets aptamer information.
    if type == 'aptamer':

        return cut_aptamer_hanging(sequence, structure)

    # Gets ribozyme information.
    if type == 'ribozyme':

        ribozyme_parts = cut_ribozyme_loops(sequence, structure, left_ribozyme)
        [loops, stem_lengths] = get_ribozyme_loops(sequence, structure, ribozyme_parts)
        return [ribozyme_parts, loops]

def cut_ribozyme_loops(sequence, structure, left_ribozyme):
    '''
    If desired, as determined by user input, removes the hairpins from the sTRSV ribozyme from the reference structure.
    If this is done, then any structure in the loop regions will be counted as a legitimate structure as long as all the
    stems are formed. If not done, then the ribozyme loops must match exactly the input structure to count as a formed
    ribozyme. Will look for a hydrogen bond at the top of loop 2 instead of a covalent bond.
    :param sequence: String denoting the sequence of the reference ribozyme.
    :param structure: Structure denoting the structure of the reference ribozyme. Only hairpins should be in the
        standard loops 1 and 2 as the natural sTRSV ribozyme.
    :param left_ribozyme: String denoting the sequence of the left side of the ribozyme. This side goes all the way up
        to the tip of loop 2, where the aptamer is normally added. Does not have to be accurate if the loops are going
        to be cut.
    :return: List of lists. Each list contains strings for sequence and structure for part of the ribozyme. If loops are
        cut, then will have a left side that includes the stem of the first loop, a top side that includes the stems of
        both loops and the catalytic core, and a right side that includes the stem of the second loop. If loops are left
        in, then will have a left and right part, split at the middle of loop 2.
    '''

    # Gets user input as to whether loops should be left in or cut out.
    cut_loops = input('Remove loops from reference structure? Will not look for formed loops to count as a ribozyme. '
                      'Will still look for stems. y/n ')

    if cut_loops == 'y':

        # If cut out, defines the loop parts as ending on the stem before and after each loop.
        [starts, ends] = Util_functions.find_hairpins(structure)

        left = [sequence[0:starts[0] + 1], structure[0:starts[0] + 1]]
        top = [sequence[ends[0]:starts[1] + 1], structure[ends[0]:starts[1] + 1]]
        right = [sequence[ends[1]:], structure[ends[1]:]]

        labels = ['left', 'top', 'right']

        # Prints out the new reference for each side for user knowledge.
        for index, label in enumerate(labels):
            print('Reference for the ' + label + ' side:')
            print([left, top, right][index][0])
            print([left, top, right][index][1])

        return [left, top, right]

    # If not cut out, define the part as left and right side, split at the middle of loop 2.
    left = [sequence[:len(left_ribozyme)], structure[:len(left_ribozyme)-2] + ".("]
    right = [sequence[len(left_ribozyme):], ")." + structure[len(left_ribozyme)+2:]]

    labels = ['left', 'right']

    # Prints out the new reference for each side for user knowledge.
    for index, label in enumerate(labels):
        print('Reference for the ' + label + ' side:')
        print([left, right][index][0])
        print([left, right][index][1])

    return [left, right]

def get_ribozyme_loops(sequence, structure, ribozyme_parts):
    '''
    Given a formed ribozyme, gets the sequence and structure of loop 1 and loop 2, as well as the lengths of stems
    leading to those loops.
    :param sequence: String denoting the sequence being evaluated.
    :param structure: String denoting the structure being evaluated, in dotbracket notation.
    :param ribozyme_parts: List of lists containing information on the different parts of the ribozyme. Each list has a
        sequence and structure as a string of the part. Must have 3 parts: A left side that includes the stem of the
        first loop, a top side that includes the stems of both loops and the catalytic core, and a right side that
        includes the stem of the second loop.
    :return: List of lists. Fist list is list of strings, one for each loop. Each string has, in alternating order the
        nucleotide at each position moving 5' to 3' and the structure of the nucleotide at that position as a dot or
        bracket. Defines the loop as the first continuous stretch of unbonded nucleotides following the stem from the
        ribozyme. Can have bonded nucleotides if they form another stem leading off the oop, but the loop returned here
        will not have anything past this loop. Second list contains the length of stems leading to each loop. Returns
        empty loop string and stem length of 0 if loop is not found.

        Example: The below structure would only return the loop closest to the ribozyme.
         ___       ___
        /   \_____/   \_____
        |    _|_|_|    _|_|_ RIBOZYME
        \___/     \___/
    '''
    # Gets information on stem lengths.
    [base_stem_lengths, length_modifications] = Util_functions.get_ribozyme_stem_length(sequence, structure, ribozyme_parts)

    if base_stem_lengths[0] + length_modifications[0] == 0 and base_stem_lengths[1] + length_modifications[1] == 0:
        return[['', ''], [0, 0]]

    # Finds the beginning and end of each loop.
    loop1_indices = [sequence.find(ribozyme_parts[0][0]) + len(ribozyme_parts[0][0]) + length_modifications[0],
             sequence.find(ribozyme_parts[1][0]) - length_modifications[0]]

    loop2_indices = [sequence.find(ribozyme_parts[1][0]) + len(ribozyme_parts[1][0]) + length_modifications[1],
             sequence.find(ribozyme_parts[2][0]) - length_modifications[1]]

    loop1 = ''
    loop2 = ''
    loops = [loop1, loop2]
    out_loops = []
    # Iterates through each loop.
    for in_loop, loop_indices in zip(loops, [loop1_indices, loop2_indices]):

        # Iterates through each nucleotide in the loop, recording sequence and structure.
        next_index = -1
        for index in range(loop_indices[0], loop_indices[1]):
            if index >= next_index:
                in_loop += sequence[index] + structure[index]

                # Stops counting if hits a bond, which indicates another stem coming off the loop. Waits until the stem
                # comes back into the loop to keep counting again.
                if structure[index] == '(':
                    next_index = Util_functions.get_index_of_bonded(structure, index)

        out_loops.append(in_loop)

    # Returns loop strings and overall stem length.
    return [out_loops, [base_stem_lengths[0] + length_modifications[0], base_stem_lengths[1] + length_modifications[1]]]

def cut_aptamer_hanging(sequence, structure):
    '''
    Removes unbonded nucleotides from the 5' or 3' end of the aptamer sequence from consideration for a formed aptamer.
    Asks user if removal is desired.
    :param sequence: String denoting the sequence of the aptamer.
    :param structure: String denoting the structure of the aptamer, in dotbracket notation.
    :return: List of strings denoting the output sequence and structure.
    '''

    # Checks to see if either end is unbonded.
    if structure[0] == '.' or structure[-1] == '.':

            # Asks for input as to whether hanging ends should eb removed.
            hanging = input('Remove hanging ends? Will not look for ends to be unbonded to count as aptamer '
                                'formed. y/n ' )
            if hanging == 'y':

                # Changes structure and sequence to remove unbonded groups at the end. Prints out new references for
                # the user.
                structure = structure[structure.find('('):structure.rfind(')')+1]
                sequence = sequence[structure.find('('):structure.rfind(')')+1]
                print('New reference structure:')
                print(sequence)
                print(structure)

    return [sequence, structure]
