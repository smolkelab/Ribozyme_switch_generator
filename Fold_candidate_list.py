import os
import Ribozyme_generation
import Util_functions
import pickle
import shutil

# Load in the list of sequences to fold
struct_file = open('seq_list.pkl', 'rb')
full_list = pickle.load(struct_file)

# Get the structure of the native ribozyme for comparison
five_HHRz = 'GCUGUCACCGGAUGUGCUUUCCGGUCUGAUGAGUCCGU'
three_HHRz = 'GAGGACGAAACAGC'
[ribozyme_parts, loops] = Ribozyme_generation.RNAStructure_get_reference_structures(five_HHRz + three_HHRz,
                                                                                 'ribozyme', five_HHRz)

print([ribozyme_parts, loops])

tuple_list = []
bar = Util_functions.ProgressBar(len(full_list))

for seq in full_list:

    if not os.path.exists("Test_ribozymes"):
        os.makedirs("Test_ribozymes")

    Ribozyme_generation.RNAStructure_minimal_generator(seq, 1, "Test_ribozymes/test")

    # Iterates through each tested ribozyme structure(teststruct) and finds ribozyme active and aptamer formed.
    testopen = open("Test_ribozymes/test.1.txt")
    testopen = testopen.readlines()

    # Gets the sequence of the loops for the sequence, if correctly folded.
    teststruct = ''
    if testopen != []:
        teststruct = testopen[2][:-1]
        [loops, stem_lengths] = Ribozyme_generation.get_ribozyme_loops(seq, teststruct, ribozyme_parts)

    else:
        [loops, stem_lengths] = [['', ''], [0, 0]]

    shutil.rmtree("Test_ribozymes")

    tuple_list.append((seq, [loops, stem_lengths], teststruct))

    bar.update()
    print(bar.get_bar())
    print(bar.get_time_remaining())

# Dumps list of sequences, folded structure, and loop sequences to a pickle file for storage and later analysis.
struct_file = open('Candidate_list_RNAs_min_structures.pkl', 'wb')
pickle.dump(tuple_list, struct_file)
struct_file.close()