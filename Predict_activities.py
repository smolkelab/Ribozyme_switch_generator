import pickle
import numpy as np
from keras.models import model_from_json
import csv

def loop_one_hot_encode(loop_seq, loop_struct):
    '''
    Encodes a ribozyme loop into a 3-dimensional array meant to represent its 3-domensional physical form. The
    dimensions are as follows:
        width - 2, to represent the nucleotides as the directly emerge from a paired stem
        height - 15, to encapsulate long loops. Most positions will be blank for shorter loops
        channels - 8, to represent the 4 nucleotides, each in a bonded or unbonded state
    :param loop_seq: String of nucleotides representing the sequence of the loop
    :param loop_struct: String of dot-bracket notation representing the structure of the loop
    :return: numpy array of size (2, 1, 15, 8)
    '''

    encoded_loop = np.zeros((2, 1, 15, 8))
    loop_codes = ['A', 'A', 'U', 'U', 'C', 'C', 'G', 'G']

    if len(loop_seq) <= 30:
        # Iterates forward from the 5' end of the loop and backwards from the 3' end simultaneously, marking at each
        # position the correct channel.
        for i in range(len(loop_seq) // 2):
            encoded_loop[0, 0, i, loop_codes.index(loop_seq[i]) + (loop_struct[i] == '(' or loop_struct[i] == ')')] = 1
            encoded_loop[1, 0, i, loop_codes.index(loop_seq[-(i + 1)]) + (loop_struct[-(i + 1)] == '(' or loop_struct[-(i + 1)] == ')')] = 1

        # If an odd number of nucleotides, put the on at the apex of the loop into the 1st column.
        if len(loop_seq) % 2 == 1:
            i = len(loop_seq) // 2
            encoded_loop[0, 0, i, loop_codes.index(loop_seq[i]) + (loop_struct[i] == '(' or loop_struct[i] == ')')] = 1

    else:
        # If the loop is too large to fit, just does the first 15 on each side.
        for i in range(15):
            encoded_loop[0, 0, i, loop_codes.index(loop_seq[i]) + (loop_struct[i] == '(' or loop_struct[i] == ')')] = 1
            encoded_loop[1, 0, i, loop_codes.index(loop_seq[-(i + 1)]) + (loop_struct[-(i + 1)] == '(' or loop_struct[-(i + 1)] == ')')] = 1

    return encoded_loop

def struct_dict_to_array(in_dict):
    '''
    Turns a dictionary of ribozyme sequences and structures into a 5-dimensional numpy array, containing 4-dimensional
    arrays representing paired ribozyme loops for each sequence, along with an array containing the basal
    gene-regulatory activity for each sequence.
    :param in_dict: Dictionary of ribozyme sequences where the keys are tuples like so: (sequence string, dot-bracket
    structure string) and the value is the basal gene-regulatory activity for that sequence.
    :return: Tuple containing 3 arrays:
        5-dimensional numpy array representing paired ribozyme loops for each sequence
        1-dimensional numpy array containing the basal gene-regulatory activity for each sequence
        1-dimensional array of tuples, containing the sequences of the 2 loops for each sequence
    '''

    out_X = []
    out_y = []
    out_loops = []
    first, mid, last = 'GCUGUC', 'CUGAUGA', 'GAAACAGC'

    # Iterates through each sequence in the dictionary
    for seq in in_dict.keys():
        test_seq = seq[0]
        test_struct = seq[1]

        # Gets the sequence and array of the first loop
        loop_seq = test_seq[test_seq.find(first) + len(first): test_seq.find(mid)]
        loop_struct = test_struct[test_seq.find(first) + len(first): test_seq.find(mid)]
        l1 = loop_one_hot_encode(loop_seq, loop_struct)
        l1_seq = loop_seq

        # Gets the sequence and array of the second loop
        loop_seq = test_seq[test_seq.find(mid) + len(mid): test_seq.find(last)]
        loop_struct = test_struct[test_seq.find(mid) + len(mid): test_seq.find(last)]
        l2 = loop_one_hot_encode(loop_seq, loop_struct)

        # Pairs the 2 loops up and stores paired array and value
        out_X.append(np.concatenate((l1, l2), axis=1))
        out_y.append(in_dict[seq][0])
        out_loops.append((l1_seq, loop_seq))

    out_X = np.array(out_X).astype('float32')
    out_y = np.expand_dims(out_y, axis=1)

    return out_X, out_y, out_loops

# Reads in test data
test_list = pickle.load(open('Candidate_list_RNAs_min_structures.pkl', 'rb'))
rest_dict = {}
test_dict = {}

for seq in test_list:
    target_dict = test_dict
    target_dict[seq[0]] = seq

test_seq_dict = {}
# Pulls test data into dictionary for later conversion to array
for i in test_dict.keys():
    test_seq_dict[(tuple(test_dict[i][1][0]), tuple(test_dict[i][1][1]), test_dict[i][0], test_dict[i][2])] = \
        [1]

# Segments training data by structure, creating a dictionaries where each sequence in that dictionary has the same
# structure
test_segmented_dict = {}
for loop in test_seq_dict:
    if (len(loop[0][0]) / 2, len(loop[0][1]) / 2, loop[1][0], loop[1][1]) not in test_segmented_dict:
        test_segmented_dict[(int(len(loop[0][0]) / 2), int(len(loop[0][1]) / 2), loop[1][0], loop[1][1])] = {}

    test_segmented_dict[(len(loop[0][0]) / 2, len(loop[0][1]) / 2, loop[1][0], loop[1][1])][(loop[2], loop[3])] = \
    test_seq_dict[loop]

# List of structures in the segmented dictionary that have over 50 entries
test_populous_loops = [i for i in test_segmented_dict.keys() if len(test_segmented_dict[i]) > 1]

all_pr = []
all_loops = []
# For each structure segment, finds the appropriate model, pulls it, and gets predictions for sequences in that segment
for te_seg in test_segmented_dict:
    try:
        teX, teY, teloops = struct_dict_to_array(test_segmented_dict[te_seg])

        json_file = open('Models/' + str(list(te_seg)) + ".json", 'r')
        loaded_model_json = json_file.read()
        json_file.close()
        loaded_model = model_from_json(loaded_model_json)
        loaded_model.load_weights('Models/' + str(list(te_seg)) + "model.h5")

        pr = loaded_model.predict(teX, batch_size=32)

        all_loops.extend(teloops)
        all_pr.extend(pr)
        print("Model for " + str(te_seg) + " found and used.")

    except:
        print("Model for " + str(te_seg) + " not found.")

# Writes the predicted values out to a csv in order of lowest predicted basal gene-regulatory activity to highest.
with open('predictions.csv', 'w', newline='') as csvfile:
    writer = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(['Loop I seq', 'Loop II seq', 'Predicted basal log10(GFP/mCh)'])
    best_pr = [i[0] for i in sorted(enumerate(all_pr), key=lambda x:x[1])]
    for i in best_pr:
        writer.writerow([all_loops[i][0], all_loops[i][1], all_pr[i][0]])