import pickle
import numpy as np
from keras.models import Sequential
from keras.models import model_from_json
from keras.layers import Dense, Dropout, Flatten
from keras.layers import Conv3D

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

# Reads in training data
rs_list = []
rs_list.extend(pickle.load(open('NGS_data/CMS_NGS1_5_RNAs_min_structures_all_activities.pkl', 'rb'), encoding = 'latin1'))
rs_list.extend(pickle.load(open('NGS_data/CMS_NGS2_1_RNAs_min_structures_all_activities_normed.pkl', 'rb'), encoding = 'latin1'))
rs_list.extend(pickle.load(open('NGS_data/CMS_NGS3_1_RNAs_min_structures_all_activities.pkl', 'rb'), encoding = 'latin1'))

rest_dict = {}
# Pulls training data into dictionary for later conversion to array
for seq in rs_list:
    target_dict = rest_dict
    target_dict[seq[0][0][0]] = seq

train_seq_dict = {}
train_dict = rest_dict
for i in train_dict.keys():
    train_seq_dict[(tuple(train_dict[i][0][0][1][0]), tuple(train_dict[i][0][0][1][1]), train_dict[i][0][0][0],
                    train_dict[i][0][0][2])] = \
        [train_dict[i][0][1], train_dict[i][1]]

# Segments training data by structure, creating a dictionaries where each sequence in that dictionary has the same
# structure
train_segmented_dict = {}
for loop in train_seq_dict:
    if (len(loop[0][0]) / 2, len(loop[0][1]) / 2, loop[1][0], loop[1][1]) not in train_segmented_dict:
        train_segmented_dict[(len(loop[0][0]) / 2, len(loop[0][1]) / 2, loop[1][0], loop[1][1])] = {}

    train_segmented_dict[(len(loop[0][0]) / 2, len(loop[0][1]) / 2, loop[1][0], loop[1][1])][(loop[2], loop[3])] = \
    train_seq_dict[loop]

# List of structures in the segmented dictionary that have over 50 entries
train_populous_loops = [i for i in train_segmented_dict.keys() if len(train_segmented_dict[i]) > 50]

# List detailing how to expand the structure search parameter if not enough sequences in training set to train model
diff_list = [(1, 1, 1, 1), (2, 1, 1, 1), (2, 2, 1, 1), (3, 2, 1, 1), (3, 3, 1, 1), (3, 3, 2, 1), (3, 3, 2, 2),
             (4, 4, 2, 2)]

# Get the loop size of model to save
te_seg = [0, 0, 0, 0]
te_seg[0] = int(input("Loop 1 size: "))
te_seg[1] = int(input("Loop 2 size: "))
te_seg[2] = int(input("Stem 1 length: "))
te_seg[3] = int(input("Stem 2 length: "))


training_dict = {}
# Iterates through the diff_list, relaxing the requirements for structural similarity until 1000 sequences are in the
# training set
for i in diff_list:
    # Updates the dictionary with sequences that have the prescribed difference in structure
    for tr_seg in train_populous_loops:
        if abs(te_seg[0] - tr_seg[0]) < i[0] and abs(te_seg[1] - tr_seg[1]) < i[1] and \
                abs(te_seg[2] - tr_seg[2]) < i[2] and abs(te_seg[3] - tr_seg[3]) < i[3]:
            training_dict.update(train_segmented_dict[tr_seg])

    # If at least 1000 sequences are in the training set, moves forward
    trX, trY, trloops = struct_dict_to_array(training_dict)
    if len(trY) > 1000:
        break

# Defines the model
model = Sequential()
layer = Conv3D(32, (2, 2, 2),
               activation='relu',
               input_shape=(2, 2, 15, 8))
model.add(layer)
model.add(Dropout(0.25))
model.add(Flatten())
model.add(Dense(64, activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(output_dim=1))

# Optimize with SGD
model.compile(loss='mean_squared_error', optimizer='adam')

# Fit model in batches
model.fit(trX, trY[:, 0], nb_epoch=100, batch_size=1000, verbose=1)

# Dumps model to json for later predictions
model_json = model.to_json()
with open('Models/' + str(te_seg) + ".json", "w") as json_file:
    json_file.write(model_json)
model.save_weights('Models/' + str(te_seg) + "model.h5")
print("Saved model to disk")