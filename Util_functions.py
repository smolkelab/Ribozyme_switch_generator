import time

def find_hairpins(structure):
    '''
    Finds the start and end indices of all the hairpins in a given dotbracket structure.
    :param structure: String denoting secondary structure in dotbracket notation.
    :return: List of lists. Each list contains integers denoting the start and end indices of each hairpin.
    '''

    starts = []
    ends = []
    # Iterates through each bond in a structure.
    for index, bond in enumerate(structure):

        # Every time a forward facing bond is hit, looks to see if a backwards facing bond is the next one. If so, then
        # at a hairpin. Records start and stop indices.
        if bond == '(' and (structure.find(')', index + 1) < structure.find('(', index + 1) or structure.find('(', index + 1) < 0):
            starts.append(index)
            ends.append(structure.find(')', index + 1))

    return [starts, ends]

def get_index_of_bonded(structure, starting_index):
    '''
    Given a dotbracket structure and an index of a nucleotide, finds the index of the nucleotide it is is bonding with.
    The structure must be valid (equal numbers of ( and ), in correct order).
    :param structure: String denoting the structure being evaluated, in dotbracket notation.
    :param starting_index: Integer denoting nucleotide to start from. Must be less than the length of the structure.
    :return: Integer denoting index of bonded nucleotide. Returns -1 for unbonded nucleotides.
    '''

    # Deals with unbonded nucleotides.
    if structure[starting_index] == '.':
        return -1

    # Looks forward through the structure until no more forward looking bonds need a mate.
    if structure[starting_index] == '(':
        bonds_to_go = 0
        for index in range(starting_index + 1, len(structure)):
            if structure[index] == ')' and bonds_to_go == 0:
                return index
            elif structure[index] == '(':
                bonds_to_go += 1
            elif structure[index] == ')':
                bonds_to_go -= 1

    # Looks backward through the structure until no more backward looking bonds need a mate.
    if structure[starting_index] == ')':
        bonds_to_go = 0
        for index in range(starting_index - 1, -1, -1):
            if structure[index] == '(' and bonds_to_go == 0:
                return index
            elif structure[index] == ')':
                bonds_to_go += 1
            elif structure[index] == '(':
                bonds_to_go -= 1

def get_ribozyme_stem_length(sequence, structure, ribozyme_parts):
    '''
    Given a formed ribozyme, gets the length of stem 1 and stem 2.
    :param sequence: String denoting the sequence being evaluated.
    :param structure: String denoting the structure being evaluated, in dotbracket notation.
    :param ribozyme_parts: List of lists containing information on the different parts of the ribozyme. Each list has a
        sequence and structure as a string of the part. Must have 3 parts: A left side that includes the stem of the
        first loop, a top side that includes the stems of both loops and the catalytic core, and a right side that
        includes the stem of the second loop.
    :return: List of lists. Fist list is a list of integers denoting the base lengths of the ribozyme stems. The second
        list contains integers denoting the tested structures' deviation from the base lengths.
    '''

    # Gets the base length of stem 1 and stem 2.
    stem1_length = 0
    # Goes to the start of the loop and works backward until there is no more bonding.
    for five_stem1, three_stem1 in zip(reversed(ribozyme_parts[0][1]), ribozyme_parts[1][1]):
        if (five_stem1 != '(' or three_stem1 != ')'):
            break
        stem1_length += 1

    stem2_length = 0
    for five_stem2, three_stem2 in zip(reversed(ribozyme_parts[1][1]), ribozyme_parts[2][1]):
        if (five_stem2 != '(' or three_stem2 != ')'):
            break
        stem2_length += 1

    stem_lengths = [stem1_length, stem2_length]
    modifications = []

    # Checks each stem for an extended or reduced stem.
    for i in range(2):

        #Starts with a change of zero.
        modification = 0
        loop_start = sequence.find(ribozyme_parts[i][0]) + len(ribozyme_parts[i][0]) - 1
        loop_end = sequence.find(ribozyme_parts[i + 1][0])

        # Checks to see if the stem is reduced.
        if structure[loop_start] == '.':

            # Runs down the stem, checking each nucleotide for a bond.
            for j in range(1, stem_lengths[i]):
                if structure[loop_start - j] == '(' and structure[loop_end + j] == ')':

                    # Makes sure that the bond is to the correct nucleotide.
                    if get_index_of_bonded(structure, loop_start - j) == loop_end + j:
                        modification = -j
                        break
                    else:
                        return [[0, 0], [0, 0]]

                elif not (structure[loop_start - j] == '.' and structure[loop_end + j]) == '.':
                    return [[0, 0], [0, 0]]

        # Checks to see if the stem is extended.
        else:
            added_length = 0

            # Runs up the stem, checking each nucleotide for a bond.
            while True:
                if structure[loop_start + added_length + 1] == '(' and structure[loop_end - added_length - 1] == ')':

                    # Makes sure that the bond is to the correct nucleotide.
                    if get_index_of_bonded(structure, loop_start + added_length) == loop_end - added_length:
                        added_length += 1
                    else:
                        modification = added_length
                        break
                else:
                    modification = added_length
                    break


        modifications.append(modification)

    # Checks that the whole first part of the ribozyme is correct, except for any reduced stem.
    end_index = sequence.find(ribozyme_parts[0][0]) + len(ribozyme_parts[0][0]) - 1
    for i in range(end_index + modifications[0],
                   sequence.find(ribozyme_parts[0][0]) - 1, -1):

        if -(end_index - i + 1) < 0:
            if structure[i] != ribozyme_parts[0][1][-(end_index - i + 1)]:
                return [[0, 0], [0, 0]]

    # Checks that the whole second part of the ribozyme is correct, except for any reduced stem.
    end_index = sequence.find(ribozyme_parts[1][0]) + len(ribozyme_parts[1][0]) - 1
    start_index = sequence.find(ribozyme_parts[1][0]) - 1
    for i in range(end_index + modifications[1],
                   start_index - modifications[0], -1):

        if -(end_index - i + 1) < 0 and end_index - i + 1 < len(ribozyme_parts[1][0]):
            if structure[i] != ribozyme_parts[1][1][-(end_index - i + 1)]:
                return [[0, 0], [0, 0]]

    # Checks that the whole third part of the ribozyme is correct, except for any reduced stem.
    end_index = sequence.find(ribozyme_parts[2][0]) + len(ribozyme_parts[2][0]) - 1
    start_index = sequence.find(ribozyme_parts[2][0]) - 1
    for i in range(end_index,
                   start_index - modifications[1], -1):

        if -(end_index - i + 1) < 0 and end_index - i + 1 < len(ribozyme_parts[2][0]):
            if structure[i] != ribozyme_parts[2][1][-(end_index - i + 1)]:
                return [[0, 0], [0, 0]]

    return [stem_lengths, modifications]

class ProgressBar:
    '''
    Class that can be used to keep track of how far along a process has gotten, and can give a text progress bar and
    estimate the amount of time left of the process.
    '''

    def __init__(self, full_count):
        '''
        Initializes with a count of 0 and an idea of how many times the process will loop before completion. Also
        records the time the process began.
        :param full_count: Integer denoting how many times the process will loop before completion.
        :return: None.
        '''

        self.full_count = full_count
        self.count = 0
        self.start_time = time.time()

    def update(self):
        '''
        Updates the count to reflect more progress.
        :return: None.
        '''

        self.count += 1

    def get_bar(self):
        '''
        Calculates how far the process has gotten to completion and gives a text bar for display.
        :return: String of a text progress bar.
        '''

        # Calculates how many progress bars to use.
        progress = int(50 * (self.count / float(self.full_count)))
        progress_bar = (progress * "=") + ((50 - progress) * "_")
        return "[" + str(progress_bar) + "] " + str(self.count) + "/" + str(self.full_count)

    def get_time_remaining(self):
        '''
        Calculates how much time the process has left by calculating the average time per increment and multiplying that
        by the number of increments left.
        :return: Float value denoting how many seconds remain.
        '''

        time_spent = time.time() - self.start_time

        # Avoids a divide by zero error if called before updating.
        if self.count > 0:
            return time_spent / self.count * (self.full_count - self.count)
        else:
            return time_spent / 1 * (self.full_count - self.count)
