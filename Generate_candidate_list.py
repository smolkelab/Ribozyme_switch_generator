import pickle

def to_base_5(n):
    '''
    Converts number to its base 4 equivalent. Used to generate all possible combinations of nucleotides.
    :param i: Current iteration
    :return: base 4 number
    '''
    s=[]
    while n:
        s.append(str(n%5))
        n=n//5
    return ''.join(s[::-1])

def N_replace(num_N, iteration):
    '''
    Replaces Ns with specific nucleotides based on which iteration it is.
    :param num_N: The number of Ns being changed.
    :param iteration: The current iteration.
    :return: Returns nucleotides to fill in place of Ns
    '''
    new_N = to_base_5(iteration)
    num = 0
    # Replaces each number with corresponding base
    for base in ['S','A','U','C','G']:
        new_N = new_N.replace(str(num), base)
        num += 1
    return new_N

# Get the upper and lower bounds on the lengths of the random loop
low_N = int(input("Smallest loop size to consider: "))
high_N = int(input("Largest loop size to consider: "))

# Create the list of random loop sequences
loop_list = []
for i in range(5 ** high_N):
    loop = N_replace(high_N, i)
    if 'S' not in loop and len(loop) > low_N - 1:
        loop_list.append(loop)

apt = input("Aptamer sequence: ")
five_HHRz = 'GCUGUCACCGGA'
mid_HHRz = 'UCCGGUCUGAUGAGUCC'
three_HHRz = 'GGACGAAACAGC'
five_insulator = 'GGGAAACAAACAAA'
three_insulator = 'AAAAAGAAAAAUAAAAA'

# Add the aptamer and random loop sequences onto each of the two ribozyme loops to create list of candidate sequences
rbz_list = []
for j in loop_list:
    rbz_list.append(five_insulator + five_HHRz + j + mid_HHRz + apt + three_HHRz + three_insulator)
    rbz_list.append(five_insulator + five_HHRz + apt + mid_HHRz + j + three_HHRz + three_insulator)

# Dump the list as a pickle file
struct_file = open('seq_list.pkl', 'wb')
pickle.dump(rbz_list, struct_file)
struct_file.close()