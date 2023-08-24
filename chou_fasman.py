import pandas as pd  # import pandas library for data manipulation and analysis

def pre_process(table_file, sequence_file):
    """
    This function reads in a CSV file containing scoring information for protein sequences,
    as well as a text file containing a protein sequence. It returns the sequence, the scoring
    information in the form of a pandas dataframe, and the length of the sequence.

    Parameters:
    - table_file: (str) file name of the CSV file containing scoring information
    - sequence_file: (str) file name of the text file containing the protein sequence
    
    Returns:
    - (tuple) sequence (str), pandas dataframe, N (int)
    """
    try:  # attempt to read the csv file using pandas library
        cf_table = pd.read_csv(table_file)
    except FileNotFoundError:  # handle file not found error
        print(f"Error: File '{table_file}' not found in the opened directory!")
        exit(1)  # exit the program with an error code of 1

    try:  # attempt to read the sequence file
        with open(sequence_file, 'r') as f:
            seq = f.read()  # read the content of the file and store it in a variable
    except FileNotFoundError:  # handle file not found error
        print(f"Error: File '{sequence_file}' not found in the opened directory!")
        exit(1)  # exit the program with an error code of 1

    N = len(seq)  # get the length of the sequence

    return seq, cf_table, N  # return the sequence, the pandas dataframe, and the length of the sequence

def extension(struct, strt, stp, cf_table, X, who):
    """
    Extend a secondary structure window by four to the left and right to include a range of residues with a specific score.

    Parameters:
    - struct (list): The secondary structure to be modified.
    - strt (int): The starting index of the range to be modified.
    - stp (int): The ending index (exclusive) of the range to be modified.
    - cf_table (pandas.DataFrame): A table of contact frequencies for each amino acid pair.
    - X (str): The character to be used for the extension (e.g. 'S' for sheet, 'H' for helix).
    - who (str): The column in `cf_table` to use for scoring the residues ('Pa' for helix, 'Pb' for sheet).

    Returns:
    The modified `struct` list.
    """
    # Set the values in the range [strt, stp) which is either a 5 or 6 window of the `struct` list to `X` depending of helix of sheet
    struct[strt:stp] = [X]*(stp - strt)

    # Initialize variables for scoring and string indices
    lft_scr = 4
    rgt_scr = 4
    #scores intialised to 4 so that they can iterate inside the while loop given the condition, set to 0 right after loop starts
    rgt_ptr = stp - 3
    lft_ptr = strt + 3

    # Extend to the right until the score is less than 4 or the end of the sequence is reached
    while(rgt_scr >= 4 and rgt_ptr + 4 < N ):
        rgt_scr = 0
        # Compute the score for the next 4 characters in the sequence
        for j in seq[rgt_ptr:rgt_ptr + 4]:
            rgt_scr += float(cf_table.loc[cf_table['code'] == j, who])
        # If the score is greater than or equal to 4, set the next character in `struct` to `X`
        if(rgt_scr >= 4):
            struct[rgt_ptr + 3] = X
        rgt_ptr += 1

    # Extend to the left until the score is less than 4 or the beginning of the sequence is reached
    while(lft_scr >= 4 and lft_ptr - 4 > 0):
        lft_scr = 0
        # Compute the score for the previous 4 characters in the sequence
        for j in seq[lft_ptr - 4:lft_ptr]:
            lft_scr += float(cf_table.loc[cf_table['code'] == j, who])
        # If the score is greater than or equal to 4, set the previous character in `struct` to `X`
        if(lft_scr >= 4):
            struct[lft_ptr - 4] = X
        lft_ptr += 1

    # Return the modified `struct` list
    return struct

def alpha_chk(N, seq, cf_table,helix):
    '''
    This function checks for alpha helixes in a protein sequence and extends them if they meet the criteria.
    
    Arguments:
    - N (int): The length of the protein sequence
    - seq (str): The protein sequence to be checked
    - cf_table (pd.DataFrame): A Pandas DataFrame containing amino acid properties
    - helix (list): A list representing the current state of the protein structure
    
    Returns:
    - str: A string representing the state of the protein structure after checking for and extending alpha helixes
    '''
    strt = 0 # starting index of current six-letter subsequence
    stp = 6 # ending index of current six-letter subsequence
    while(stp <= N):
        score = 0 # score of current six-letter subsequence
        for i in seq[strt:stp]:
            score += min(int(cf_table.loc[cf_table['code'] == i, 'Pa']),1) # adding score of each amino acid to the total score
        if(score >=4): # if the score is greater than or equal to 4, extend the helix
            extension(helix, strt, stp, cf_table, 'H', 'Pa')
        # move to the next six-letter subsequence
        strt += 1
        stp += 1
    return ''.join(helix) # return the final protein structure as a string


def beta_chk(N, seq, cf_table,sheet):
    """
    This function checks for the presence of beta sheets in the protein sequence and updates the sheet structure accordingly.

    Arguments:
    - N (int): The length of the protein sequence.
    - seq (str): The protein sequence.
    - cf_table (pd.DataFrame): A pandas dataframe containing the conformational free energies for each amino acid.
    - sheet (List[str]): A list of characters representing the secondary structure of the protein. Each character can be either 'H' (for helix), 'S' (for sheet), or '_' (for no defined structure).

    Returns:
    - str: The updated sheet structure of the protein.
    """
    # Set the start and end positions for the sliding window
    strt = 0
    stp = 5
    
    # Slide the window along the protein sequence
    while(stp <= N):
        # Initialize the score to 0
        score = 0
        # Check the score for each amino acid in the window
        for i in seq[strt:stp]:
            # Add the minimum of the conformational free energy for beta sheets or 1
            score += min(int(cf_table.loc[cf_table['code'] == i, 'Pb']),1)
        # If the score meets the threshold, update the sheet structure
        if(score >=3):
            extension(sheet, strt, stp, cf_table, 'S', 'Pb')
        # Slide the window by one position
        strt += 1
        stp += 1
    
    # Return the updated sheet structure
    return ''.join(sheet)

def conflict_resolution(N, seq,helix, sheet, cf_table, scnd_struct):
    """
    Resolves conflicts between predicted secondary structures (helix or sheet) at each residue in a protein sequence.

    Arguments:
    - N (int): length of protein sequence
    - seq (str): protein sequence
    - helix (str): predicted helix secondary structure for each residue
    - sheet (str): predicted sheet secondary structure for each residue
    - cf_table (pandas DataFrame): table of amino acid pair conflict scores
    - scnd_struct (str): list of secondary structure predictions after conflict resolution

    Returns:
    - scnd_struct (str): list of secondary structure predictions after conflict resolution
    """

    i = 0
    while(i < N):
        #store the starting value in case i shifted due to constinuous string of conflicting structures
        j = i
        score_a = 0
        score_b = 0
        while(i < N and helix[i] == 'H' and sheet[i] == 'S'):
            # calculate conflict scores for helix and sheet
            score_a += float(cf_table.loc[cf_table['code'] == seq[i], 'Pa'])
            score_b += float(cf_table.loc[cf_table['code'] == seq[i], 'Pb'])
            i += 1

        # choose secondary structure with the lower conflict score
        cnflct = 'S' if score_a < score_b else 'H'

        # assign the chosen secondary structure to all residues in the conflicting region
        for f in range(j,i):
            scnd_struct[f] = cnflct

        # if there are no conflicts, assign the secondary structure predicted by helix or sheet
        if(i < N):
            scnd_struct[i] = 'H' if helix[i] != '_' else ('S' if sheet[i] != '_' else '_')
            i += 1

    return ''.join(scnd_struct)


# Main Code

# Call the pre_process function to read in the Chou-Fasman table and protein sequence from files
seq, cf_table, N = pre_process("ChouFas.csv", 'Protein_seq.txt')

# Create empty lists of length N to store predicted secondary structure using the alpha and beta check methods
predict_a =['_']*N
predict_b =['_']*N

# Create an empty list of length N to store the final predicted secondary structure
predict = ['_']*N



# Use the alpha check method to predict helical regions and store the predictions in predict_a
alpha_chk(N,seq,cf_table,predict_a)

# Use the beta check method to predict sheet regions and store the predictions in predict_b
beta_chk(N,seq,cf_table,predict_b)

# Use the conflict resolution method to combine the predictions from predict_a and predict_b and store the final predictions in predict
conflict_resolution(N, seq, predict_a, predict_b, cf_table, predict)

# Print the protein sequence in a easily viewable 50 string sliced output
i = 0
while(i<N):
    print()
    end = min(i+50, N)
    info_str = str(i)+':'+str(end)
    spces = ' '*(20-len(info_str))
    info_str += spces
    print(info_str, end = '')
    for j in range(i, end):
        print(seq[j],end = '')
    print(end = "\n"+' '*20)
    for j in range(i, end):
        print(predict[j],end = '')
    i+=50
    print()
print()
