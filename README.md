# Chou-Fasman Secondary Structure Prediction

This README provides an overview of the Chou-Fasman algorithm and its implementation for predicting the secondary structure of a protein sequence.

## Project Overview

This project implements the Chou-Fasman algorithm, a method for predicting the secondary structure of a protein sequence. The algorithm utilizes a combination of amino acid propensities to form alpha-helices and beta-sheets.

### Code Files

- `chou_fasman.py`: Contains the implementation of the Chou-Fasman algorithm.
- `ChouFas.csv`: CSV file containing prediction parameters for the algorithm.
- `Protein_seq.txt`: Text file containing the input protein sequence.

## Output

When the code in `chou_fasman.py` is executed correctly, it returns the following output:

```plaintext
0:50  SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDTVYCPRHVICTAEDML _HHHHHHHHHHHSSSSSSSSSSSSSSHHHHHHHHHSSSSHHHHHHHHHHH
50:100 NPNYEDLLIRKSNHSFLVQAGNVQLRVIGHSMQNCLLRLKVDTSNPKTPK H__HHHHHHHHHHSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS____HS
100:150 YKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNHTIKGSFLNGSCGSVGF SSSSSSSSSSSSSSSSSSSS__SSSSSSSSSS___HHHHHH_________
