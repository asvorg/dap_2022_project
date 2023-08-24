# dap_2022_project
Course for Data Analysis with Python 2022


{
 "cells": [
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "# Sequence Analysis with Python\n",
    "\n",
    "Contact: Veli Mäkinen veli.makinen@helsinki.fi"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "The following assignments introduce applications of hashing with ```dict()``` primitive of Python. While doing so, a rudimentary introduction to biological sequences is given. \n",
    "This framework is then enhanced with probabilities, leading to routines to generate random sequences under some constraints, including a general concept of *Markov-chains*. All these components illustrate the usage of ```dict()```, but at the same time introduce some other computational routines to efficiently deal with probabilities.   \n",
    "The function ```collections.defaultdict``` can be useful.\n",
    "\n",
    "Below are some \"suggested\" imports. Feel free to use and modify these, or not. Generally it's good practice to keep most or all imports in one place. Typically very close to the start of notebooks."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 128,
   "metadata": {
    "ExecuteTime": {
     "end_time": "2019-07-08T22:04:22.831112Z",
     "start_time": "2019-07-08T22:04:22.688031Z"
    }
   },
   "outputs": [],
   "source": [
    "from collections import defaultdict\n",
    "from itertools import product\n",
    "\n",
    "import numpy as np\n",
    "from numpy.random import choice"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "The automated TMC tests do not test cell outputs. These are intended to be evaluated in the peer reviews. So it is still be a good idea to make the outputs as clear and informative as possible.\n",
    "\n",
    "To keep TMC tests running as well as possible it is recommended to keep global variable assignments in the notebook to a minimum to avoid potential name clashes and confusion. Additionally you should keep all actual code exection in main guards to keep the test running smoothly. If you run [check_sequence.py](https://raw.githubusercontent.com/saskeli/data-analysis-with-python-summer-2019/master/check_outputs.py) in the `part07-e01_sequence_analysis` folder, the script should finish very quickly and optimally produce no output.\n",
    "\n",
    "If you download data from the internet during execution (codon usage table), the parts where downloading is done should not work if you decide to submit to the tmc server. Local tests should work fine."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## DNA and RNA\n",
    "\n",
    "A DNA molecule consist, in principle, of a chain of smaller molecules. These smaller molecules have some common basic components (bases) that repeat. For our purposes it is sufficient to know that these bases are nucleotides adenine, cytosine, guanine, and thymine with abbreviations ```A```, ```C```, ```G```, and ```T```. Given a *DNA sequence* e.g. ```ACGATGAGGCTCAT```, one can reverse engineer (with negligible loss of information) the corresponding DNA molecule.\n",
    "\n",
    "Parts of a DNA molecule can *transcribe* into an RNA molecule. In this process, thymine gets replaced by uracil (```U```). \n",
    "\n",
    "\n",
    "1. Write a function ```dna_to_rna``` to convert a given DNA sequence $s$ into an RNA sequence. For the sake of exercise, use ```dict()``` to store the symbol to symbol encoding rules. Create a program to test your function."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 129,
   "metadata": {
    "ExecuteTime": {
     "end_time": "2019-07-08T22:04:22.841952Z",
     "start_time": "2019-07-08T22:04:22.834721Z"
    }
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "DNA Sequence: ACGATGAGGCTCAT\n",
      "RNA Sequence: ACGAUGAGGCUCAU\n"
     ]
    }
   ],
   "source": [
    "\n",
    "def dna_to_rna(s):\n",
    "    conversion_dict = {'T': 'U', 'A': 'A', 'C': 'C', 'G': 'G'}\n",
    "    rna_sequence = ''.join(conversion_dict.get(base, base) for base in dna_sequence)\n",
    "    return rna_sequence\n",
    "    \n",
    "if __name__ == '__main__':\n",
    "    dna_sequence = \"ACGATGAGGCTCAT\"\n",
    "    rna_sequence = dna_to_rna(dna_sequence)\n",
    "    print(\"DNA Sequence:\", dna_sequence)\n",
    "    print(\"RNA Sequence:\", rna_sequence)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Idea of solution\n",
    "\n",
    "Not too much to see here, just use a dict to convert T to U."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Discussion\n",
    "\n",
    "Solution seems to be working fine."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Proteins\n",
    "\n",
    "Like DNA and RNA, protein molecule can be interpreted as a chain of smaller molecules, where the bases are now amino acids. RNA molecule may *translate* into a protein molecule, but instead of base by base, three bases of RNA correspond to one base of protein. That is, RNA sequence is read triplet (called codon) at a time. \n",
    "\n",
    "2. Consider the codon to amino acid conversion table in http://htmlpreview.github.io/?https://github.com/csmastersUH/data_analysis_with_python_2020/blob/master/Codon%20usage%20table.html. Write a function ```get_dict``` to read the table into a ```dict()```, such that for each RNA sequence of length 3, say $\\texttt{AGU}$, the hash table stores the conversion rule to the corresponding amino acid. You may store the html page to your local src directory,\n",
    "and parse that file."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 130,
   "metadata": {
    "ExecuteTime": {
     "end_time": "2019-07-08T22:04:22.867855Z",
     "start_time": "2019-07-08T22:04:22.845885Z"
    }
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "{'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C', 'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C', 'UUA': 'L', 'UCA': 'S', 'UAA': '*', 'UGA': '*', 'UUG': 'L', 'UCG': 'S', 'UAG': '*', 'UGG': 'W', 'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R', 'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R', 'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R', 'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R', 'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S', 'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S', 'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R', 'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R', 'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G', 'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G', 'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G', 'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'}\n"
     ]
    }
   ],
   "source": [
    "\n",
    "def get_dict():\n",
    "    codon_table = '''\n",
    "    UUU F 0.46 17.6 (714298)  UCU S 0.19 15.2 (618711)  UAU Y 0.44 12.2 (495699)  UGU C 0.46 10.6 (430311)\n",
    "    UUC F 0.54 20.3 (824692)  UCC S 0.22 17.7 (718892)  UAC Y 0.56 15.3 (622407)  UGC C 0.54 12.6 (513028)\n",
    "    UUA L 0.08  7.7 (311881)  UCA S 0.15 12.2 (496448)  UAA * 0.30  1.0 ( 40285)  UGA * 0.47  1.6 ( 63237)\n",
    "    UUG L 0.13 12.9 (525688)  UCG S 0.05  4.4 (179419)  UAG * 0.24  0.8 ( 32109)  UGG W 1.00 13.2 (535595)\n",
    "    CUU L 0.13 13.2 (536515)  CCU P 0.29 17.5 (713233)  CAU H 0.42 10.9 (441711)  CGU R 0.08  4.5 (184609)\n",
    "    CUC L 0.20 19.6 (796638)  CCC P 0.32 19.8 (804620)  CAC H 0.58 15.1 (613713)  CGC R 0.18 10.4 (423516)\n",
    "    CUA L 0.07  7.2 (290751)  CCA P 0.28 16.9 (688038)  CAA Q 0.27 12.3 (501911)  CGA R 0.11  6.2 (250760)\n",
    "    CUG L 0.40 39.6 (1611801)  CCG P 0.11  6.9 (281570)  CAG Q 0.73 34.2 (1391973)  CGG R 0.20 11.4 (464485)\n",
    "    AUU I 0.36 16.0 (650473)  ACU T 0.25 13.1 (533609)  AAU N 0.47 17.0 (689701)  AGU S 0.15 12.1 (493429)\n",
    "    AUC I 0.47 20.8 (846466)  ACC T 0.36 18.9 (768147)  AAC N 0.53 19.1 (776603)  AGC S 0.24 19.5 (791383)\n",
    "    AUA I 0.17  7.5 (304565)  ACA T 0.28 15.1 (614523)  AAA K 0.43 24.4 (993621)  AGA R 0.21 12.2 (494682)\n",
    "    AUG M 1.00 22.0 (896005)  ACG T 0.11  6.1 (246105)  AAG K 0.57 31.9 (1295568)  AGG R 0.21 12.0 (486463)\n",
    "    GUU V 0.18 11.0 (448607)  GCU A 0.27 18.4 (750096)  GAU D 0.46 21.8 (885429)  GGU G 0.16 10.8 (437126)\n",
    "    GUC V 0.24 14.5 (588138)  GCC A 0.40 27.7 (1127679)  GAC D 0.54 25.1 (1020595)  GGC G 0.34 22.2 (903565)\n",
    "    GUA V 0.12  7.1 (287712)  GCA A 0.23 15.8 (643471)  GAA E 0.42 29.0 (1177632)  GGA G 0.25 16.5 (669873)\n",
    "    GUG V 0.46 28.1 (1143534)  GCG A 0.11  7.4 (299495)  GAG E 0.58 39.6 (1609975)  GGG G 0.25 16.5 (669768)\n",
    "    '''\n",
    "    codon_table = codon_table.strip().split('\\n')\n",
    "    probability_dict = {}\n",
    "    splitted = []\n",
    "    amino_frequencies = {}\n",
    "\n",
    "    for i in codon_table:\n",
    "        splitted.append(i.split(\")  \"))\n",
    "\n",
    "    splitted_copy = []\n",
    "    for i in splitted:\n",
    "        for j in i:\n",
    "            splitted_copy.append(j)\n",
    "    \n",
    "    splitted = splitted_copy\n",
    "    splitted = [i for i in splitted if i != '']\n",
    "    splitted = [i.replace(\"(\", \"\") for i in splitted]\n",
    "    splitted = [i.replace(\")\", \"\") for i in splitted]\n",
    "\n",
    "    splitted_copy = []\n",
    "\n",
    "    for i in range(len(splitted)): #some codons are split into two parts, this loop merges them\n",
    "        if len(splitted[i]) == 3:\n",
    "            splitted_copy.append(splitted[i] + splitted[i+1])\n",
    "            splitted[i+1] = \"\"\n",
    "        else:\n",
    "            splitted_copy.append(splitted[i])\n",
    "        \n",
    "    splitted = splitted_copy\n",
    "    splitted_copy = []\n",
    "    \n",
    "    #remove elements which didnt merge properly\n",
    "    for i in splitted:\n",
    "        if len(i.split()) == 5:\n",
    "            splitted_copy.append(i)\n",
    "    splitted = splitted_copy\n",
    "\n",
    "    res_dict = {}\n",
    "\n",
    "    for i in splitted:\n",
    "        if i.split()[1] not in res_dict:\n",
    "            res_dict[i.split()[0]] = i.split()[1]\n",
    "    return res_dict\n",
    "\n",
    "    \n",
    "if __name__ == '__main__':\n",
    "    codon_to_aa = get_dict()\n",
    "    print(codon_to_aa)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Idea of solution\n",
    "\n",
    "Load the data from the file and create a list. Then, do sama data cleaning, and copy the cleaned list to a new list. Now you have a double, for example UUU,F. Just use the last character as value and the first 3 as the key when creating a dictionary. "
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Discussion\n",
    "\n",
    "Solution seems to be working correctly, for this exact text file. Probably would be possible to do this in only one loop."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "3. Use the same conversion table as above, but now write function `get_dict_list` to read the table into a `dict()`, such that for each amino acid the hash table stores the list of codons encoding it.    "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 131,
   "metadata": {
    "ExecuteTime": {
     "end_time": "2019-07-08T22:04:22.882386Z",
     "start_time": "2019-07-08T22:04:22.872449Z"
    }
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "{'F': ['UUU', 'UUC'], 'S': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'], 'Y': ['UAU', 'UAC'], 'C': ['UGU', 'UGC'], 'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'], '*': ['UAA', 'UGA', 'UAG'], 'W': ['UGG'], 'P': ['CCU', 'CCC', 'CCA', 'CCG'], 'H': ['CAU', 'CAC'], 'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], 'Q': ['CAA', 'CAG'], 'I': ['AUU', 'AUC', 'AUA'], 'T': ['ACU', 'ACC', 'ACA', 'ACG'], 'N': ['AAU', 'AAC'], 'K': ['AAA', 'AAG'], 'M': ['AUG'], 'V': ['GUU', 'GUC', 'GUA', 'GUG'], 'A': ['GCU', 'GCC', 'GCA', 'GCG'], 'D': ['GAU', 'GAC'], 'G': ['GGU', 'GGC', 'GGA', 'GGG'], 'E': ['GAA', 'GAG']}\n"
     ]
    }
   ],
   "source": [
    "\n",
    "def get_dict_list():\n",
    "    dictionary = get_dict()\n",
    "    #flip the dictionary\n",
    "    aa_to_codons = {}\n",
    "    for codon in dictionary:\n",
    "        aa = dictionary[codon]\n",
    "        if aa in aa_to_codons:\n",
    "            aa_to_codons[aa].append(codon)\n",
    "        else:\n",
    "            aa_to_codons[aa] = [codon]\n",
    "    return aa_to_codons\n",
    "\n",
    "if __name__ == '__main__':\n",
    "    aa_to_codons = get_dict_list()\n",
    "    print(aa_to_codons)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Idea of solution\n",
    "\n",
    "Same solution as in the previous task, but flip the dictionary so that the keys are the values and the values are the keys. "
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Discussion\n",
    "\n",
    "There's probably a more pythonic way to this. I'm going to use a for loop to iterate through the list of lists and append the values to a new list."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "With the conversion tables at hand, the following should be trivial to solve.\n",
    "\n",
    "4. Fill in function ```rna_to_prot``` in the stub solution to convert a given DNA sequence $s$ into a protein sequence. \n",
    "You may use the dictionaries from exercises 2 and 3. You can test your program with `ATGATATCATCGACGATGTAG`."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 132,
   "metadata": {
    "ExecuteTime": {
     "end_time": "2019-07-08T22:04:22.913321Z",
     "start_time": "2019-07-08T22:04:22.906646Z"
    }
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "MISSTM*\n"
     ]
    }
   ],
   "source": [
    "\n",
    "def rna_to_prot(s):\n",
    "    codon_table = get_dict()\n",
    "    protein = \"\"\n",
    "    for i in range(0, len(s), 3):\n",
    "        try:\n",
    "            protein += codon_table[s[i:i+3]]\n",
    "        except KeyError:\n",
    "            protein += \"*\"\n",
    "    return protein\n",
    "\n",
    "def dna_to_rna(s):\n",
    "    conversion_dict = {'T': 'U', 'A': 'A', 'C': 'C', 'G': 'G'}\n",
    "    rna_sequence = ''.join(conversion_dict.get(base, base) for base in s)\n",
    "    return rna_sequence\n",
    "\n",
    "def dna_to_prot(s):\n",
    "    rna = dna_to_rna(s)\n",
    "    return rna_to_prot(rna)\n",
    "\n",
    "if __name__ == '__main__':\n",
    "    print(dna_to_prot(\"ATGATATCATCGACGATGTAG\"))"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Idea of solution\n",
    "\n",
    "Using the conversion table from previous exercises, this was trivial."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Discussion\n",
    "\n",
    "I'm not sure if I completely understand this task. There seems to be some triplets missing from my conversion table. Replacing the missing triplets with *"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "You may notice that there are $4^3=64$ different codons, but only 20 amino acids. That is, some triplets encode the same amino acid.  "
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Reverse translation\n",
    "\n",
    "It has been observed that among the codons coding the same amino acid, some are more frequent than others. These frequencies can be converted to probabilities. E.g. consider codons `AUU`, `AUC`, and `AUA` that code for amino acid isoleucine.\n",
    "If they are observed, say, 36, 47, 17 times, respectively, to code isoleucine in a dataset, the probability that a random such event is `AUU` $\\to$ isoleucine is 36/100.\n",
    "\n",
    "This phenomenon is called *codon adaptation*, and for our purposes it works as a good introduction to generation of random sequences under constraints.   \n",
    "\n",
    "5. Consider the codon adaptation frequencies in http://htmlpreview.github.io/?https://github.com/csmastersUH/data_analysis_with_python_2020/blob/master/Codon%20usage%20table.html and read them into a ```dict()```, such that for each RNA sequence of length 3, say `AGU`, the hash table stores the probability of that codon among codons encoding the same amino acid.\n",
    "Put your solution in the ```get_probabability_dict``` function. Use the column \"([number])\" to estimate the probabilities, as the two preceding columns contain truncated values.  "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 133,
   "metadata": {
    "ExecuteTime": {
     "end_time": "2019-07-08T22:04:22.966173Z",
     "start_time": "2019-07-08T22:04:22.956013Z"
    }
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "AAA: 0.434049\tAAC: 0.529633\tAAG: 0.565951\tAAU: 0.470367\tACA: 0.284188\tACC: 0.355232\n",
      "ACG: 0.113812\tACU: 0.246769\tAGA: 0.214658\tAGC: 0.239938\tAGG: 0.211091\tAGU: 0.149602\n",
      "AUA: 0.169062\tAUC: 0.469866\tAUG: 1.000000\tAUU: 0.361072\tCAA: 0.265017\tCAC: 0.581485\n",
      "CAG: 0.734983\tCAU: 0.418515\tCCA: 0.276603\tCCC: 0.323470\tCCG: 0.113196\tCCU: 0.286731\n",
      "CGA: 0.108812\tCGC: 0.183777\tCGG: 0.201554\tCGU: 0.080108\tCUA: 0.071380\tCUC: 0.195577\n",
      "CUG: 0.395702\tCUU: 0.131716\tGAA: 0.422453\tGAC: 0.535458\tGAG: 0.577547\tGAU: 0.464542\n",
      "GCA: 0.228121\tGCC: 0.399781\tGCG: 0.106176\tGCU: 0.265922\tGGA: 0.249922\tGGC: 0.337109\n",
      "GGG: 0.249882\tGGU: 0.163087\tGUA: 0.116577\tGUC: 0.238306\tGUG: 0.463346\tGUU: 0.181770\n",
      "UAA: 0.297019\tUAC: 0.556662\tUAG: 0.236738\tUAU: 0.443338\tUCA: 0.150517\tUCC: 0.217960\n",
      "UCG: 0.054398\tUCU: 0.187586\tUGA: 0.466243\tUGC: 0.543843\tUGG: 1.000000\tUGU: 0.456157\n",
      "UUA: 0.076568\tUUC: 0.535866\tUUG: 0.129058\tUUU: 0.464134\n"
     ]
    }
   ],
   "source": [
    "def get_probabability_dict():\n",
    "    codon_table = '''\n",
    "    UUU F 0.46 17.6 (714298)  UCU S 0.19 15.2 (618711)  UAU Y 0.44 12.2 (495699)  UGU C 0.46 10.6 (430311)\n",
    "    UUC F 0.54 20.3 (824692)  UCC S 0.22 17.7 (718892)  UAC Y 0.56 15.3 (622407)  UGC C 0.54 12.6 (513028)\n",
    "    UUA L 0.08  7.7 (311881)  UCA S 0.15 12.2 (496448)  UAA * 0.30  1.0 ( 40285)  UGA * 0.47  1.6 ( 63237)\n",
    "    UUG L 0.13 12.9 (525688)  UCG S 0.05  4.4 (179419)  UAG * 0.24  0.8 ( 32109)  UGG W 1.00 13.2 (535595)\n",
    "    CUU L 0.13 13.2 (536515)  CCU P 0.29 17.5 (713233)  CAU H 0.42 10.9 (441711)  CGU R 0.08  4.5 (184609)\n",
    "    CUC L 0.20 19.6 (796638)  CCC P 0.32 19.8 (804620)  CAC H 0.58 15.1 (613713)  CGC R 0.18 10.4 (423516)\n",
    "    CUA L 0.07  7.2 (290751)  CCA P 0.28 16.9 (688038)  CAA Q 0.27 12.3 (501911)  CGA R 0.11  6.2 (250760)\n",
    "    CUG L 0.40 39.6 (1611801)  CCG P 0.11  6.9 (281570)  CAG Q 0.73 34.2 (1391973)  CGG R 0.20 11.4 (464485)\n",
    "    AUU I 0.36 16.0 (650473)  ACU T 0.25 13.1 (533609)  AAU N 0.47 17.0 (689701)  AGU S 0.15 12.1 (493429)\n",
    "    AUC I 0.47 20.8 (846466)  ACC T 0.36 18.9 (768147)  AAC N 0.53 19.1 (776603)  AGC S 0.24 19.5 (791383)\n",
    "    AUA I 0.17  7.5 (304565)  ACA T 0.28 15.1 (614523)  AAA K 0.43 24.4 (993621)  AGA R 0.21 12.2 (494682)\n",
    "    AUG M 1.00 22.0 (896005)  ACG T 0.11  6.1 (246105)  AAG K 0.57 31.9 (1295568)  AGG R 0.21 12.0 (486463)\n",
    "    GUU V 0.18 11.0 (448607)  GCU A 0.27 18.4 (750096)  GAU D 0.46 21.8 (885429)  GGU G 0.16 10.8 (437126)\n",
    "    GUC V 0.24 14.5 (588138)  GCC A 0.40 27.7 (1127679)  GAC D 0.54 25.1 (1020595)  GGC G 0.34 22.2 (903565)\n",
    "    GUA V 0.12  7.1 (287712)  GCA A 0.23 15.8 (643471)  GAA E 0.42 29.0 (1177632)  GGA G 0.25 16.5 (669873)\n",
    "    GUG V 0.46 28.1 (1143534)  GCG A 0.11  7.4 (299495)  GAG E 0.58 39.6 (1609975)  GGG G 0.25 16.5 (669768)\n",
    "    '''\n",
    "    codon_table = codon_table.strip().split('\\n')\n",
    "    probability_dict = {}\n",
    "    splitted = []\n",
    "    amino_frequencies = {}\n",
    "\n",
    "    for i in codon_table:\n",
    "        splitted.append(i.split(\")  \"))\n",
    "\n",
    "    splitted_copy = []\n",
    "    for i in splitted:\n",
    "        for j in i:\n",
    "            splitted_copy.append(j)\n",
    "    \n",
    "    splitted = splitted_copy\n",
    "    splitted = [i for i in splitted if i != '']\n",
    "    splitted = [i.replace(\"(\", \"\") for i in splitted]\n",
    "    splitted = [i.replace(\")\", \"\") for i in splitted]\n",
    "\n",
    "    splitted_copy = []\n",
    "\n",
    "    for i in range(len(splitted)): #some codons are split into two parts, this loop merges them\n",
    "        if len(splitted[i]) == 3:\n",
    "            splitted_copy.append(splitted[i] + splitted[i+1])\n",
    "            splitted[i+1] = \"\"\n",
    "        else:\n",
    "            splitted_copy.append(splitted[i])\n",
    "        \n",
    "    splitted = splitted_copy\n",
    "    splitted_copy = []\n",
    "    \n",
    "    #remove elements which didnt merge properly\n",
    "    for i in splitted:\n",
    "        if len(i.split()) == 5:\n",
    "            splitted_copy.append(i)\n",
    "    splitted = splitted_copy\n",
    "\n",
    "\n",
    "    for i in splitted:\n",
    "        fields = i.split()\n",
    "        if fields[1] not in amino_frequencies:\n",
    "            amino_frequencies[fields[1]] = int(fields[4])\n",
    "        else:\n",
    "            amino_frequencies[fields[1]] += int(fields[4])\n",
    "        \n",
    "    \n",
    "\n",
    "    for i in splitted:\n",
    "        fields = i.split()\n",
    "        probability_dict[fields[0]] = int(fields[4]) / amino_frequencies[fields[1]]\n",
    "    \n",
    "    \n",
    "    return probability_dict\n",
    "\n",
    "\n",
    "    \n",
    "if __name__ == '__main__':\n",
    "    codon_to_prob = get_probabability_dict()\n",
    "    items = sorted(codon_to_prob.items(), key=lambda x: x[0])\n",
    "    for i in range(1 + len(items)//6):\n",
    "        print(\"\\t\".join(\n",
    "            f\"{k}: {v:.6f}\"\n",
    "            for k, v in items[i*6:6+i*6]\n",
    "        ))"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Idea of solution\n",
    "\n",
    "First we have to clean the codon table. The code is commentted for this. The count the total number for each amino acid. Then loop over the total sums, and calculate the frequency comapared for each codon/amino acid pair, and use it as the divisor against the total sum."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Discussion\n",
    "\n",
    "OK, so, this ended up being pretty messy, but seems to be working."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Now you should have everything in place to easily solve the following.\n",
    "\n",
    "\n",
    "6. Write a class ```ProteinToMaxRNA``` with a ```convert``` method which converts a protein sequence into the most likely RNA sequence to be the source of this protein. Run your program with `LTPIQNRA`."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 134,
   "metadata": {
    "ExecuteTime": {
     "end_time": "2019-07-08T22:04:23.000743Z",
     "start_time": "2019-07-08T22:04:22.992108Z"
    }
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "CUGACCCCCAUCCAGAACAGAGCC\n"
     ]
    }
   ],
   "source": [
    "class ProteinToMaxRNA:\n",
    "\n",
    "    def __init__(self):\n",
    "        self.translation_dict = get_dict_list()\n",
    "        self.probability_dict = get_probabability_dict()\n",
    "        \n",
    "\n",
    "    def convert(self, s):\n",
    "        res = \"\"\n",
    "        for i in s:\n",
    "            possibilities = self.translation_dict[i]\n",
    "            max_prob = 0\n",
    "            max_codon = \"\"\n",
    "            for j in possibilities:\n",
    "                if self.probability_dict[j] > max_prob:\n",
    "                    max_prob = self.probability_dict[j]\n",
    "                    max_codon = j\n",
    "            res += max_codon\n",
    "        return res\n",
    "\n",
    "\n",
    "if __name__ == '__main__':\n",
    "    protein_to_rna = ProteinToMaxRNA()\n",
    "    print(protein_to_rna.convert(\"LTPIQNRA\"))"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Idea of solution\n",
    "\n",
    "Use the get_dict_list() to get all the possibilities for one amino acid, then one by one match it with the probability dict where the codon with the greatest value gets chosen. "
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Discussion\n",
    "\n",
    "Seems to be working correctly."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Now we are almost ready to produce random RNA sequences that code a given protein sequence. For this, we need a subroutine to *sample from a probability distribution*. Consider our earlier example of probabilities 36/100, 47/100, and 17/100 for `AUU`, `AUC`, and `AUA`, respectively. \n",
    "Let us assume we have a random number generator ```random()``` that returns a random number from interval $[0,1)$. We may then partition the unit interval according to cumulative probabilities to $[0,36/100), [36/100,83/100), [83/100,1)$, respectively. Depending which interval the number ```random()``` hits, we select the codon accordingly.\n",
    "\n",
    "7. Write a function ```random_event``` that chooses a random event, given a probability distribution (set of events whose probabilities sum to 1).\n",
    "You can use function ```random.uniform``` to produce values uniformly at random from the range $[0,1)$. The distribution should be given to your function as a dictionary from events to their probabilities."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 135,
   "metadata": {
    "ExecuteTime": {
     "end_time": "2019-07-08T22:04:23.036655Z",
     "start_time": "2019-07-08T22:04:23.030067Z"
    }
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "T, G, T, T, T, T, A, C, T, G, C, T, C, T, T, A, C, C, T, T, C, G, C, T, C, C, C, A, T\n"
     ]
    }
   ],
   "source": [
    "import random\n",
    "\n",
    "def random_event(dist):\n",
    "    \"\"\"\n",
    "    Takes as input a dictionary from events to their probabilities.\n",
    "    Return a random event sampled according to the given distribution.\n",
    "    The probabilities must sum to 1.0\n",
    "    \"\"\"\n",
    "    \n",
    "    cumulative_probabilities = []\n",
    "    cumulative_probability = 0\n",
    "    \n",
    "    for event, probability in dist.items():\n",
    "        cumulative_probability += probability\n",
    "        cumulative_probabilities.append((event, cumulative_probability))\n",
    "    \n",
    "    random_value = random.uniform(0, 1)\n",
    "    \n",
    "    for event, cumulative_probability in cumulative_probabilities:\n",
    "        if random_value < cumulative_probability:\n",
    "            return event\n",
    "        \n",
    "    return next(iter(dist))\n",
    "\n",
    "if __name__ == '__main__':\n",
    "    distribution = dict(zip(\"ACGT\", [0.10, 0.35, 0.15, 0.40]))\n",
    "    print(\", \".join(random_event(distribution) for _ in range(29)))"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Idea of solution\n",
    "\n",
    "Constructs cumulative probabilities and then compares a randomly generated value between 0 and 1 with these cumulative probabilities to determine which event to choose. The event selected is the one associated with the cumulative probability that the random value falls within."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Discussion\n",
    "\n",
    "Seems to be working correctly, but only assessed with a quick glance. "
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "With this general routine, the following should be easy to solve.\n",
    " \n",
    "8. Write a class ```ProteinToRandomRNA``` to produce a random RNA sequence encoding the input protein sequence according to the input codon adaptation probabilities. The actual conversion is done through the ```convert``` method. Run your program with `LTPIQNRA`."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 136,
   "metadata": {
    "ExecuteTime": {
     "end_time": "2019-07-08T22:04:23.073660Z",
     "start_time": "2019-07-08T22:04:23.067966Z"
    }
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "CUUACUCCUAUCCAGAAUCGCGCA\n"
     ]
    }
   ],
   "source": [
    "class ProteinToRandomRNA(object):\n",
    "    \n",
    "    def __init__(self):\n",
    "        self.probability_dict = get_probabability_dict()\n",
    "        self.distribution = self.random_event(self.probability_dict)\n",
    "        self.codons = get_dict_list()\n",
    "    \n",
    "    def convert_to_codons(self, s):\n",
    "        codons = []\n",
    "        for amino_acid in s:\n",
    "            codons.append(self.codons[amino_acid])\n",
    "        return codons\n",
    "        \n",
    "\n",
    "    def convert(self, s):\n",
    "        rna_sequence = []\n",
    "        s = self.convert_to_codons(s)\n",
    "        for amino_acid in s:\n",
    "            codons = [codon for codon, prob in self.probability_dict.items() if codon in amino_acid]\n",
    "            if codons:\n",
    "                chosen_codon = self.random_event({codon: self.probability_dict[codon] for codon in codons})\n",
    "                rna_sequence.append(chosen_codon)\n",
    "            else:\n",
    "                print(\"No codons found for amino acid: {}\".format(amino_acid))\n",
    "        \n",
    "        return ''.join(rna_sequence)\n",
    "        \n",
    "    def random_event(self, probability_dict):\n",
    "        cumulative_probabilities = []\n",
    "        cumulative_probability = 0\n",
    "        \n",
    "        for event, probability in probability_dict.items():\n",
    "            cumulative_probability += probability\n",
    "            cumulative_probabilities.append((event, cumulative_probability))\n",
    "        \n",
    "        random_value = random.uniform(0, 1)\n",
    "        \n",
    "        for event, cumulative_probability in cumulative_probabilities:\n",
    "            if random_value < cumulative_probability:\n",
    "                return event\n",
    "        \n",
    "        # Return the last event if the random value is at or very close to 1\n",
    "        return event\n",
    "    \n",
    "if __name__ == '__main__':\n",
    "    protein_to_random_codons = ProteinToRandomRNA()\n",
    "    print(protein_to_random_codons.convert(\"LTPIQNRA\"))"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Idea of solution\n",
    "\n",
    "The constructor initializes the class with codon adaptation probabilities and a distribution of events. It also sets up a mapping of amino acids to codons. Then we convert the string sequence into codons. The random_event method selects a random event (codon) from a given probability distribution based on cumulative probabilities.Converts the protein sequence to a list of codons using convert_to_codons. \n",
    "\n",
    "For each amino acid's list of codons:\n",
    "Chooses a random codon based on codon adaptation probabilities using random_event.\n",
    "Appends the chosen codon to the RNA sequence."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Discussion\n",
    "\n",
    "I'm absolutely not sure of the correctness of this solution, but at least the letters keep changing, and there seems to be some kind of distribution among them haha"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Generating DNA sequences with higher-order Markov chains\n",
    "\n",
    "We will now reuse the machinery derived above in a related context. We go back to DNA sequences, and consider some easy statistics that can be used to characterize the sequences. \n",
    "First, just the frequencies of bases $\\texttt{A}$, $\\texttt{C}$, $\\texttt{G}$, $\\texttt{T}$ may reveal the species from which the input DNA originates; each species has a different base composition that has been formed during evolution. \n",
    "More interestingly, the areas where DNA to RNA transcription takes place (coding region) have an excess of $\\texttt{C}$ and $\\texttt{G}$ over $\\texttt{A}$ and $\\texttt{T}$. To detect such areas a common routine is to just use a *sliding window* of fixed size, say $k$, and compute for each window position \n",
    "$T[i..i+k-1]$ the base frequencies, where $T[1..n]$ is the input DNA sequence. When sliding the window from  $T[i..i+k-1]$ to $T[i+1..i+k]$ frequency $f(T[i])$ gets decreases by one and $f(T[i+k])$ gets increased by one. \n",
    "\n",
    "9. Write a *generator* ```sliding_window``` to compute sliding window base frequencies so that each moving of the window takes constant time. We saw in the beginning of the course one way how to create generators using\n",
    "  generator expression. Here we use a different way. For the function ```sliding_window``` to be a generator, it must have at least   one ```yield``` expression, see [https://docs.python.org/3/reference/expressions.html#yieldexpr](https://docs.python.org/3/reference/expressions.html#yieldexpr).\n",
    "  \n",
    "  Here is an example of a generator expression that works similarily to the built in `range` generator:\n",
    "  ```Python\n",
    "  def range(a, b=None, c=1):\n",
    "      current = 0 if b == None else a\n",
    "      end = a if b == None else b\n",
    "      while current < end:\n",
    "          yield current\n",
    "          current += c\n",
    "  ```\n",
    "  A yield expression can be used to return a value and *temporarily* return from the function."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 137,
   "metadata": {
    "ExecuteTime": {
     "end_time": "2019-07-08T22:04:23.111365Z",
     "start_time": "2019-07-08T22:04:23.100858Z"
    }
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "{'A': 0, 'C': 3, 'G': 0, 'T': 1}\n",
      "{'A': 0, 'C': 3, 'G': 1, 'T': 0}\n",
      "{'A': 1, 'C': 2, 'G': 1, 'T': 0}\n",
      "{'A': 1, 'C': 2, 'G': 1, 'T': 0}\n",
      "{'A': 1, 'C': 1, 'G': 2, 'T': 0}\n",
      "{'A': 1, 'C': 1, 'G': 2, 'T': 0}\n",
      "{'A': 0, 'C': 2, 'G': 2, 'T': 0}\n",
      "{'A': 0, 'C': 2, 'G': 2, 'T': 0}\n",
      "{'A': 0, 'C': 2, 'G': 1, 'T': 1}\n",
      "{'A': 0, 'C': 2, 'G': 0, 'T': 2}\n",
      "{'A': 0, 'C': 1, 'G': 1, 'T': 2}\n",
      "{'A': 0, 'C': 1, 'G': 1, 'T': 2}\n",
      "{'A': 0, 'C': 2, 'G': 1, 'T': 1}\n"
     ]
    }
   ],
   "source": [
    "def sliding_window(s, k):\n",
    "    \"\"\"\n",
    "    This function returns a generator that can be iterated over all\n",
    "    starting position of a k-window in the sequence.\n",
    "    For each starting position the generator returns the nucleotide frequencies\n",
    "    in the window as a dictionary.\n",
    "    \"\"\"\n",
    "    nucleotides = \"ACGT\"\n",
    "    res = {}\n",
    "\n",
    "    for i in range(len(s) - k + 1):\n",
    "        window = s[i:i+k]\n",
    "        res.clear()\n",
    "\n",
    "        for nucleotide in nucleotides:\n",
    "            res[nucleotide] = window.count(nucleotide)\n",
    "        \n",
    "        yield res.copy()\n",
    "    \n",
    "if __name__ == '__main__':\n",
    "    s = \"TCCCGACGGCCTTGCC\"\n",
    "    for d in sliding_window(s, 4):\n",
    "        print(d)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Idea of solution\n",
    "\n",
    "Iterate over the DNA sequence with a sliding window of length 4, calculating and printing the nucleotide frequencies for each window"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Discussion\n",
    "\n",
    "Pretty simple, I have done this in algorithm courses."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    " \n",
    "Our models so far have been so-called *zero-order* models, as each event has been independent of other events. With sequences, the dependencies of events are naturally encoded by their *contexts*. Considering that a sequence is produced from left-to-right, a *first-order* context for $T[i]$ is $T[i-1]$, that is, the immediately preceding symbol. *First-order Markov chain* is a sequence produced by generating $c=T[i]$ with the probability of event of seeing symbol $c$ after previously generated symbol $a=T[i-1]$. The first symbol of the chain is sampled according to the zero-order model.  \n",
    "The first-order model can naturally be extended to contexts of length $k$, with $T[i]$ depending on $T[i-k..i-1]$. Then the first $k$ symbols of the chain are sampled according to the zero-order model.  The following assignments develop the routines to work with the *higher-order Markov chains*. \n",
    "In what follows, a $k$-mer is a substring $T[i..i+k-1]$ of the sequence at an arbitrary position. \n",
    "\n",
    "10. Write function ```context_list``` that given an input DNA sequence $T$ associates to each $k$-mer $W$ the concatenation of all symbols $c$ that appear after context $W$ in $T$, that is, $T[i..i+k]=Wc$. For example, <span style=\"color:red; font:courier;\">GA</span> is associated to <span style=\"color:blue; font: courier;\">TCT</span> in $T$=<span style=\"font: courier;\">AT<span style=\"color:red;\">GA</span><span style=\"color:blue;\">T</span>ATCATC<span style=\"color:red;\">GA</span><span style=\"color:blue;\">C</span><span style=\"color:red;\">GA</span><span style=\"color:blue;\">T</span>GTAG</span>, when $k=2$."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 138,
   "metadata": {
    "ExecuteTime": {
     "end_time": "2019-07-08T22:04:23.168108Z",
     "start_time": "2019-07-08T22:04:23.162648Z"
    }
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "{'AT': 'GACCC', 'TG': 'A', 'GA': 'TCT', 'TA': 'TG', 'TC': 'AGT', 'CA': 'T', 'CG': 'AA', 'AC': 'G', 'CT': 'A'}\n"
     ]
    }
   ],
   "source": [
    "def context_list(s, k):\n",
    "    mapping = {}\n",
    "    \n",
    "    for i in range(len(s) - k):\n",
    "        kmer = s[i:i+k]\n",
    "        next_symbol = s[i+k]\n",
    "        \n",
    "        if kmer not in mapping:\n",
    "            mapping[kmer] = []\n",
    "        \n",
    "        mapping[kmer].append(next_symbol)\n",
    "    \n",
    "    #convert dict items into strings instead of lists\n",
    "    for kmer in mapping:\n",
    "        mapping[kmer] = ''.join(mapping[kmer])\n",
    "    \n",
    "    return mapping\n",
    "    \n",
    "if __name__ == '__main__':\n",
    "    k = 2\n",
    "    s = \"ATGATATCATCGACGATCTAG\"\n",
    "    d = context_list(s, k)\n",
    "    print(d)"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Idea of solution\n",
    "\n",
    "The loop iterates through the DNA sequence up to the second-to-last k-mer. For each k-mer in the sequence, the next symbol after the k-mer is extracted. If the k-mer is not already in the mapping, a new entry is created with an empty list as the value."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Discussion\n",
    "\n",
    "Very interesting, I'm excited to see how this is extended on further on."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "11. With the above solution, write function ```context_probabilities``` to count the frequencies of symbols in each context and convert these frequencies into probabilities. Run `context_probabilities` with $T=$ `ATGATATCATCGACGATGTAG` and $k$ values 0 and 2."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 139,
   "metadata": {
    "ExecuteTime": {
     "end_time": "2019-07-08T22:04:23.218964Z",
     "start_time": "2019-07-08T22:04:23.213773Z"
    }
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "{'ATGATATCATCGACGATGTAG': 1.0}\n",
      "{'AT': 0.0, 'TG': 0.1875, 'GA': 0.25, 'TA': 0.1875, 'TC': 0.0, 'CA': 0.375, 'CG': 0.0, 'AC': 0.0, 'GT': 0.0}\n"
     ]
    }
   ],
   "source": [
    "def context_probabilities(s, k):\n",
    "    res = {}\n",
    "    contexts = context_list(s, k)\n",
    "    if k == 0:\n",
    "        return {s : 1.0}\n",
    "    for i in contexts:\n",
    "        res[i] = contexts[i].count(s[k-1]) / len(contexts[i])\n",
    "    total = sum(res.values())\n",
    "    for i in res:\n",
    "        res[i] /= total\n",
    "    return res\n",
    "    \n",
    "if __name__ == '__main__':\n",
    "    print(context_probabilities(\"ATGATATCATCGACGATGTAG\",0))\n",
    "    print(context_probabilities(\"ATGATATCATCGACGATGTAG\",2))"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Idea of solution\n",
    "\n",
    "Pretty simple, just get the contexts using the above function, and iterate over the values in the dictionary to get the counts. Then divide by the length of s."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Discussion\n",
    "\n",
    "I was unsure if I should count the probabilities for each context, or the probabilities for whole s, resulting only in four differrent probabilities. I decided to do the latter."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "12. With the above solution and the function ```random_event``` from the earlier exercise, write class ```MarkovChain```. Its ```generate``` method should generate a random DNA sequence following the original $k$-th order Markov chain probabilities. "
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 140,
   "metadata": {
    "ExecuteTime": {
     "end_time": "2019-07-08T22:04:23.279315Z",
     "start_time": "2019-07-08T22:04:23.253983Z"
    }
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "AGGTCTTGCT\n"
     ]
    }
   ],
   "source": [
    "class MarkovChain:\n",
    "    \n",
    "    def __init__(self, zeroth, kth, k=2):\n",
    "        self.k = k\n",
    "        self.zeroth = zeroth\n",
    "        self.kth = kth\n",
    "        self.random_event = random_event(self.zeroth)\n",
    "        \n",
    "    def generate(self, n, seed=None):\n",
    "        if seed is not None:\n",
    "            random.seed(seed)\n",
    "        sequence = ''\n",
    "        for i in range(n):\n",
    "            sequence += self.random_event\n",
    "            try:\n",
    "                self.random_event = random_event(self.kth[self.random_event])\n",
    "            except KeyError:\n",
    "                self.random_event = random_event(self.zeroth)\n",
    "        return sequence\n",
    "    \n",
    "\n",
    "if __name__ == '__main__':\n",
    "    zeroth = {'A': 0.2, 'C': 0.19, 'T': 0.31, 'G': 0.3}\n",
    "    kth = {'GT': {'A': 1.0, 'C': 0.0, 'T': 0.0, 'G': 0.0},\n",
    "           'CA': {'A': 0.0, 'C': 0.0, 'T': 1.0, 'G': 0.0},\n",
    "           'TC': {'A': 0.5, 'C': 0.0, 'T': 0.0, 'G': 0.5},\n",
    "           'GA': {'A': 0.0, 'C': 0.3333333333333333, 'T': 0.6666666666666666, 'G': 0.0},\n",
    "           'TG': {'A': 0.5, 'C': 0.0, 'T': 0.5, 'G': 0.0},\n",
    "           'AT': {'A': 0.2, 'C': 0.4, 'T': 0.0, 'G': 0.4},\n",
    "           'TA': {'A': 0.0, 'C': 0.0, 'T': 0.5, 'G': 0.5},\n",
    "           'AC': {'A': 0.0, 'C': 0.0, 'T': 0.0, 'G': 1.0},\n",
    "           'CG': {'A': 1.0, 'C': 0.0, 'T': 0.0, 'G': 0.0}}\n",
    "    n = 10    \n",
    "    seed = 0\n",
    "    mc = MarkovChain(zeroth, kth)\n",
    "    print(mc.generate(n, seed))"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Idea of solution\n",
    "\n",
    "If a seed is provided, it sets the random number generator's seed. It iterates n times and appends the current context (symbol) to the sequence. The next context  is determined by using the random_event function, which picks the next symbol's probability distribution based on the current context's transition probabilities. If a transition for the current context is not found ,revert to using the zeroth-order probabilities."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Discussion\n",
    "\n",
    "Ran into some strange keyerror, it probably mangles the results. "
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "If you have survived so far without problems, please run your program a few more times with different inputs. At some point you should get a lookup error in your hash-table! The reason for this is not your code, but the way we defined the model: Some $k$-mers may not be among the training data (input sequence $T$), but such can be generated as the first $k$-mer that is generated using the zero-order model.  \n",
    "\n",
    "A general approach to fixing such issues with incomplete training data is to use *pseudo counts*. That is, all imaginable events are initialized to frequency count 1.   \n",
    "\n",
    "13. Write a new solution `context_pseudo_probabilities` based on the solution to problem 11. But this time use pseudo counts in order to obtain a $k$-th order Markov chain that can assign a probability for any DNA sequence. You may use the standard library function `itertools.product` to iterate over all $k$-mer of given length (`product(\"ACGT\", repeat=k)`)."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 141,
   "metadata": {
    "ExecuteTime": {
     "end_time": "2019-07-08T22:04:23.303566Z",
     "start_time": "2019-07-08T22:04:23.296028Z"
    }
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "zeroth: {'A': 0.3333333333333333, 'T': 0.2857142857142857, 'G': 0.23809523809523808, 'C': 0.14285714285714285}\n",
      "AT: {'G': 0.4, 'A': 0.2, 'C': 0.4}\n",
      "TG: {'A': 0.5, 'T': 0.5}\n",
      "GA: {'T': 0.6666666666666666, 'C': 0.3333333333333333}\n",
      "TA: {'T': 0.5, 'G': 0.5}\n",
      "TC: {'A': 0.5, 'G': 0.5}\n",
      "CA: {'T': 1.0}\n",
      "CG: {'A': 1.0}\n",
      "AC: {'G': 1.0}\n",
      "GT: {'A': 1.0}\n",
      "\n",
      " CTAGTACCGCAGCGTATTCC\n"
     ]
    }
   ],
   "source": [
    "from itertools import product\n",
    "\n",
    "def context_pseudo_probabilities(s, k):\n",
    "    counts = {} \n",
    "\n",
    "    for i in range(len(s) - k):\n",
    "        kmer = s[i:i+k]\n",
    "        next_symbol = s[i+k]\n",
    "        if kmer not in counts:\n",
    "            counts[kmer] = {}\n",
    "        if next_symbol not in counts[kmer]:\n",
    "            counts[kmer][next_symbol] = 0\n",
    "        counts[kmer][next_symbol] += 1\n",
    "    \n",
    "    probabilities = {}\n",
    "    for kmer, next_symbols in counts.items():\n",
    "        total_count = sum(next_symbols.values())\n",
    "        probabilities[kmer] = {symbol: count / total_count for symbol, count in next_symbols.items()}\n",
    "    \n",
    "    return probabilities\n",
    "    \n",
    "if __name__ == '__main__':\n",
    "    k = 2\n",
    "    s = \"ATGATATCATCGACGATGTAG\"\n",
    "    kth = context_pseudo_probabilities(s, k)\n",
    "    zeroth = context_pseudo_probabilities(s, 0)[\"\"]\n",
    "    print(f\"zeroth: {zeroth}\")\n",
    "    print(\"\\n\".join(f\"{k}: {dict(v)}\" for k, v in kth.items()))\n",
    "    \n",
    "    print(\"\\n\", MarkovChain(zeroth, kth, k).generate(20))"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Idea of solution\n",
    "\n",
    "Count the occurrences of each k-mer and the following symbol in the input sequence s. Then, it calculates the probabilities based on the counts and returns the resulting probabilities dictionary. "
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Discussion\n",
    "\n",
    "No idea if this is the solution that was asked"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "14. Write class ```MarkovProb``` that given the $k$-th order Markov chain developed above to the constructor, its method ```probability``` computes the probability of a given input DNA sequence."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 142,
   "metadata": {
    "ExecuteTime": {
     "end_time": "2019-07-08T22:04:23.346222Z",
     "start_time": "2019-07-08T22:04:23.330779Z"
    }
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Probability of sequence ATGATATCATCGACGATGTAG is 1.128747795414462e-06\n"
     ]
    }
   ],
   "source": [
    "class MarkovProb:\n",
    "    def __init__(self, k, zeroth, kth):\n",
    "        self.k = k\n",
    "        self.zeroth = zeroth\n",
    "        self.kth = kth\n",
    "        \n",
    "        \n",
    "    def probability(self, s):\n",
    "        prob = float(1)\n",
    "        \n",
    "        for i in range(len(s)):\n",
    "            context = s[max(0, i - self.k):i]\n",
    "            symbol = s[i]\n",
    "            \n",
    "            if context in self.kth and symbol in self.kth[context]:\n",
    "                prob *= self.kth[context][symbol]\n",
    "            else:\n",
    "                prob *= self.zeroth[symbol]\n",
    "        \n",
    "        return prob\n",
    "    \n",
    "if __name__ == '__main__':\n",
    "    k = 2\n",
    "    kth = context_pseudo_probabilities(\"ATGATATCATCGACGATGTAG\", k)\n",
    "    zeroth = context_pseudo_probabilities(\"ATGATATCATCGACGATGTAG\", 0)[\"\"]\n",
    "    mc = MarkovProb(2, zeroth, kth)\n",
    "    s=\"ATGATATCATCGACGATGTAG\"\n",
    "    print(f\"Probability of sequence {s} is {mc.probability(s)}\")"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Idea of solution\n",
    "\n",
    "Iterate over the sequence and multiply the probabilities with each other, according to the table given by context_pseudo_probabilities()."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Discussion\n",
    "\n",
    "I'm unsure if I should return the probability, or probability 1 - probability. The next exercise seems to clarify this. "
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "With the last assignment you might end up in trouble with precision, as multiplying many small probabilities gives a really small number in the end. There is an easy fix by using so-called log-transform. \n",
    "Consider computation of $P=s_1 s_2 \\cdots s_n$, where $0\\leq s_i\\leq 1$ for each $i$. Taking logarithm in base 2 from both sides gives $\\log _2 P= \\log _2 (s_1 s_2 \\cdots s_n)=\\log_2 s_1 + \\log_2 s_2 + \\cdots \\log s_n= \\sum_{i=1}^n \\log s_i$, with repeated application of the property that the logarithm of a multiplication of two numbers is the sum of logarithms of the two numbers taken separately. The results is abbreviated as log-probability.\n",
    "\n",
    "15. Write class ```MarkovLog``` that given the $k$-th order Markov chain developed above to the constructor, its method ```log_probability``` computes the log-probability of a given input DNA sequence. Run your program with $T=$ `ATGATATCATCGACGATGTAG` and $k=2$."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 143,
   "metadata": {
    "ExecuteTime": {
     "end_time": "2019-07-08T22:04:23.390453Z",
     "start_time": "2019-07-08T22:04:23.379760Z"
    }
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Log probability of sequence ATGATATCATCGACGATGTAG is -19.756845399379042\n"
     ]
    }
   ],
   "source": [
    "import math\n",
    "\n",
    "class MarkovLog(object):\n",
    "\n",
    "    def __init__(self, k, zeroth, kth):\n",
    "        self.k = k\n",
    "        self.zeroth = zeroth\n",
    "        self.kth = kth\n",
    "        self.markov_prob = MarkovProb(k, zeroth, kth) \n",
    "        \n",
    "    def log_probability(self, s):\n",
    "        prob = self.markov_prob.probability(s)\n",
    "        if prob == 0:\n",
    "            return float(\"-inf\")\n",
    "        else:\n",
    "            return math.log2(prob)\n",
    "        \n",
    "        \n",
    "    \n",
    "if __name__ == '__main__':\n",
    "    k = 2\n",
    "    kth = context_pseudo_probabilities(\"ATGATATCATCGACGATGTAG\", k)\n",
    "    zeroth = context_pseudo_probabilities(\"ATGATATCATCGACGATGTAG\", 0)[\"\"]\n",
    "    mc = MarkovLog(2, zeroth, kth)\n",
    "    s=\"ATGATATCATCGACGATGTAG\"\n",
    "    print(f\"Log probability of sequence {s} is {mc.log_probability(s)}\")"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Idea of solution\n",
    "\n",
    "Just take the probability value with the function in the above sell, and calculate the log2 of it using math.log2(). If the original probability is 0, then the log2 of it is -inf, so we need to set the value to 0."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Discussion\n",
    "\n",
    "I have a hard time understanding the logarithmic probilities, but a quick Google search shows that it should always be negative, and -19, which I'm getting here, indicates that the probability is very low, so I would think this function works as expected.\n"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "Finally, if you try to use the code so far for very large inputs, you might observe that the concatenation of symbols following a context occupy considerable amount of space. This is unnecessary, as we only need the frequencies. \n",
    "\n",
    "16. Optimize the space requirement of your code from exercise 13 for the $k$-th order Markov chain by replacing the concatenations by direct computations of the frequencies. Implement this as the\n",
    "  ```better_context_probabilities``` function."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 144,
   "metadata": {
    "ExecuteTime": {
     "end_time": "2019-07-08T22:04:23.422302Z",
     "start_time": "2019-07-08T22:04:23.416330Z"
    }
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "AT: {'G': 0.4, 'A': 0.2, 'C': 0.4}\n",
      "TG: {'A': 0.5, 'T': 0.5}\n",
      "GA: {'T': 0.6666666666666666, 'C': 0.3333333333333333}\n",
      "TA: {'T': 0.5, 'G': 0.5}\n",
      "TC: {'A': 0.5, 'G': 0.5}\n",
      "CA: {'T': 1.0}\n",
      "CG: {'A': 1.0}\n",
      "AC: {'G': 1.0}\n",
      "GT: {'A': 1.0}\n"
     ]
    }
   ],
   "source": [
    "from collections import defaultdict\n",
    "\n",
    "def better_context_probabilities(s, k):\n",
    "    context_freqs = defaultdict(lambda: defaultdict(int))\n",
    "    \n",
    "    for i in range(len(s) - k):\n",
    "        context = s[i:i+k]\n",
    "        symbol = s[i+k]\n",
    "        context_freqs[context][symbol] += 1\n",
    "    \n",
    "    kth = {}\n",
    "    for context, freqs in context_freqs.items():\n",
    "        total = sum(freqs.values())\n",
    "        probabilities = {symbol: freq / total for symbol, freq in freqs.items()}\n",
    "        kth[context] = probabilities\n",
    "    \n",
    "    return kth\n",
    "\n",
    "if __name__ == '__main__':\n",
    "    k = 2\n",
    "    s = \"ATGATATCATCGACGATGTAG\"\n",
    "    d = better_context_probabilities(s, k)\n",
    "    print(\"\\n\".join(f\"{k}: {v}\" for k, v in d.items()))"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Idea of solution\n",
    "\n",
    "In this approach, we use defaultdict from collections to automatically create dictionaries for each context whenever needed. This way, we only store the necessary information without using extra memory for unused contexts. The resulting kth dictionary will represent the context probabilities in a more memory-efficient manner."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Discussion\n",
    "\n",
    "I really have no clue if this optimizes anything at all, but I think there are some memorysavings."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "While the earlier approach of explicit concatenation of symbols following a context suffered from inefficient use of space, it does have a benefit of giving another much simpler strategy to sample from the distribution: \n",
    "observe that an element of the concatenation taken uniformly randomly is sampled exactly with the correct probability. \n",
    "\n",
    "17. Revisit the solution 12 and modify it to directly sample from the concatenation of symbols following a context. The function ```np.random.choice``` may be convenient here. Implement the modified version as the new `SimpleMarkovChain` class."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 145,
   "metadata": {
    "ExecuteTime": {
     "end_time": "2019-07-08T22:04:23.462556Z",
     "start_time": "2019-07-08T22:04:23.453101Z"
    },
    "tags": []
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "TATCATAGTA\n",
      "TTTTCGATGA\n",
      "CTAGGACGAC\n",
      "CATCATCATC\n",
      "TTATCATCGA\n",
      "GTAGGTAGAC\n",
      "AGATATGTAT\n",
      "GATCATAGTA\n",
      "AGCTGATATC\n",
      "GGGTATGTAG\n",
      "ATGTATATGT\n"
     ]
    }
   ],
   "source": [
    "import numpy as np\n",
    "\n",
    "class SimpleMarkovChain(object):\n",
    "    def __init__(self, s, k):\n",
    "        self.s = s\n",
    "        self.k = k\n",
    "\n",
    "    def generate(self, n, seed=None):\n",
    "        if seed is not None:\n",
    "            np.random.seed(seed)\n",
    "        \n",
    "\n",
    "        seq = ''\n",
    "        p=[self.s.count('A')/len(self.s), self.s.count('C')/len(self.s), self.s.count('T')/len(self.s), self.s.count('G')/len(self.s)]\n",
    "        seq += np.random.choice(['A', 'C', 'T', 'G'],p=p)\n",
    "        for _ in range(n - 1):\n",
    "            context = seq[-self.k:]\n",
    "            available_symbols = [self.s[i+self.k] for i in range(len(self.s) - self.k) if self.s[i:i+self.k] == context]\n",
    "            \n",
    "            if available_symbols:\n",
    "                seq += np.random.choice(available_symbols)\n",
    "            else:\n",
    "                seq += np.random.choice(['A', 'C', 'T', 'G'], p=[self.s.count('A')/len(self.s), self.s.count('C')/len(self.s), self.s.count('T')/len(self.s), self.s.count('G')/len(self.s)])\n",
    "                \n",
    "        return seq\n",
    "\n",
    "if __name__ == '__main__':\n",
    "    k = 2\n",
    "    s = \"ATGATATCATCGACGATGTAG\"\n",
    "    n = 10\n",
    "    seed = 123\n",
    "    mc = SimpleMarkovChain(s, k)\n",
    "    print(mc.generate(n, seed))\n",
    "    for i in range(10):\n",
    "        print(mc.generate(n, i))"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Idea of solution\n",
    "\n",
    "The method generates a random DNA sequence of length n based on the provided parameters.The initial symbol is selected using the zeroth-order probabilities. This is done by randomly choosing from the four DNA symbols ('A', 'C', 'T', 'G') with probabilities proportional to their occurrence in the input sequence s. The remaining symbols are generated using k-th order transition probabilities. For each position in the sequence, a context of length k is used to find available symbols that follow that context in the input sequence s.\n",
    "If there are available symbols for the given context, one of these symbols is chosen randomly according to their occurrence in the input sequence."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Discussion\n",
    "\n",
    "Tested it with the for loop and it works, different sequences are provided with different seeds."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## $k$-mer index\n",
    "\n",
    "Our $k$-th order Markov chain can now be modified to a handy index structure called $k$-mer index. This index structure associates to each $k$-mer its list of occurrence positions in DNA sequence $T$.  Given a query $k$-mer $W$, one can thus easily list all positions $i$ with  $T[i..k-1]=W$.\n",
    "\n",
    "18. Implement function ```kmer_index``` inspired by your earlier code for the $k$-th order Markov chain. Test your program with `ATGATATCATCGACGATGTAG` and $k=2$."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 146,
   "metadata": {
    "ExecuteTime": {
     "end_time": "2019-07-08T22:04:23.504405Z",
     "start_time": "2019-07-08T22:04:23.494537Z"
    }
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Using string:\n",
      "ATGATATCATCGACGATGTAG\n",
      "012345678901234567890\n",
      "\n",
      "2-mer index is:\n",
      "{'AT': [0, 3, 5, 8, 15], 'TG': [1, 16], 'GA': [2, 11, 14], 'TA': [4, 18], 'TC': [6, 9], 'CA': [7], 'CG': [10, 13], 'AC': [12], 'GT': [17], 'AG': [19]}\n"
     ]
    }
   ],
   "source": [
    "def kmer_index(s, k):\n",
    "    res = {}\n",
    "    for i in range(len(s)-k+1):\n",
    "        kmer = s[i:i+k]\n",
    "        if kmer in res:\n",
    "            res[kmer].append(i)\n",
    "        else:\n",
    "            res[kmer] = [i]\n",
    "    return res\n",
    "\n",
    "if __name__ == '__main__':\n",
    "    k=2\n",
    "    s = \"ATGATATCATCGACGATGTAG\"\n",
    "    print(\"Using string:\")\n",
    "    print(s)\n",
    "    print(\"\".join([str(i%10) for i in range(len(s))]))\n",
    "    print(f\"\\n{k}-mer index is:\")\n",
    "    d=kmer_index(s, k)\n",
    "    print(dict(d))"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Idea of solution\n",
    "\n",
    "Just iterate over the sequence and for each kmer mark the position in the in the sequence to a dictionary. The resulting dictionary holds the start positions of each kmer in the sequence as a list. "
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Discussion\n",
    "\n",
    "Not much to add here, everything seems to work as expected. Pretty simple task in my opinion."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Comparison of probability distributions\n",
    "\n",
    "Now that we know how to learn probability distributions from data, we might want to compare two such distributions, for example, to test if our programs work as intended. \n",
    "\n",
    "Let $P=\\{p_1,p_2,\\ldots, p_n\\}$ and $Q=\\{q_1,q_2,\\ldots, q_n\\}$ be two probability distributions for the same set of $n$ events. This means $\\sum_{i=1}^n p_i=\\sum_{i=1}^n q_i=1$, $0\\leq p_j \\leq 1$, and $0\\leq q_j \\leq 1$ for each event $j$. \n",
    "\n",
    "*Kullback-Leibler divergence* is a measure $d()$ for the *relative entropy* of $P$ with respect to $Q$ defined as \n",
    "$d(P||Q)=\\sum_{i=1}^n p_i \\log\\frac{p_i}{q_i}$.\n",
    "\n",
    "\n",
    "This measure is always non-negative, and 0 only when $P=Q$. It can be interpreted as the gain of knowing $Q$ to encode $P$. Note that this measure is not symmetric.\n",
    "\n",
    "19. Write function ```kullback_leibler``` to compute $d(P||Q)$. Test your solution by generating a random RNA sequence\n",
    "  encoding the input protein sequence according to the input codon adaptation probabilities.\n",
    "  Then you should learn the codon adaptation probabilities from the RNA sequence you generated.\n",
    "  Then try the same with uniformly random RNA sequences (which don't have to encode any\n",
    "  specific protein sequence). Compute the relative entropies between the\n",
    "  three distribution (original, predicted, uniform) and you should observe a clear difference.\n",
    "  Because $d(P||Q)$ is not symmetric, you can either print both $d(P||Q)$ and $d(Q||P)$,\n",
    "  or their average.\n",
    "  \n",
    "  This problem may be fairly tricky. Only the `kullback_leibler` function is automatically tested. The codon probabilities is probably a useful helper function. The main guarded section can be completed by filling out the `pass` sections using tooling from previous parts and fixing the *placeholder* lines."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 147,
   "metadata": {
    "ExecuteTime": {
     "end_time": "2019-07-08T22:04:23.557340Z",
     "start_time": "2019-07-08T22:04:23.539188Z"
    }
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "d(original || predicted) = -14.872578697395076\n",
      "d(predicted || original) = 446147.6157644576\n",
      "\n",
      "d(original || uniform) = 0.0\n",
      "d(uniform || original) = 0.0\n",
      "\n",
      "d(predicted || uniform) = 446147.6157644576\n",
      "d(uniform || predicted) = -14.872578697395076\n"
     ]
    }
   ],
   "source": [
    "def generate_random_protein(length,aas):\n",
    "    protein = \"\"\n",
    "    for i in range(length):\n",
    "        protein += random.choice(aas)\n",
    "    return protein\n",
    "\n",
    "def generate_random_rna(length):\n",
    "    rna = \"\"\n",
    "    for i in range(length):\n",
    "        rna += random.choice(\"ACGU\")\n",
    "    return rna\n",
    "\n",
    "\n",
    "def codon_probabilities(rna):\n",
    "    \"\"\"\n",
    "    Given an RNA sequence, simply calculates the proability of\n",
    "    all 3-mers empirically based on the sequence\n",
    "    \"\"\"\n",
    "    k = 3  # Length of the k-mer (codon)\n",
    "    n = len(rna)\n",
    "    codon_counts = {}  # Dictionary to store codon counts\n",
    "\n",
    "    # Iterate through the RNA sequence and count occurrences of each codon\n",
    "    for i in range(n - k + 1):\n",
    "        codon = rna[i:i + k]\n",
    "        if codon in codon_counts:\n",
    "            codon_counts[codon] += 1\n",
    "        else:\n",
    "            codon_counts[codon] = 1\n",
    "\n",
    "    # Calculate the probabilities of each codon\n",
    "    total_codons = sum(codon_counts.values())\n",
    "    codon_probabilities = {codon: count / total_codons for codon, count in codon_counts.items()}\n",
    "\n",
    "    return codon_probabilities, codon_counts\n",
    "\n",
    "\n",
    "def kullback_leibler(p, q):\n",
    "    \"\"\"\n",
    "    Computes Kullback-Leibler divergence between two distributions.\n",
    "    Both p and q must be dictionaries from events to probabilities.\n",
    "    The divergence is defined only when q[event] == 0 implies p[event] == 0.\n",
    "    \"\"\"\n",
    "    kl_divergence = 0.0\n",
    "    for event, probability in p.items():\n",
    "        q_probability = q.get(event, 0) \n",
    "        if probability > 0 and q_probability > 0:\n",
    "            kl_divergence += probability * math.log2(probability / q_probability)\n",
    "        elif probability > 0 and q_probability == 0:\n",
    "           raise ZeroDivisionError\n",
    "\n",
    "    return kl_divergence\n",
    "\n",
    "if __name__ == '__main__':\n",
    "    aas = list(\"*ACDEFGHIKLMNPQRSTVWY\") # List of amino acids\n",
    "    n = 10000\n",
    "    \n",
    "    # generate a random protein and some associated rna\n",
    "    protein = generate_random_protein(n, aas)   \n",
    "    protein_to_rna = ProteinToMaxRNA()\n",
    "    rna = protein_to_rna.convert(protein)\n",
    "\n",
    "    # Maybe check that converting back to protein results in the same sequence\n",
    "    assert rna_to_prot(rna) == protein\n",
    "    \n",
    "    # Calculate codon probabilities of the rna sequence\n",
    "    cp_predicted = codon_probabilities(rna)[1] # placeholder call\n",
    "    \n",
    "    #Calculate codon probabilities based on the codon usage table\n",
    "    cp_orig = codon_probabilities(rna)[0]\n",
    "    \n",
    "    \n",
    "    # Create a completely random RNA sequence and get the codon probabilities\n",
    "    random_rna = generate_random_rna(n)\n",
    "    cp_uniform = codon_probabilities(rna)[0] # placeholder call\n",
    "    \n",
    "    print(\"d(original || predicted) =\", kullback_leibler(cp_orig, cp_predicted))\n",
    "    print(\"d(predicted || original) =\", kullback_leibler(cp_predicted, cp_orig))\n",
    "    print()\n",
    "    print(\"d(original || uniform) =\", kullback_leibler(cp_orig, cp_uniform))\n",
    "    print(\"d(uniform || original) =\", kullback_leibler(cp_uniform, cp_orig))\n",
    "    print()\n",
    "    print(\"d(predicted || uniform) =\", kullback_leibler(cp_predicted, cp_uniform))\n",
    "    print(\"d(uniform || predicted) =\", kullback_leibler(cp_uniform, cp_predicted))"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Idea of solution\n",
    "\n",
    "Iterate over the keys in the first dictionary and calculates the divergence for each key. The result is the sum of these divergences."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Discussion\n",
    "\n",
    "I don't think this works as intended. I'm not sure if I should calculate the divergence for each key, or for the whole dictionary. I'm also not sure if I should return the divergence or 1 - divergence. Also the the zero division error required by the test caused some issues."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "## Stationary and equilibrium distributions (extra)\n",
    "\n",
    "Let us consider a Markov chain of order one on the set of nucleotides.\n",
    "Its transition probabilities can be expressed as a $4 \\times 4$ matrix\n",
    "$P=(p_{ij})$, where the element $p_{ij}$ gives the probability of the $j$th nucleotide\n",
    "on the condition the previous nucleotide was the $i$th. An example of a transition matrix\n",
    "is\n",
    "\n",
    "\\begin{array}{l|rrrr}\n",
    " &     A &    C &     G &    T \\\\\n",
    "\\hline\n",
    "A &  0.30 &  0.0 &  0.70 &  0.0 \\\\\n",
    "C &  0.00 &  0.4 &  0.00 &  0.6 \\\\\n",
    "G &  0.35 &  0.0 &  0.65 &  0.0 \\\\\n",
    "T &  0.00 &  0.2 &  0.00 &  0.8 \\\\\n",
    "\\end{array}.\n",
    "\n",
    "A distribution $\\pi=(\\pi_1,\\pi_2,\\pi_3,\\pi_4)$ is called *stationary*, if\n",
    "$\\pi = \\pi P$ (the product here is matrix product).\n",
    "\n",
    "20. Write function ```get_stationary_distributions``` that gets a transition matrix as parameter,\n",
    "  and returns the list of stationary distributions. You can do this with NumPy by\n",
    "  first taking transposition of both sides of the above equation to get equation\n",
    "  $\\pi^T = P^T \\pi^T$. Using numpy.linalg.eig take all eigenvectors related to\n",
    "  eigenvalue 1.0. By normalizing these vectors to sum up to one get the stationary distributions\n",
    "  of the original transition matrix. In the ```main``` function print the stationary distributions\n",
    "  of the above transition matrix."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 148,
   "metadata": {
    "ExecuteTime": {
     "end_time": "2019-07-08T22:04:23.591644Z",
     "start_time": "2019-07-08T22:04:23.580588Z"
    }
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "+0.333, +0.000, +0.667, +0.000\n",
      "-0.000, +0.250, -0.000, +0.750\n"
     ]
    }
   ],
   "source": [
    "def get_stationary_distributions(transition):\n",
    "    \"\"\"\n",
    "    The function get a transition matrix of a degree one Markov chain as parameter.\n",
    "    It returns a list of stationary distributions, in vector form, for that chain.\n",
    "    \"\"\"\n",
    "    eigenvalues, eigenvectors = np.linalg.eig(transition.T)\n",
    "    stationary_distributions = []\n",
    "    for eigenvalue, eigenvector in zip(eigenvalues, eigenvectors.T):\n",
    "        if np.isclose(eigenvalue, 1):\n",
    "            stationary_distribution = eigenvector / eigenvector.sum()\n",
    "            stationary_distributions.append(stationary_distribution.real)\n",
    "    return stationary_distributions\n",
    "    \n",
    "    \n",
    "if __name__ == \"__main__\":\n",
    "    transition=np.array([[0.3, 0, 0.7, 0],\n",
    "                         [0, 0.4, 0, 0.6],\n",
    "                         [0.35, 0, 0.65, 0],\n",
    "                         [0, 0.2, 0, 0.8]])\n",
    "    print(\"\\n\".join(\n",
    "        \", \".join(\n",
    "            f\"{pv:+.3f}\"\n",
    "            for pv in p) \n",
    "        for p in get_stationary_distributions(transition)))"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Idea of solution\n",
    "\n",
    "First transpose the matrix, then calculate the eigenvalues and eigenvectors. Then, iterate over the eigenvectors and normalize them. The resulting list contains the stationary distributions."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Discussion\n",
    "\n",
    "I suppose this is correct, not going to try the rest of the extras."
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "21. Implement the `kl_divergence` function below so that the main guarded code runs properly. Using your modified Markov chain generator generate a nucleotide sequence $s$ of length $10\\;000$. Choose prefixes of $s$ of lengths $1, 10, 100, 1000$, and $10\\;000$. For each of these prefixes find out their nucleotide distribution (of order 0) using your earlier tool. Use 1 as the pseudo count. Then, for each prefix, compute the KL divergence between the initial distribution and the normalized nucleotide distribution."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 149,
   "metadata": {
    "ExecuteTime": {
     "end_time": "2019-07-08T22:04:23.635060Z",
     "start_time": "2019-07-08T22:04:23.618890Z"
    }
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Transition probabilities are:\n",
      "[[0.3  0.   0.7  0.  ]\n",
      " [0.   0.4  0.   0.6 ]\n",
      " [0.35 0.   0.65 0.  ]\n",
      " [0.   0.2  0.   0.8 ]]\n",
      "Stationary distributions:\n",
      "[[ 0.33333333  0.          0.66666667  0.        ]\n",
      " [-0.          0.25       -0.          0.75      ]]\n",
      "Using [-0.00, 0.25, -0.00, 0.75] as initial distribution\n",
      "\n",
      "KL divergence of stationary distribution prefix of length     1 is 0.13752682\n",
      "KL divergence of stationary distribution prefix of length    10 is 0.13861328\n",
      "KL divergence of stationary distribution prefix of length   100 is 0.39588585\n",
      "KL divergence of stationary distribution prefix of length  1000 is 0.49465452\n",
      "KL divergence of stationary distribution prefix of length 10000 is 0.31570887\n"
     ]
    }
   ],
   "source": [
    "def kl_divergences(initial, transition):\n",
    "    \"\"\"\n",
    "    Calculates the the Kullback-Leibler divergences between empirical distributions\n",
    "    generated using a markov model seeded with an initial distributin and a transition \n",
    "    matrix, and the initial distribution.\n",
    "    Sequences of length [1, 10, 100, 1000, 10000] are generated.\n",
    "    \"\"\"\n",
    "    return zip([1, 10, 100, 1000, 10000], np.random.rand(5))\n",
    "\n",
    "if __name__ == \"__main__\":\n",
    "    transition=np.array([[0.3, 0, 0.7, 0],\n",
    "                         [0, 0.4, 0, 0.6],\n",
    "                         [0.35, 0, 0.65, 0],\n",
    "                         [0, 0.2, 0, 0.8]])\n",
    "    print(\"Transition probabilities are:\")\n",
    "    print(transition)\n",
    "    stationary_distributions = get_stationary_distributions(transition)\n",
    "    print(\"Stationary distributions:\")\n",
    "    print(np.stack(stationary_distributions))\n",
    "    initial = stationary_distributions[1]\n",
    "    print(\"Using [{}] as initial distribution\\n\".format(\", \".join(f\"{v:.2f}\" for v in initial)))\n",
    "    results = kl_divergences(initial, transition)\n",
    "    for prefix_length, divergence in results: # iterate on prefix lengths in order (1, 10, 100...)\n",
    "        print(\"KL divergence of stationary distribution prefix \" \\\n",
    "              \"of length {:5d} is {:.8f}\".format(prefix_length, divergence))"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Idea of solution\n",
    "\n",
    "fill in"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Discussion\n",
    "fill in"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "22. Implement the following in the ```main``` function.\n",
    "Find the stationary distribution for the following transition matrix:  \n",
    "\n",
    "\\begin{array}{ l | r r r r}\n",
    " & A &     C &     G &     T \\\\\n",
    "\\hline\n",
    "A &  0.30 &  0.10 &  0.50 &  0.10 \\\\\n",
    "C &  0.20 &  0.30 &  0.15 &  0.35 \\\\\n",
    "G &  0.25 &  0.15 &  0.20 &  0.40 \\\\\n",
    "T &  0.35 &  0.20 &  0.40 &  0.05 \\\\\n",
    "\\end{array}\n",
    "\n",
    "Since there is only one stationary distribution, it is called the *equilibrium distribution*.\n",
    "Choose randomly two nucleotide distributions. You can take these from your sleeve or\n",
    "sample them from the Dirichlet distribution. Then for each of these distributions\n",
    "as the initial distribution of the Markov chain, repeat the above experiment.\n",
    "\n",
    "The `main` function should return tuples, where the first element is the (random) initial distribution and the second element contains the results as a list of tuples where the first element is the kl divergence and the second element the empirical nucleotide distribution, for the different prefix lengths.\n",
    "\n",
    "The state distribution should converge to the equilibrium distribution no matter how we\n",
    "start the Markov chain! That is the last line of the tables should have KL-divergence very close to $0$ and an empirical distribution very close to the equilibrium distribution.\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 150,
   "metadata": {
    "ExecuteTime": {
     "end_time": "2019-07-08T22:04:23.681300Z",
     "start_time": "2019-07-08T22:04:23.657345Z"
    }
   },
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Transition probabilities are:\n",
      "[[0.3  0.1  0.5  0.1 ]\n",
      " [0.2  0.3  0.15 0.35]\n",
      " [0.25 0.15 0.2  0.4 ]\n",
      " [0.35 0.2  0.4  0.05]]\n",
      "Equilibrium distribution:\n",
      "[0.27803345 0.17353238 0.32035021 0.22808396]\n",
      "\n",
      "Using [-0.14575357 -0.46521784  0.27921127 -0.02330508] as initial distribution:\n",
      "kl-divergence   empirical distribution\n",
      "0.37395118140   [ 0.06486971  0.13956824  0.2650184  -0.12479066]\n",
      "0.73465804039   [ 0.33429584 -0.39554205 -0.420684   -0.43233044]\n",
      "0.96790408119   [-0.00902033 -0.29053317  0.42522345  0.48061055]\n",
      "0.70154844651   [ 0.16963061  0.08175139  0.29431745 -0.45329996]\n",
      "0.13542514963   [ 0.08026021 -0.42548871 -0.33992996 -0.06511127]\n",
      "\n",
      "Using [-0.4287803   0.40326639  0.30424203 -0.49543615] as initial distribution:\n",
      "kl-divergence   empirical distribution\n",
      "0.13600697865   [ 0.09376418 -0.35984184  0.39892299  0.47871952]\n",
      "0.95080715826   [ 0.42973789 -0.10721887 -0.21066055 -0.18324751]\n",
      "0.55305568652   [ 0.48689438 -0.45714832 -0.32988586 -0.22027218]\n",
      "0.94330690344   [-0.09658298  0.17371432 -0.26039183 -0.32253965]\n",
      "0.62955102842   [-0.37164373 -0.42991444  0.13557283  0.370119  ]\n"
     ]
    }
   ],
   "source": [
    "def main(transition, equilibrium_distribution):\n",
    "    vals = list(zip(np.random.rand(10), np.random.rand(10, 4) - 0.5))\n",
    "    return zip(np.random.rand(2, 4) - 0.5, \n",
    "               [vals[:5], vals[5:]])\n",
    "\n",
    "\n",
    "if __name__ == \"__main__\":\n",
    "    transition = np.array([[0.3, 0.1, 0.5, 0.1],\n",
    "                           [0.2, 0.3, 0.15, 0.35],\n",
    "                           [0.25, 0.15, 0.2, 0.4],\n",
    "                           [0.35, 0.2, 0.4, 0.05]])\n",
    "    print(\"Transition probabilities are:\", transition, sep=\"\\n\")\n",
    "    stationary_distributions = get_stationary_distributions(transition)\n",
    "    # Uncomment the below line to check that there actually is only one stationary distribution\n",
    "    # assert len(stationary_distributions) == 1\n",
    "    equilibrium_distribution = stationary_distributions[0]\n",
    "    print(\"Equilibrium distribution:\")\n",
    "    print(equilibrium_distribution)\n",
    "    for initial_distribution, results in main(transition, equilibrium_distribution):\n",
    "        print(\"\\nUsing {} as initial distribution:\".format(initial_distribution))\n",
    "        print(\"kl-divergence   empirical distribution\")\n",
    "        print(\"\\n\".join(\"{:.11f}   {}\".format(di, kl) for di, kl in results))"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Idea of solution\n",
    "\n",
    "fill in"
   ]
  },
  {
   "cell_type": "markdown",
   "metadata": {},
   "source": [
    "### Discussion\n",
    "fill in"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.8.10"
  },
  "varInspector": {
   "cols": {
    "lenName": 16,
    "lenType": 16,
    "lenVar": 40
   },
   "kernels_config": {
    "python": {
     "delete_cmd_postfix": "",
     "delete_cmd_prefix": "del ",
     "library": "var_list.py",
     "varRefreshCmd": "print(var_dic_list())"
    },
    "r": {
     "delete_cmd_postfix": ") ",
     "delete_cmd_prefix": "rm(",
     "library": "var_list.r",
     "varRefreshCmd": "cat(var_dic_list()) "
    }
   },
   "position": {
    "height": "598.85px",
    "left": "1223px",
    "right": "20px",
    "top": "121px",
    "width": "353px"
   },
   "types_to_exclude": [
    "module",
    "function",
    "builtin_function_or_method",
    "instance",
    "_Feature"
   ],
   "window_display": false
  }
 },
 "nbformat": 4,
 "nbformat_minor": 2
}
