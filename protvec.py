import pandas as pd


# read protvec file
protvec_file = "protVec_100d_3grams.csv"
protvec_data = pd.read_csv(protvec_file, sep="\t", header=0,skiprows=1)

# dictionary to hold all 3-grams
protvec = {}


rows_protvec = protvec_data.shape[0]


for row in range(0,1): #TODO: change it back to rows_protvec to get everything
    dict_key = protvec_data.iloc[[row],[0]]
    #TODO: parsing the right type and data
    dict_key = dict_key[3:6]
    dict_value = []
    for column in range(1,101):
        dict_value.append(protvec_data.iloc[[row],[column]])
    protvec[dict_key] = dict_value


# protvec["AAA"]
#protvec.iloc[[0],[0]]


# test sequence
seq = "MSRFSIEGKSLKLDAITTEDEKSVFAVLLEDDSVKEIVLSGNTIGTEAARWLSENIASKKDLEIAEFSDIFTGRVKDEIPEALRLLLQALLKCPKLHTVRLSDNAFGPTAQEPLIDFLSKHTPLEHLYLHNNGLGPQAGAKIARALQELAVNKKAKNAPPLRSIICGRNRLENGSMKEWAKTFQSHRLLHTVKMVQNGIRPEGIEHLLLEGLAYCQELKVLDLQDNTFTHLGSSALAIALKSWPNLRELGLNDCLLSARGAAAVVDAFSKLENIGLQTLRLQYNEIELDAVRTLKTVIDEKMPDLLFLELNGNRFSEEDDVVDEIREVFSTRGRGELDELDDMEELTDEEEEDEEEEAESQSPEPETSEEEKEDKELADELSKAHI"

# creating a list which is then filled with 100 zero values
sequence_vec = []
for i in range(0,100):
    sequence_vec.append(0)


# iterate overall overlapping 3-grams (e.g. MAK, AKL, KLP, LPQ,...)
for i in range(0, len(seq)-2, 1):
    threegram = seq[i:i+3]
    #prot_vec_3_gram = protvec[threegram]
    sequence_vec = sequence_vec #+ prot_vec_3_gram

print(sequence_vec)




# playground

#list = [[[1,2],2],[3,4]]
#dic = {"three" : [1,2]}
#dic["three"]
#print(list[[0]])