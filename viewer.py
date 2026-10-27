from mendeleev import element
import sys

non_metal_series = [ #series that make up the non metals (used to exclude metals)
    "Noble gases",
    "Halogens",
    "Nonmetals", #(not currently including metaloids due to diffrent rules)
]

smiles_input = input("please input a molecule in SMILES format: -Note metals are not currently supported in this version")
elements_included = [] #list of element symbols in order of mention in smiles input
valences = [] #list of element valences (in alignment order to "elements_included")
for i in range(len(smiles_input)):
    if smiles_input[i].isalpha() and smiles_input[i].isupper(): #if current string position is alphabetical and uppercase (signaling start of a periodic symbol)
        if len(smiles_input) > 1: #if there actually are more than 1 chars
            if smiles_input[i+1].isalpha() and smiles_input[i+1].islower(): #if next character is again alphabetical but lowercase (meaning a 2 char periodic symbol)
                symbol = smiles_input[i] + smiles_input[i+1] #combine strings to form full periodic symbol
                elements_included.append(symbol)
            else: #if only a 1 char symbol
                symbol = smiles_input[i]
                elements_included.append(symbol) #only append the 1 char as it is the full symbol
        else: #if there is only 1 char input
            print("please include more than one element")
            sys.exit()

print(elements_included)
for i in elements_included: #warn about metals before exit
    try: #catch error at very first instance of caling data about element (indicating no such element exists)
        if element(i).series not in non_metal_series: #check element is non metal
            print("Metals where detected in your SMILES input -please only use non metals (currently not including metaloids)")
            sys.exit()
        elif element(i).series == "Noble gases": #check element is not a noble gas
            print("Noble gases where detected in your SMILES input -these do not form bonds")
            sys.exit()
    except:
        print(f"please check you correctly spelled or entered your elements correctly, or that they actually exist: {i} is not recognised")
        sys.exit()

for i in elements_included: #number of valence electrons and number of electrons needed to respect octate rule (for each element)
    valence = element(i).nvalence()
    remaining = 8 - valence

arrangment = {} #bonds dictionary with each included element haaving a list of bonds and locations and strengths

count = 0 #keeps track of how many elements so far
for i in range(len(smiles_input)): #for every symbol in smiles input

    if smiles_input[i].isalpha() and smiles_input[i].isupper(): #if current string position is alphabetical and uppercase (signaling start of a periodic symbol)
        count += 1 #element found = increase count

        if smiles_input[i+1].isalpha() and smiles_input[i+1].islower(): #if next character is again alphabetical but lowercase (meaning a 2 char periodic symbol)
            arrangment[f"{smiles_input[i]+smiles_input[i+1]}_{count}"] = { #create element in dictionary and define its location in the smiles input and its symbol and leave space for its bonds
                "location": [i,i+1],
                "symbol": f"{smiles_input[i]+smiles_input[i+1]}",
                "bonds": []
            }
            #begin seraching for bonds in front
            for j in range(i+2,len(smiles_input)): #from next char after end of atomic symbol, to end of smiles input
                if smiles_input[j].isalpha() and smiles_input[j].isupper():
                    if j != i+2: #if theres a new element but its not directly after the current element
                        break #all bonds are already collected
            
                    else: #if there is a new element but it is directly after the current element
                        if smiles_input[j+1].isalpha() and smiles_input[i+1].islower(): #if next character is again alphabetical but lowercase (meaning a 2 char periodic symbol)
                            arrangment[f"{smiles_input[i]+smiles_input[i+1]}_{count}"]["bonds"].append({ #append to the current elements bonds
                                "location": [j,j+1], #bounds of element symbol
                                "symbol": f"{smiles_input[j]+smiles_input[j+1]}",
                                "type": "-" #single bonds are infered
                            })
                        else: #if 1 char symbol
                            arrangment[f"{smiles_input[i]+smiles_input[i+1]}_{count}"]["bonds"].append({ #append to the current elements bonds
                                "location": [j],
                                "symbol": f"{smiles_input[j]}",
                                "type": "-" #single bonds are infered
                            })
                elif smiles_input[j] == "=": #if curent char is an open bracket (indicating a double bond)

                elif smiles_input[j] == "(": #if curent char is an open bracket (indicating start of a branch)
                    b = j+1 #start of branch
                    branch = "" #string of collected branch elements
                    while smiles_input[b] != ")": #loop through until you reach the end of the branch
                        branch += smiles_input[b]
                        b+=1
                    #after collecting branch terms
                    arrangment[f"{smiles_input[i]+smiles_input[i+1]}_{count}"]["bonds"].append({ #append to the current elements bonds
                        "location": [j,b], #bounds of branch
                        "symbol": f"{branch}", #collection of symbols inside branch (to be parsed later)
                        "type": "-" #single bonds are infered
                    })
                    j=b #skip past branch