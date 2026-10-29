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
            if i + 1 < len(smiles_input) and  smiles_input[i+1].isalpha() and smiles_input[i+1].islower(): #if next character is again alphabetical but lowercase (meaning a 2 char periodic symbol) (i + 1 < len(smiles_input) prevents checkign for 2 char if it is the last value in input)
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

        if i + 1 < len(smiles_input) and  smiles_input[i+1].isalpha() and smiles_input[i+1].islower(): #if next character is again alphabetical but lowercase (meaning a 2 char periodic symbol) (i + 1 < len(smiles_input) prevents checkign for 2 char if it is the last value in input)
            symbol = smiles_input[i]+smiles_input[i+1]
        else: #1 char element
            symbol = smiles_input[i]
        arrangment[f"{symbol}_{count}"] = { #create element in dictionary and define its location in the smiles input and its symbol and leave space for its bonds
            "location": list(range(i, i + len(symbol))), #lists positions it occupies
            "symbol": f"{symbol}",
            "bonds": []
        }

        #begin seraching for bonds in front
        j = i + len(symbol)
        branch_skipped = False #checks if a branch has been skipped (to ensure the bond after is still included despite not being adjacent)
        ring_skipped = False #checks if a ring has been skipped (to ensure the bond after is still included despite not being adjacent)
        while j < len(smiles_input): #from next char after end of atomic symbol, to end of smiles input

            if smiles_input[j].isalpha() and smiles_input[j].isupper():
                if j != i+len(symbol) and not branch_skipped and not ring_skipped: #if theres a new element but its not directly after the current element (doesnt apply if you skipped over a branch or ring)
                    branch_skipped = False #prevents infinate bonding after branch
                    ring_skipped = False #prevents infinate bonding after ring
                    break #all bonds are already collected
        
                else: #if there is a new element but it is directly after the current element
                    if j + 1 < len(smiles_input) and smiles_input[j+1].isalpha() and smiles_input[j+1].islower(): #if next character is again alphabetical but lowercase (meaning a 2 char periodic symbol) (j + 1 < len(smiles_input) prevents checking for 2 char if it is the last value in input)
                        arrangment[f"{symbol}_{count}"]["bonds"].append({ #append to the current elements bonds
                            "location": [j,j+1], #bounds of element symbol
                            "symbol": f"{smiles_input[j]+smiles_input[j+1]}",
                            "type": "-" #single bonds are infered
                        })
                    else: #if 1 char symbol
                        arrangment[f"{symbol}_{count}"]["bonds"].append({ #append to the current elements bonds
                            "location": [j],
                            "symbol": f"{smiles_input[j]}",
                            "type": "-" #single bonds are infered
                        })
                    break #all bonds are collected

            elif smiles_input[j] == "=" or smiles_input[j] == "#": #if curent char is an equals or a hashtag (indicating a double bond or triple bond)
                if smiles_input[j+1].isalpha() and smiles_input[j+1].isupper():
                    if j + 2 < len(smiles_input) and smiles_input[j+2].isalpha() and smiles_input[j+2].islower(): #(j + 2 < len(smiles_input) prevents checking for 2 char if it is the last value in input)
                        arrangment[f"{symbol}_{count}"]["bonds"].append({ #append to the current elements bonds
                            "location": [j+1,j+2],
                            "symbol": f"{smiles_input[j+1]+smiles_input[j+2]}", #2 char symbol in front of bond
                            "type": f"{smiles_input[j]}" #representing the double or triple bond
                        })
                    else: #1 char symbol after bond
                        arrangment[f"{symbol}_{count}"]["bonds"].append({ #append to the current elements bonds
                            "location": [j+1],
                            "symbol": f"{smiles_input[j+1]}", #1 char symbol in front of bond
                            "type": f"{smiles_input[j]}" #representing the double or triple bond
                        })
                elif smiles_input[j+1] == "(": #if double or triple bonded to branch
                    b = j+2 #start of branch
                    closes = 1 #number of ) to find before stopping (avoids stopping at the end of a nested branch instead of main branch)
                    branch = "" #string of collected branch elements
                    while True: #loop through until you reach the end of the branch
                        if smiles_input[b] == "(":
                            closes += 1 #an additionally close bracket is required to complete branch if a nested branch is detected
                        elif smiles_input[b] == ")":
                            closes -= 1 #1 closer to completing branch
                            if closes == 0:
                                break #avoids including the end bracket if it is the true end of the main branch
                        branch += smiles_input[b]
                        b+=1
                    #after collecting branch terms
                    arrangment[f"{symbol}_{count}"]["bonds"].append({ #append to the current elements bonds
                        "location": [j,b], #bounds of branch
                        "symbol": f"{branch}", #collection of symbols inside branch (to be parsed later)
                        "type": f"{smiles_input[j]}" #representing the double or triple bond
                    })
                    j=b+1 #skip past branch
                    branch_skipped = True
                    continue #return to top

                elif smiles_input[j+1].isnumeric(): #e.g. 1,2,3,4 ect (indicating start of a ring)
                    ring_start = smiles_input[j+1] #record number of ring start to find the end of the ring later
                    r = j+2 #start of ring
                    ring = "" #string of collected ring elements
                    while smiles_input[r] != ring_start: #loop through until you reach the end of the ring
                        ring += smiles_input[r]
                        r+=1
                    #after collecting branch terms
                    arrangment[f"{symbol}_{count}"]["bonds"].append({ #append to the current elements bonds
                        "location": [j,r], #bounds of ring
                        "symbol": f"{ring}", #collection of symbols inside ring (to be parsed later)
                        "type": f"{smiles_input[j]}" #representing the double or triple bond
                    })
                    j=r+1 #skip past ring
                    ring_skipped = True
                    continue #return to top

            elif smiles_input[j] == "(": #if curent char is an open bracket (indicating start of a branch)
                b = j+1 #start of branch
                closes = 1 #number of ) to find before stopping (avoids stopping at the end of a nested branch instead of main branch)
                branch = "" #string of collected branch elements
                while True: #loop through until you reach the end of the branch
                    if smiles_input[b] == "(":
                        closes += 1 #an additionally close bracket is required to complete branch if a nested branch is detected
                    elif smiles_input[b] == ")":
                        closes -= 1 #1 closer to completing branch
                        if closes == 0:
                            break #avoids including the end bracket if it is the true end of the main branch
                    branch += smiles_input[b]
                    b+=1
                #after collecting branch terms
                arrangment[f"{symbol}_{count}"]["bonds"].append({ #append to the current elements bonds
                    "location": [j,b], #bounds of branch
                    "symbol": f"{branch}", #collection of symbols inside branch (to be parsed later)
                    "type": "-" #single bonds are infered
                })
                j=b+1 #skip past branch
                branch_skipped = True
                continue #return to top

            elif smiles_input[j].isnumeric(): #e.g. 1,2,3,4 ect (indicating start of a ring)
                ring_start = smiles_input[j] #record number of ring start to find the end of the ring later
                r = j+1 #start of ring
                ring = "" #string of collected ring elements
                while smiles_input[r] != ring_start: #loop through until you reach the end of the ring
                    ring += smiles_input[r]
                    r+=1
                #after collecting branch terms
                arrangment[f"{symbol}_{count}"]["bonds"].append({ #append to the current elements bonds
                    "location": [j,r], #bounds of ring
                    "symbol": f"{ring}", #collection of symbols inside ring (to be parsed later)
                    "type": "-" #single bonds are infered
                })
                j=r+1 #skip past ring
                ring_skipped = True
                continue #return to top
            j += 1 #increment while loop

                
        #begin seraching for bonds behind (despite having branch logic branches are are not recorded as bonds but just skipped over to find the actuall bond on th other side) - (as branches can only be bonded with on the left side)
        j = i -1 #doesnt need to adjust for symbol length
        while j > -1: #from the char before the atomic symbol, to start of smiles input

            if smiles_input[j].isalpha(): #doesnt check if its upper becuse it still indicates new element 1 char or 2 char
                if j != i-1 and not branch_skipped and not ring_skipped: #if theres a new element but its not directly before the current element (doesnt apply if you skipped over a branch or ring)
                    branch_skipped = False #prevents infinate bonding after branch
                    ring_skipped = False #prevents infinate bonding after ring
                    break #all bonds are already collected
        
                else: #if there is a new element but it is directly before the current element
                    if smiles_input[j].isalpha() and smiles_input[j].islower(): #if previous character is again alphabetical but lowercase (meaning a 2 char periodic symbol)
                        arrangment[f"{symbol}_{count}"]["bonds"].append({ #append to the current elements bonds
                            "location": [j-1,j], #bounds of element symbol
                            "symbol": f"{smiles_input[j-1]+smiles_input[j]}",
                            "type": "-" #single bonds are infered
                        })
                    else: #if 1 char symbol
                        arrangment[f"{symbol}_{count}"]["bonds"].append({ #append to the current elements bonds
                            "location": [j],
                            "symbol": f"{smiles_input[j]}",
                            "type": "-" #single bonds are infered
                        })
                    break #all bonds are collected

            elif smiles_input[j] == "=" or smiles_input[j] == "#": #if curent char is an equals or a hashtag (indicating a double bond or triple bond)
                if smiles_input[j-1].isalpha() and smiles_input[j-1].isupper():
                    if smiles_input[j-2].isalpha() and smiles_input[j-2].islower():
                        arrangment[f"{symbol}_{count}"]["bonds"].append({ #append to the current elements bonds
                            "location": [j-2,j-1],
                            "symbol": f"{smiles_input[j-2]+smiles_input[j-1]}", #2 char symbol behind of bond
                            "type": f"{smiles_input[j]}" #representing the double or triple bond
                        })
                    else: #1 char symbol after bond
                        arrangment[f"{symbol}_{count}"]["bonds"].append({ #append to the current elements bonds
                            "location": [j-1],
                            "symbol": f"{smiles_input[j-1]}", #1 char symbol behind of bond
                            "type": f"{smiles_input[j]}" #representing the double or triple bond
                        })
                elif smiles_input[j-1] == ")": #if double or triple bonded to branch
                    b = j-2 #end of branch
                    closes = 1 #number of ) to find before stopping (avoids stopping at the start of a nested branch instead of main branch)
                    while True: #loop through until you reach the start of the branch
                        if smiles_input[b] == ")":
                            closes += 1 #an additionally open bracket is required to complete branch if a nested branch is detected
                        elif smiles_input[b] == "(":
                            closes -= 1 #1 closer to completing branch
                            if closes == 0:
                                break #avoids including the start bracket if it is the true start of the main branch
                        b-=1
                    j=b-1 #skip past branch
                    branch_skipped = True
                    continue #return to top

                elif smiles_input[j-1].isnumeric(): #e.g. 1,2,3,4 ect (indicating end of a ring)
                    ring_start = smiles_input[j-1] #record number of ring end to find the start of the ring later
                    r = j-2 #end of ring
                    ring = "" #string of collected ring elements
                    while smiles_input[r] != ring_start: #loop through until you reach the start of the ring
                        ring += smiles_input[r]
                        r-=1
                    #after collecting branch terms
                    arrangment[f"{symbol}_{count}"]["bonds"].append({ #append to the current elements bonds
                        "location": [r,j], #bounds of ring
                        "symbol": f"{ring}", #collection of symbols inside ring (to be parsed later)
                        "type": f"{smiles_input[j]}" #representing the double or triple bond
                    })
                    j=r-1 #skip past ring
                    ring_skipped = True
                    continue #return to top

            elif smiles_input[j] == ")": #if curent char is an open bracket (indicating end of a branch)
                b = j-1 #end of branch
                closes = 1 #number of ) to find before stopping (avoids stopping at the start of a nested branch instead of main branch)
                while True: #loop through until you reach the start of the branch
                    if smiles_input[b] == ")":
                        closes += 1 #an additionally open bracket is required to complete branch if a nested branch is detected
                    elif smiles_input[b] == "(":
                        closes -= 1 #1 closer to completing branch
                        if closes == 0:
                            break #avoids including the start bracket if it is the true start of the main branch
                    b-=1
                j=b-1 #skip past branch
                branch_skipped = True
                continue #return to top

            elif smiles_input[j].isnumeric(): #e.g. 1,2,3,4 ect (indicating end of a ring)
                ring_start = smiles_input[j] #record number of ring end to find the start of the ring later
                r = j-1 #end of ring
                ring = "" #string of collected ring elements
                while smiles_input[r] != ring_start: #loop through until you reach the startof the ring
                    ring += smiles_input[r]
                    r-=1
                #after collecting branch terms
                arrangment[f"{symbol}_{count}"]["bonds"].append({ #append to the current elements bonds
                    "location": [r, j], #bounds of ring
                    "symbol": f"{ring}", #collection of symbols inside ring (to be parsed later)
                    "type": "-" #single bonds are infered
                })
                j=r-1 #skip past ring
                ring_skipped = True
                continue #return to top
            j -= 1 #increment while loop

print(arrangment)