from mendeleev import element #for extracting infomation regarding elements e.g. valence
import sys #for ending program
import urllib.request #for url retrieval from pubchem api
import urllib.parse
import json
from OpenGL.GL import *
from OpenGL.GLU import *
from OpenGL.GLUT import *
import math

vertex_count_lat = 8 #vertical lines
vertex_count_lng = 8 #horozontal lines
radius = 0.5 #row of sphere

#first previous mouse values
previous_x = 0
previous_y = 0

#initial polar coordiantes for camera
camera_theta = 0 
camera_phi =  math.pi / 2 #level spheres equator
camera_r = 15 #semi zoomed out

#initial cartesian coordiantes for camera (perfectly line up wuith init polar coordinates)
camera_x = 3 
camera_y = 0
camera_z = 0

click = False

sensitivity = 0.003 #how sensitive the camera rotation is to mouse movement

def update(Value): #update time and calculate mouse accelaration (value is for the manditory data argument from glut timer func)
    glutPostRedisplay() #redraw the display
    glEnable(GL_DEPTH_TEST) #allows closer surfaces to hide further/obstructed ones
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
    camera()

    glutPassiveMotionFunc(mouse_movement) #calls mouse_movement function with x and y values included
    glutMouseFunc(mouse_click) #detects button presses, same idea as above

    glutTimerFunc(16, update, 0) #schedule update in 16 millisceonds (about 60fps -  1/60 = 0.0166..)

def mouse_movement(x,y): #handles mouse input and calculates mouse displacement used to rotate camera
    global previous_x, previous_y, camera_phi, camera_theta, camera_r #the last recorded mouse position used to calculate mouse displacement and the phi and theta and r values for the camera
    global camera_x, camera_y, camera_z #camera vlaues in cartesian co ordinates
    if click == True:
        x_displacement = x-previous_x #changes in x and y every frame
        y_displacement = y-previous_y
        print(x_displacement, y_displacement)

        camera_theta -= x_displacement * sensitivity #ccumulating mouse changes to camera angles multiplied by a sensitivity dampener
        camera_phi -= y_displacement * sensitivity

        camera_x = camera_r * math.cos(camera_theta) * math.sin(camera_phi) #uses same principles as for drawring the sphere for finding position of camera
        camera_y = camera_r * math.sin(camera_theta) * math.sin(camera_phi)
        camera_z = camera_r * math.cos(camera_phi)

    #update previous recorded mouse locations now camera calculations are done
    previous_x = x
    previous_y = y

def mouse_click(button, state, x, y): #only allows passive camera movement following left clicks
    global click #updates click condition
    if button == 0 and state == 1: #left mouse click button is clicked
        if click == False: #and if click is currently stll false
            click = True
        else:
            click = False

    print(button, state)

def camera():
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluPerspective(45, 1, 0.1, 100)

    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()

    gluLookAt(
        camera_x, camera_y, camera_z,    #cameras real 3d position
        0, 0, 0,  #what position the camera looks at
        0, 0, 1     #witch direction is up (z -defualt)
    )

def draw_sphere(offset): #defining a opengl display function

    cx,cy,cz = offset

    glColor3f(1, 0, 0)

    #solid sphere(triangles)
    for i in range (vertex_count_lng): #iterating through vertical lines
        phi = (math.pi) * (i/vertex_count_lng) #angle of vertical steps for any given vertex count by dividing vertex percentage multiplyer by radians in semi circle
        next_phi = math.pi * (i + 1) / vertex_count_lng #step angle for next vertex down
        for j in range (vertex_count_lat): #iterating through horozontal lines
            theta = 2*(math.pi) * (j/vertex_count_lat) #angle of horozontal steps for any given vertex count by dividing vertex percentage multiplyer by radians in circle
            next_theta = 2 * math.pi * (j + 1) / vertex_count_lat #step angle for next vertex across

            #(v0)----(v1)
            # |       |
            # |       |
            #(v2)----(v3)

            #vertices of the quad
            #open gl can only render lines points and triangle sos the quad neds to be broken down into 2 triangles

            #covert polar verticies to cartesian co ordiantes
            v0 = (cx + radius * math.cos(theta) * math.sin(phi), #finding opposite side (width wize length) using cos(on xy plane) while additionally multiplying by sin(phi) (proportion between hypot and opposite side on z plane - scales radius with polar angle) +offsets
                cy + radius * math.sin(theta) * math.sin(phi), #finding adjacent side (depth wize length) using sin(on xy plane) while additionally multiplying by sin(phi) (proportion between hypot and opposite side on z plane - scales radius with polar angle) +offsets
                cz + radius * math.cos(phi)) #finding adjacent side (height length) using cos(on z plane) +offsets
            
            v1 = (cx + radius * math.cos(next_theta) * math.sin(phi), #repeated steps for further angle increments to form all four verticies in grid pattern - to later construct into 2 triangles
                cy + radius * math.sin(next_theta) * math.sin(phi),
                cz + radius * math.cos(phi))
            
            v2 = (cx + radius * math.cos(theta) * math.sin(next_phi),
                cy + radius * math.sin(theta) * math.sin(next_phi),
                cz + radius * math.cos(next_phi))
            
            v3 = (cx + radius * math.cos(next_theta) * math.sin(next_phi),
                cy + radius * math.sin(next_theta) * math.sin(next_phi),
                cz + radius * math.cos(next_phi))
            
            glBegin(GL_TRIANGLES) #every 3 vertices stated past this declaration form 1 triangle

            # triangle 1
            glVertex3f(*v0) #unpacking x, y and z for each vertex
            glVertex3f(*v2)
            glVertex3f(*v1)

            # triangle 2
            glVertex3f(*v1)
            glVertex3f(*v2)
            glVertex3f(*v3)

            glEnd()

    #sphere wireframe(line loop)

    glColor3f(0, 0, 0)

    #latitude lines
    for i in range(vertex_count_lng): #uses same math for fidnding phi and theta increments and converting polar to cartesian co ordinates but instead has each vertice added to a wireframe loop (repeated below for longtitude)
        phi = math.pi * (i / vertex_count_lng)
        glBegin(GL_LINE_STRIP)
        for j in range(vertex_count_lat+1): #inner loops contain +1 to let them overlap (complete full latitudinal ring as gl line strip requires overlap to the start position)
            theta = 2 * math.pi * (j / vertex_count_lat)
            x = cx + radius * math.cos(theta) * math.sin(phi)
            y = cy + radius * math.sin(theta) * math.sin(phi)
            z = cz + radius * math.cos(phi)
            glVertex3f(x, y, z)
        glEnd()

    #longitude lines
    for j in range(vertex_count_lat):
        theta = 2 * math.pi * (j / vertex_count_lat)
        glBegin(GL_LINE_STRIP)
        for i in range(vertex_count_lng+1): #still contains +1 to reach final semi circle position
            phi = math.pi * (i / vertex_count_lng)
            x = cx + radius * math.cos(theta) * math.sin(phi)
            y = cy + radius * math.sin(theta) * math.sin(phi)
            z = cz + radius * math.cos(phi)
            glVertex3f(x, y, z)
        glEnd()

#<=======================input=======================>


smiles_input = input("Please input a molecule's name (all caps) or its code in SMILES format: -Note metals are not currently supported in this version and molecular names have to be entered in all capitals").replace(" ", "") # takes in input and removes all white space throughout

if smiles_input.isalpha() and smiles_input.isupper() and len(smiles_input) > 1: #checks if whole string is all caps and all aplphabetical and greater than 1 charachter (meaning molecular name isntead of smiles code)
    smiles_input =  urllib.parse.quote(smiles_input) #makes the smile input safe for url
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{smiles_input}/property/ConnectivitySMILES/JSON" #url to retrieve smiles code for molecule name 
    try:
        response = urllib.request.urlopen(url, timeout=10) #send api request to url(pubchem) -(wait maximum of 10s before raising timeout error)
        data = json.loads(response.read().decode("utf-8")) #read and decode api data
        print(data)
    except: #if error fetching smiles data
        print("Molecule name not recognised or found in database, please check your spelling or variation of the molecular name and try again, or enter manually using SMILES formatting")
        sys.exit()

    smiles_input = data["PropertyTable"]["Properties"][0]["ConnectivitySMILES"] #path to the smiles code from api result

#<=======================validation=======================>

print("retrieving moecular data...")

non_metal_series = [ #series that make up the non metals (used to exclude metals)
    "Noble gases",
    "Halogens",
    "Nonmetals", #(not currently including metaloids due to diffrent rules)
]

elements_included = [] #list of element symbols in order of mention in smiles input
valences = [] #list of element valences (in alignment order to "elements_included")
for i in range(len(smiles_input)):
    if smiles_input[i].isalpha() and smiles_input[i].isupper(): #if current string position is alphabetical and uppercase (signaling start of a periodic symbol)
        if i + 1 < len(smiles_input) and  smiles_input[i+1].isalpha() and smiles_input[i+1].islower(): #if next character is again alphabetical but lowercase (meaning a 2 char periodic symbol) (i + 1 < len(smiles_input) prevents checkign for 2 char if it is the last value in input)
            symbol = smiles_input[i] + smiles_input[i+1] #combine strings to form full periodic symbol
            elements_included.append(symbol)
        else: #if only a 1 char symbol
            symbol = smiles_input[i]
            elements_included.append(symbol) #only append the 1 char as it is the full symbol

print(elements_included)
valence_data = {} #store valence data during validation to minimize later queries on validated elements
for i in elements_included: #warn about metals before exit
    current_element = element(i) #used to minimize mendeleev queries - thus
    try: #catch error at very first instance of caling data about element (indicating no such element exists)
        if current_element.series not in non_metal_series: #check element is non metal
            print("Metals where detected in your SMILES input -please only use non metals (currently not including metaloids)")
            sys.exit()
        elif current_element.series == "Noble gases": #check element is not a noble gas
            print("Noble gases where detected in your SMILES input -these do not form bonds")
            sys.exit()
        valence_data[i] = {
            "valence": current_element.nvalence()
        }
    except:
        print(f"Please check you correctly spelled or entered your elements correctly, or that they actually exist: {i} is not recognised")
        sys.exit()

for i in elements_included: #number of valence electrons and number of electrons needed to respect octate rule (for each element)
    valence = valence_data[i]["valence"]
    remaining = 8 - valence
    valences.append(valence)

#<=======================bonding arrangment=======================>

print("decoding molecule...")

def add_bond(a, b, bond_type): #function to add a bond to both parties envolved in bond
    for x, y in ((a,b), (b,a)): #bond element a and element b from both atoms perspectives, e.g. carbon bonds with oxygen and oxygen bonds with carbon
        arrangment[x]["bonds"].append({ #append other atom y to atoms x's bonds
            "atom": y,
            "location": arrangment[y]["location"], #works as we only move forward so the atom we are bonding to has already been created
            "symbol": arrangment[y]["symbol"],
            "bond_type": bond_type
        })

arrangment = {} #bonds dictionary with each included element haaving a list of bonds and locations and strengths
stack = [] #list of elements before the start of a branch to bond with at the end of the branch, in order of discovary to allow for nesting
table = {} #dictionary for elements before the start of rings to bond with at the end of the ring, addressed by number
pending = "-" #unless otherwise stated single bonds are infered
last = None #no previous atom yet
count = 0 #atom number 0

i = 0
while i < len(smiles_input): #loop through the smiles input forwards
    char = smiles_input[i] #current character
    if char.isalpha() and char.isupper(): #if an uppercase alphabetical symbol (indicating start of an atomic symbol)
        if i+1 < len(smiles_input) and smiles_input[i+1].isalpha() and smiles_input[i+1].islower(): #if smiles input is long enough to contain another char and that char is lowecase alphabetical (meaning 2 char atomic symbol)
            symbol = char + smiles_input[i+1] #full 2char atomic symbol string
        else:
            symbol = char #1 char atomic symbol

        count += 1 #atom discovered
        atom = f"{symbol}_{count}" #name of the atom made up of its count and atomic symbol
        arrangment[atom] = {
            "symbol": symbol,
            "location": list(range(i,i+len(symbol))), #all position the atom covers/exists in (in the string)
            "valence": valences[count-1], #valences lsit is synced with order of atoms in smiles input
            "bonds": [], #empty until add bond function is called
            "lone_pairs": None, #pairs of electrons not involved in bonding that occupy an electron domain
            "azimuthal": None, #azimuthal angle from center of bonded-to element
            "polar": None #polar angle from center of bonded-to element
        }
        if not last == None: #if there is a previous element to bind with
            add_bond(atom, last, pending) #add bond between the two atoms
        last = atom #set the alst atom the the atom just proccesed
        pending = "-" #reset the bond type to single as single bonds can be infered so are default
        i += len(symbol) #move i forward by the length of the symbol
        continue #return to top to avoid incrementing i additonally

    elif char in "-=#": #current character is an explicit bond
        pending = char #set the current bond type to the char
    elif char == "(": #begining of a branch
        stack.append(last) #append the last verified atom to the stack (waiting for bonds outside of branch)
    elif char == ")": #end of a branch
        last = stack.pop() #remove the last atom from the stack as the branch it preceeded is now complete (additonally set the last element to the removded stack item, as it is able to form bonds outside the branch again)
    elif char.isnumeric() and char not in table: #if current character is a number (indicating ring start) and is not already recorded in the dictionary
        table[char] = {
            "start": last, #sets the rings start value to bond with the end later
            "start-end_bond": pending #bond before the ring number is the bond for the end ring atom
        }
        pending = "-" #reset pending to aovid start-end bond interprited as bond type with next ring atom
    elif char.isnumeric() and char in table: #already recorded ring start indicating ring closure
        start = table[char]["start"] #set start to that rings start value
        if table[char]["start-end_bond"] != "-":
            ring_bond = table[char]["start-end_bond"]
        else:
            ring_bond = pending #overwrite it with bond before ring closure (but if not mentioned it will be "-" by defualt so check if there was an explicit bond before the ring start)
        add_bond(start, last, ring_bond) #adds bond between start and end atoms of the ring
        table.pop(char) #removes dictionary item to allow that number to be reused later in the smiles input
        pending = "-" #reset pending to aovid start-end bond interprited as bond type with next atom

    i += 1 #increment i by 1 as all mentioned charachters where onyl 1 char long

#<=======================add implicit hydrogens + calculate lone pairs=======================>

for i in list(arrangment.keys()): #iterates through every key in dictionary via a copied list of the dictionary so that hwen appending to the original dictioanry the size of the iteration list doesnt change
    bond_types = [] #lsit of bond types per atom
    for j in arrangment[i]["bonds"]: #iterates through every bond in the atom
        bond_types.append(j["bond_type"]) #bond type for each bond


    value = 0
    for j in bond_types:
        bonds = "-=#"
        #total of all bond strengths (how many electrons it has gained from its covelant bonds)
        value += bonds.index(j) + 1 #retrieves current bond type from "bond_types" list then find the index of that bond in the string "bonds" and records its index (as the string is ordered by bond strength; index+1 = what level bond it is)

    if arrangment[i]["symbol"]  == "H": #explicit hydrogens are only edge case to octet rule as hydrogen and helium are the only non metals with 1 shell and helium is a noble gas so cannot be considered in bonding
        target = 2
    else:
        target = 8
    remaining_bonds = (target-arrangment[i]["valence"]-value) #how many electrons it still needs to gain from hydrogen bonds

    #//2 as lone "pairs" not lone electrons, as only lone pairs can occupy a electron domain
    lone_pairs = (arrangment[i]["valence"]-(target-arrangment[i]["valence"]))//2 #how many valence elctrons it has minus how many it needs to gain (as every electron you need to gain you need to give one away to share) - thus how many electrons are not involved in bonding
    arrangment[i]["lone_pairs"] = (lone_pairs) #set number of lone pairs

    j=0
    while j < remaining_bonds:
        symbol = "H"
        atom = f"{symbol}_{count}"
        arrangment[atom] = {
            "symbol": symbol,
            "location": "-", #no location in smiles input
            "valence": 1, 
            "bonds": [], #empty until add bond function is called
            "lone_pairs": 0, #hydrogen has 0 lone pairs always
            "azimuthal": None, #azimuthal angle from center of bonded-to element
            "polar": None #polar angle from center of bonded-to element
        }
        count += 1
        add_bond(i, atom, "-") #implicit hydrogens alwyas have single bonds
        j+=1

#<=======================bond angles=======================>

print("rendering...")

for i in arrangment:
    if arrangment[i]["symbol"] == "H": #hydrogens only form 1 electron domain(their only single bond) so bond at 180 degress regardless (this includes both explicit hydeogens and implicit hydrogens that dont have lone pair attributes to allow for accurate bond angle calculations)
        bond_angles = [[0, 0]] #only a single bond
    else:
        electron_domains = len(arrangment[i]["bonds"]) + arrangment[i]["lone_pairs"] #sum of number of bonds and lone pairs

        #decide the furthest possible locations for all points around the sphere (changes for number of points(electron domains)) (position 0 is azimuthal and position 1 is polar) -hardcoded for now just for simplicity
        #doesnt take into account lone pairs other than their existance as a electron domain (e.g. no electron bullying)
        if electron_domains == 2:
            bond_angles = [[0, 0], [0, 180]] 
        elif electron_domains == 3:
            bond_angles = [[0, 0], [0, 120], [0, 240]] 
        elif electron_domains == 4:
            bond_angles = [[0, 0], [109.5, 0], [109.5, 120], [109.5, 240]]
        elif electron_domains == 5:
            bond_angles = [[0, 0], [0, 180], [0, 90], [120, 90], [240, 90]]
        elif electron_domains == 6:
            bond_angles = [[0, 0], [0, 180], [0, 90], [90, 90], [180, 90], [270, 90]]


        position = 0 #which bond angle is chosen per bond
        for j in arrangment[i]["bonds"]:
            j["polar"] = bond_angles[position][0] #azimuthal angle from center
            j["azimuthal"] = bond_angles[position][1] #azimuthal angle from center
            position += 1 #move to a new bond angle

print(arrangment)

##<=======================3D render=======================>

print("molecule render complete!")

glutInit() #initializes GLUT
glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH) #condigures window behaviour eg colour format
glutInitWindowSize(800, 800) #sets initial window size
glutCreateWindow(b"Solid Sphere") #creates window

def display(): #all displaying has to be in a single function so glut display func can display them
    bond_length = 2 #length of bonds (starting and ending at the center of each sphere)
    offset_stack = [(0,0,0)] #stack of offsets to iterate through (used to spawn each atom) (starting with 0,0)
    atom_stack = [list(arrangment.keys())[0]] #stack of atom to iterate through (starting with the first atom (via lsiting keys in dictionary and taking the first one)) -used to retrieve information about the atom
    molecule_rendered = False
    a=0
    while not molecule_rendered:

        atom = atom_stack[a] #always pick the first atom in the stack
        draw_sphere(offset_stack[a]) #draw sphere at first offset in stack
        for i in arrangment[atom]["bonds"]:

            if i["atom"] in atom_stack: #checks if that bond has already been rendered
                continue
            polar = i["polar"] * (math.pi/180) #converts the bonds polar angle to radians
            azimuthal = i["azimuthal"] * (math.pi/180) #converts the bonds azimuthal angle to radians

            current_x, current_y, current_z = offset_stack[a] #offsets of current atom

            #coverting bond angles away into cartesian offsets (azimuth and polar are equivelant to theta and phi an dbond length is equivelent to row/radius)
            x = current_x + bond_length * math.cos(azimuthal) * math.sin(polar)
            y = current_y + bond_length * math.sin(azimuthal) * math.sin(polar)
            z = current_z + bond_length * math.cos(polar)

            offset_stack.append((x,y,z)) #offsets for new atom
            atom_stack.append(i["atom"]) #append the new atom bonded to current atom

        if a+1 == len(atom_stack): #if all atoms have been rendered (by checking if there is another atom to be rendered)
            molecule_rendered = True

        a+=1 #render the next atom

    glutSwapBuffers()


glutDisplayFunc(display) #calls function when window needs to be drawn

glutTimerFunc(16, update, 0) #starts update function for first time

glutMainLoop() #starts event loop

#To-do:
#we currently handle bond angles locally to their bonded atom but that doesnt take into account the rotaion of the atom we are bonding to itself 1st
#add cylinders in between atoms to represent bonds
#add ability to start new molecules to define mixtures using "." 2nd
#add more validation e.g. for missing brackets ot in future ring closes 3rd