import numpy as np

def opentext(filename):
    lijst = [filename, '.txt']
    filename = ''.join(lijst)
    document = open(filename, 'r')
    text = document.read()
    document.close()
    return text

def cut_text(text, cutsym):
    listperline = text.split(cutsym)
    arrayperline = np.array(listperline)
    return arrayperline
def plotxyz(array):
    return
text1 = opentext('nodes')
nodeline = cut_text(text1, '\n')
nodeline0 = nodeline[0]

allnodes = []
for i in range(len(nodeline)):
    nodeline[i] = ''.join(nodeline[i].split())
    nodelist = cut_text(nodeline[i], ',')
    for i in range(len(nodelist)):
        nodelist[i] = float(nodelist[i])
    allnodes.append(nodelist)

text1 = opentext('elements')
nodeline = cut_text(text1, '\n')
nodeline0 = nodeline[0]
allelements = []
for i in range(len(nodeline)):
    nodeline[i] = ''.join(nodeline[i].split())
    nodelist = cut_text(nodeline[i], ',')
    for i in range(len(nodelist)):
        nodelist[i] = float(nodelist[i])
    allelements.append(nodelist)

map(float, allelements)

#print('allnodes', allnodes)
#print('allelements', allelements)

# starting to couple every element with 4 points.
elementarray = []
print(len(allelements))
for i in range(len(allelements)):
    print('elementnumber', i)

    # first node: find in the elementlist the right index and find its corresponding nodes.
    print('1', allnodes[int(allelements[i][1]) - 1])
    # second node
    print('2', allnodes[int(allelements[i][2]) - 1])
    # third node
    print('3', allnodes[int(allelements[i][3]) - 1])
    # fourth node
    print('4', allnodes[int(allelements[i][4]) - 1])

    # find values for element
    elementlist = [allnodes[int(allelements[i][1]) - 1], allnodes[int(allelements[i][2]) - 1],
                   allnodes[int(allelements[i][3]) - 1], allnodes[int(allelements[i][4]) - 1]]
    print(elementlist)
    elementarray.append(elementlist)
# print (elementarray)