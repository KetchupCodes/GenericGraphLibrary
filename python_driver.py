import graphlibrary

undG = graphlibrary.Graph_string_int(False)
dG = graphlibrary.Graph_string_int(True)

for i in ["A", "B", "C", "D", "E"]:
    undG.addNode(i)
    dG.addNode(i)

s = "ABCDE"
for i in range(len(s) - 1):
    undG.addEdge(s[i], s[i+1], i + 30)
    dG.addEdge(s[i], s[i+1], i + 15)

print("\n", undG.hasCycle())
print("\n", dG.hasCycle())

undG.addEdge("A", "C", 1)
dG.addEdge("A", "C", 1000)

print("\n", undG.hasCycle())
print("\n", dG.hasCycle())

dG.addEdge("C", "A", 2)

print("\n", dG.hasCycle())

undG.addEdge("A", "D", 200)
undG.addEdge("A", "E", 100)
undG.addEdge("B", "D", 10)

undGC = undG.nodeColoring()

print("\n\nNode Coloring: ")
for i in range(len(undGC)):
    print(undGC[i], end = " ")

dGC = undG.nodeColoring()

print("\n\nNode Coloring: ")
for i in range(len(dGC)):
    print(dGC[i], end = " ")

dGC = undG.edgeColoring()

print("\n\nEdge Coloring: ")
for i in range(len(dGC)):
    print(dGC[i], end = " ")

dGC = dG.nodeColoring()

print("\n\nNode Coloring: ")
for i in range(len(dGC)):
    print(dGC[i], end = " ")

dGC = dG.edgeColoring()

print("\n\nEdge Coloring: ")
for i in range(len(dGC)):
    print(dGC[i], end = " ")

print("\nconnected components: ")
cmps = undG.connectedComponents()
for j in cmps:
    for k in j:
        print(k, end=" ")
    print()

print("\nstrongly connected components: ")
cmps = dG.connectedComponents()
for j in cmps:
    for k in j:
        print(k, end=" ")
    print()

vec = undG.primMST()

print("\nprimMST\n")
for edg in vec:
    print(edg.getSource(), "-", edg.getDestination(), "\t", edg.getWeight())

vec = undG.kruskalMST()

print("\nkruskalMST\n")
for edg in vec:
    print(edg.getSource(), "-", edg.getDestination(), "\t", edg.getWeight())

print("\niterativeDFS")
for i in undG.iterativeDFS("A"):
    print(i, end=" ")
print()

vec1 = undG.uniformCostSearch("A", "E")
print("\nuniformCostSearch")
for i in vec1:
    print(i, end=' ')
print()

def heuristic(a, b):
    ai = 0
    bi = 0
    for x in a:
        ai += ord(x)
    for x in b:
        bi += ord(x)    
    return bi - ai

vec1 = dG.aStarSearch("A", "E", heuristic)
print("\naStarSearch")
for i in vec1:
    print(i, end=' ')
print()


vec2 = undG.katzCentrality(0.5, 0.25)
print("\nkatzCentrality")
for i in vec2:
    print(i, end = ' ')
print()

vec2 = dG.katzCentrality(0.5, 0.25)
print("\nkatzCentrality")
for i in vec2:
    print(i, end = ' ')
print()