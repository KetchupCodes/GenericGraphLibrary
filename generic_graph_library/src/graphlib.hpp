#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <stack>
#include <queue>
#include <functional>
#include <numeric>
#include <math.h>
#include <stdexcept>

template <typename T, typename W> 
class Edge {
private:
    T source, destination;
    W weight;

public:
    Edge(const T& source, const T& destination, const W& weight);
    T getSource() const;
    T getDestination() const;
    W getWeight() const;
    bool operator<(Edge<T, W> &edge2);
};

template <typename T, typename W>
class Graph {
private:
    bool directed;
    std::vector<T> nodes;
    std::vector<Edge<T, W> > edges;
    std::vector<std::vector<std::pair<int, W> > > adj_list;
    std::set<T> nodes_set;
    std::map<T, int> node_index;
    int graph_size;
    
public:
    Graph(bool directed);
    void addNode(const T& node_name);
    void addEdge(const T& source, const T& destination, const W& weight = W());
    bool hasCycle() const;

    // Algorithm implementations
    std::vector<int> nodeColoring() const;
    std::vector<int> edgeColoring() const;
    void completeEdges();
    std::vector<std::vector<T> > connectedComponents() const;
    std::vector<double> katzCentrality(double alpha, double beta) const;
    std::vector<Edge<T, W> > primMST() const;
    std::vector<Edge<T, W> > kruskalMST() const;
    std::vector<T> iterativeDFS(const T& start) const;
    std::vector<T> uniformCostSearch(const T& start, const T& goal) const;
    std::vector<T> aStarSearch(const T& start, const T& goal, std::function<double(T, T)> heuristic) const;
};