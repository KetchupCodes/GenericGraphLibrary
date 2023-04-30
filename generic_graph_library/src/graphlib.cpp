#include "graphlib.hpp"

using namespace std;

template <typename T, typename W>
Edge<T, W>::Edge(const T& source, const T& destination, const W& weight) {
    this->source = source;
    this->destination = destination;
    this->weight = weight;
}

template <typename T, typename W>
T Edge<T, W>::getSource() const {
    return this->source;
}
template <typename T, typename W>
T Edge<T, W>::getDestination() const {
    return this->destination;
}

template <typename T, typename W>
W Edge<T, W>::getWeight() const {
    return this->weight;
}

template <typename T, typename W>
bool Edge<T, W>::operator<(Edge<T, W>& edge2) {
    return this->weight < edge2.weight;
}

template <typename T, typename W>
Graph<T, W>::Graph(bool directed) {
    graph_size = 0;
    this->directed = directed;
}

template <typename T, typename W>
void Graph<T, W>::addNode(const T& data) {
    node_index[data] = graph_size;
    nodes.push_back(data);
    nodes_set.insert(data);

    vector<pair<int, W> > new_row;
    adj_list.push_back(new_row);
    graph_size++;
}

template <typename T, typename W>
void Graph<T, W>::addEdge(const T& source, const T& destination, const W& weight) {
    if (nodes_set.find(source) == nodes_set.end() || nodes_set.find(destination) == nodes_set.end()) {
        throw invalid_argument("Source/Destination doesn't exist");
    }
    
    Edge<T, W> new_edge(source, destination, weight);
    edges.push_back(new_edge);

    int source_index = node_index.at(source), destination_index = node_index.at(destination);
    adj_list[source_index].push_back({destination_index, weight});
    if (!directed) {
        adj_list[destination_index].push_back({source_index, weight});
    }
}

static bool findCycle(vector<vector<int>> &adj_list, vector<bool> &in_stack, vector<bool> &visited, int cur) {
    if (in_stack[cur]) {
        return true;
    } else if (visited[cur]) {
        return false;
    } 
    visited[cur] = in_stack[cur] = true;
    bool ans = false;
    for (auto next: adj_list[cur]) {
        ans |= findCycle(adj_list, in_stack, visited, next);
    }

    in_stack[cur] = false;
    return ans;
}

static bool findCycle(vector<vector<int>> &adj_list, vector<bool> &in_stack, vector<bool> &visited, int cur, int prev) {
    if (in_stack[cur]) {
        return true;
    } else if (visited[cur]) {
        return false;
    } 
    visited[cur] = in_stack[cur] = true;
    bool ans = false;
    for (auto next: adj_list[cur]) {
        if (next != prev)
            ans |= findCycle(adj_list, in_stack, visited, next, cur);
    }

    in_stack[cur] = false;
    return ans;
}

template <typename T, typename W>
bool Graph<T, W>::hasCycle() const {
    vector<vector<int>> adj_list_(graph_size);
    vector<bool> visited(graph_size, false), in_stack(graph_size, false);

    for (int i = 0; i < graph_size; i++) {
        for (auto pair_: adj_list[i]) {
            adj_list_[i].push_back(pair_.first);
        }
    }
    
    for (int cur = 0;  cur < graph_size; cur++) {
        if (directed && !visited[cur] && findCycle(adj_list_, in_stack, visited, cur)) {
            return true;
        }
        if (!directed && !visited[cur] && findCycle(adj_list_, in_stack, visited, cur, -1)) {
            return true;
        }
    }

    return false;
}

template <typename T, typename W>
vector<int> Graph<T, W>::nodeColoring() const {
    vector<int> colors(graph_size, 0);
    vector<vector<int>> adjacent_nodes(graph_size);
    for (int i = 0; i < graph_size; i++) {
        for (auto j: adj_list[i]) {
            adjacent_nodes[i].push_back(j.first);
            if (directed) {
                adjacent_nodes[j.first].push_back(i);
            }
        }
    } 

    for (int i = 0; i < graph_size; i++) {
        vector<bool> adj_colors(graph_size, false);
        for (auto j: adjacent_nodes[i]) {
            adj_colors[colors[j]] = true;
        }
        while (adj_colors[colors[i]]) {
            colors[i]++;
        }
    }

    return colors;
}

template <typename T, typename W>
vector<int> Graph<T, W>::edgeColoring() const {
    vector<int> colors(edges.size(), 0);
    set<pair<T, int>> node_color;
    for (int i = 0; i < edges.size(); i++) {
        while (node_color.find({edges[i].getSource(), colors[i]}) != node_color.end() || node_color.find({edges[i].getDestination(), colors[i]}) != node_color.end()) {
            colors[i]++;
        }
        node_color.insert({edges[i].getSource(), colors[i]});
        node_color.insert({edges[i].getDestination(), colors[i]});
    }

    return colors;
}

template <typename T, typename W>
void Graph<T, W>::completeEdges() {
    for (int i = 0; i < graph_size; i++) {
        vector<bool> visited(graph_size, false);
        for (auto j: adj_list[i]) {
            visited[j.first] = true;
        }
        for (int j = 0; j < graph_size; j++) {
            if (!visited[j] && j != i) {
                addEdge(nodes[i], nodes[j]);
            }
        }
    }
}


template <typename T, typename W>
vector<vector<T> > Graph<T, W>::connectedComponents() const {
    vector<vector<T> > components;
    vector<bool> visited(graph_size, false); 
    vector<T> component;
    function<void(int)> dfs = [&](int cur) {
        if (visited[cur]) {
            return;
        }
        visited[cur] = true;
        component.push_back(nodes[cur]);
        for (auto j: adj_list[cur]) {
            dfs(j.first);
        }
    };
    for (int i = 0; i < graph_size; i++) {
        component.clear();
        dfs(i);
        if (component.size()) {
            components.push_back(component);
        }
    }
    return components;
}

template <typename T, typename W>
vector<double> Graph<T, W>::katzCentrality(double alpha, double beta) const {
    vector<vector<bool> > A(graph_size, vector<bool>(graph_size, false));
    for (int i = 0; i < graph_size; i++) {
        for (auto j: adj_list[i]) {
            A[i][j.first] = true;
        }
    }
    vector<double> x(graph_size, 1.0 / graph_size);
    double diff;
    do {
        vector<double> x_1(graph_size, 0.0);
        for (int i = 0; i < graph_size; i++) {
            for (int j = 0; j < graph_size; j++) {
                x_1[i] += alpha * x[j] * (A[i][j] ? 1.0: 0.0);
            }
            x_1[i] += beta;
        }

        double sum = accumulate(x_1.begin(), x_1.end(), 0.0);
        for (auto &i: x_1) {
            i /= sum;
        }
        diff = 0.0;
        for (int i = 0; i < graph_size; i++) {
            diff += (x[i] - x_1[i]) * (x[i] - x_1[i]);
        }
        diff = sqrt(diff);

        x = x_1;
        
    } while (diff > 1e-15);
    return x;
}

template <typename T, typename W>
vector<Edge<T, W> > Graph<T, W>::primMST() const {
    if (directed) {
        throw invalid_argument("MST inapplicable to directed graphs");
    }
    vector<Edge<T, W> > ans;
    vector<bool> visited(graph_size, false);
    struct cmp {
        bool operator() (Edge<T, W> &a, Edge<T, W> &b) {
            return a.getWeight() > b.getWeight();
        }
    };
    priority_queue<Edge<T, W>, vector<Edge<T, W> >, cmp> pq;                
    visited[0] = true;
    for (auto j: adj_list[0]) {
        pq.push(Edge<T, W>(nodes[0], nodes[j.first], j.second));
    }
    while (!pq.empty()) {
        Edge<T, W> cur = pq.top();
        pq.pop();
        int v = node_index.at(cur.getDestination());
        if (visited[v]) {
            continue;
        }
        visited[v] = 1;
        ans.push_back(cur);
        for (auto j: adj_list[v]) {
            pq.push(Edge<T, W>(nodes[v], nodes[j.first], j.second));
        }
    }
    if (ans.size() != graph_size - 1) {
        throw invalid_argument("graph disconnected");
    }
    return ans;
}

template <typename T, typename W>
vector<Edge<T, W> > Graph<T, W>::kruskalMST() const {
    if (directed) {
        throw invalid_argument("MST inapplicable to directed graphs");
    }

    vector<Edge<T, W> > ans;
    vector<Edge<T, W> > edges_1 = edges; // new vector to avoid messing with `edges` vector's chronological ordering
    sort(edges_1.begin(), edges_1.end());

    int n = graph_size;
    vector<int> p(n);
    iota(p.begin(), p.end(), 0);
    vector<int> setSize(n, 1);

    function<int(int)> find_p = [&](int x) {
        if (p[x] == x) {
            return x;
        }
        return p[x] = find_p(p[x]);
    };

    function<bool(int, int)> un = [&](int a, int b) {
        a = find_p(a); b = find_p(b);
        if (a == b) {
            return false;
        }
        if (setSize[a] < setSize[b]) {
            swap(a, b);
        }
        setSize[a] += setSize[b];
        p[b] = a;
        return true;
    };

    for (auto edge: edges_1) {
        int a = node_index.at(edge.getSource()) ;
        int b = node_index.at(edge.getDestination()) ;

        if (un(a, b)) {
            ans.push_back(edge);
        } else {
            // 
        }
    }

    if (ans.size() != n - 1) {
        throw invalid_argument("graph disconnected");
    }

    return ans;
}

template <typename T, typename W>
vector<T> Graph<T, W>::iterativeDFS(const T& start) const {

    vector<T> traversal = {start};
    int n = graph_size;
    vector<bool> visited(n, false);
    vector<int> chk(n, 0);
    stack<int> st;

    auto strt = node_index.find(start);
    if (strt == node_index.end()) {
        throw invalid_argument("iterativeDFS(): start node invalid");
    }
    st.push((*strt).second); visited[(*strt).second] = true;

    while (!st.empty()) {
        int cur = st.top();
        if (chk[cur] < adj_list[cur].size()) {
            int nxt = adj_list[cur][chk[cur]].first;
            if (!visited[nxt]) {
                visited[nxt] = true;
                traversal.push_back(nodes[nxt]);
                st.push(nxt);
            }
            chk[cur]++;
        } else {
            st.pop();
        }
    }

    return traversal;
}

template <typename T, typename W>
vector<T> Graph<T, W>::uniformCostSearch(const T& start, const T& goal) const {

    return aStarSearch(start, goal, [](T a, T b) { return 0.0; });

}

template <typename T, typename W>
vector<T> Graph<T, W>::aStarSearch(const T& start, const T& goal, function<double(T, T)> heuristic) const {
    auto st = node_index.find(start), gl = node_index.find(goal);
    if (st == node_index.end() || gl == node_index.end()) {
        throw invalid_argument("aStarSearch(): start or goal node invalid");
    }

    int n = graph_size;
    vector<bool> visited(n, false);
    vector<bool> in_front(n, false);
    vector<int> prev(n, -1);
    vector<double> min_dist(n);

    struct custom {
        double f_x;
        double g_x;
        int indx;
        custom(double arg1, double arg2, int arg3): f_x(arg1), g_x(arg2), indx(arg3) {};
    };
    struct custom_cmpr {
        bool operator() (custom &arg1, custom &arg2) {
            return arg1.f_x > arg2.f_x;
        }
    };

    priority_queue<custom, vector<custom>, custom_cmpr> pq;

    int start_index = (*st).second, goal_index = (*gl).second;
    pq.push(custom(heuristic(start, goal) + W(), (double)W(), start_index)); in_front[start_index] = true;

    while (!visited[goal_index] && !pq.empty()) {
        int cur = pq.top().indx;
        double g_cur = pq.top().g_x;
        pq.pop();
        if (visited[cur]) continue;
        visited[cur] = true;

        for (auto j: adj_list[cur]) {
            int nxt = j.first;
            if (visited[nxt]) continue;

            double g_nxt = g_cur + j.second;
            double f_nxt = g_nxt + heuristic(nodes[nxt], goal);
            if (!in_front[nxt]) {
                in_front[nxt] = true;
                min_dist[nxt] = g_nxt;
                prev[nxt] = cur;
            }

            pq.push(custom(f_nxt, g_nxt, nxt));

            if (g_nxt < min_dist[nxt]) {
                min_dist[nxt] = g_nxt;
                prev[nxt] = cur;
            }
        }
    }
    
    if (!visited[goal_index]) {
        throw invalid_argument("aStarSearch(): goal unreachable from start");
    }

    vector<T> path;
    while (goal_index != -1) {
        path.push_back(nodes[goal_index]);
        goal_index = prev[goal_index];
    }
    reverse(path.begin(), path.end());

    return path;
}


template class Edge<int, float>;
template class Edge<int, int>;

template class Edge<char, float>;
template class Edge<char, int>;

template class Edge<float, float>;
template class Edge<float, int>;

template class Edge<std::string, float>;
template class Edge<std::string, int>;

template class Graph<int, float>;
template class Graph<int, int>;

template class Graph<char, float>;
template class Graph<char, int>;

template class Graph<float, float>;
template class Graph<float, int>;

template class Graph<std::string, float>;
template class Graph<std::string, int>;
