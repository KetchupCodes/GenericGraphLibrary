#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <iostream>
#include "graphlib.hpp"

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)

namespace py = pybind11;

template <typename T, typename W> 
void declare_edge(py::module &m, const std::string &typestr1, const std::string &typestr2) {
    std::string pyclass_name = "Edge_" + typestr1 + "_" + typestr2;
    py::class_<Edge<T, W> >(m, pyclass_name.c_str(), py::buffer_protocol(), py::dynamic_attr())
    .def(py::init<const T&, const T&, const W&>())
    .def("getSource", &Edge<T, W>::getSource)
    .def("getDestination", &Edge<T, W>::getDestination)
    .def("getWeight", &Edge<T, W>::getWeight)
    ;
}

template <typename T, typename W> 
void declare_graph(py::module &m, const std::string &typestr1, const std::string &typestr2) {

    declare_edge<T, W>(m, typestr1, typestr2);

    std::string pyclass_name = "Graph_" + typestr1 + "_" + typestr2;
    py::class_<Graph<T, W> >(m, pyclass_name.c_str(), py::buffer_protocol(), py::dynamic_attr())
    .def(py::init<bool>())
    .def("addNode", &Graph<T, W>::addNode)
    .def("addEdge", &Graph<T, W>::addEdge)
    .def("hasCycle", &Graph<T, W>::hasCycle)
    .def("nodeColoring", &Graph<T, W>::nodeColoring)
    .def("edgeColoring", &Graph<T, W>::edgeColoring)
    .def("completeEdges", &Graph<T, W>::completeEdges)
    .def("connectedComponents", &Graph<T, W>::connectedComponents)
    .def("katzCentrality", &Graph<T, W>::katzCentrality)
    .def("primMST", &Graph<T, W>::primMST)
    .def("kruskalMST", &Graph<T, W>::kruskalMST)
    .def("iterativeDFS", &Graph<T, W>::iterativeDFS)
    .def("uniformCostSearch", &Graph<T, W>::uniformCostSearch)
    .def("aStarSearch", &Graph<T, W>::aStarSearch)
    ;
}

PYBIND11_MODULE(graphlibrary, m) {
    m.doc() = R"pbdoc(
        Generic Graph Library

    )pbdoc";

    declare_graph<int, float>(m, "int", "float");
    declare_graph<int, int>(m, "int", "int");

    declare_graph<char, float>(m, "char", "float");
    declare_graph<char, int>(m, "char", "int");

    declare_graph<float, float>(m, "float", "float");
    declare_graph<float, int>(m, "float", "int");

    declare_graph<std::string, float>(m, "string", "float");
    declare_graph<std::string, int>(m, "string", "int");



#ifdef VERSION_INFO
    m.attr("__version__") = MACRO_STRINGIFY(VERSION_INFO);
#else
    m.attr("__version__") = "dev";
#endif
}
