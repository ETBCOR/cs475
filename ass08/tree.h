#ifndef TREEH
#define TREEH

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>

// // // // // // // // // // // // // // // // // // // // 
// 
// Class Tree: super simple tree class with both labeled nodes and edges.
// There is NOT a parent pointer in this simple tree.
// 
// version: 2.0  Apr 30, 2019
//
class Tree {
protected:
    std::string name;                      // name of the node
    std::vector<Tree *> children;          // pointers to children
    std::vector<std::string> edgeNames;    // names of the edges to each child ("" is no name)

// constructors
// some of these constructors allow the adding of edges at construction time for convenience
// otherwise use addChild method
public:
    Tree();                                               // empty tree
    Tree(const std::string n);                            // named root node
    Tree(const std::string n, const std::string ln, Tree *c);   // named root and edge name + one child
    Tree(const std::string n,                             // named root + two children (more can be added)
         const std::string leftn, Tree *leftc,
         const std::string rightn, Tree *rightc);
    ~Tree();

protected:
    void printAux(bool edge, int level);                  

public:
    bool isLeaf() { return children.size()==0; }    // is this a leaf (no children)
    int numChildren() { if(children.empty()) return 0; return children.size(); }   // how many children
    Tree *addChild(std::string ln, Tree *c);        // add a child with name and tree
    Tree *addChild(std::string leftn, Tree *leftc, std::string  rightn, Tree *rightc); // add children with names and trees
    void printNode();                               // dump contents of a node (useful for debugging)
    void printWithEdges() { printAux(true, 0); }         // print with edge names
    void print() { printAux(false, 0); }            // print tree without edge names
    std::string getName() { return name; }          // get the name of the tree
    Tree *getChild(const std::string edgeName);     // return tree pointed to by edge with the given name.  NULL if not there
//    void printChildren(){
//        for(Tree* child : children)
//            printf("child name: %s\n", child->getName().c_str());    
//    }
};

#endif
