#include "tree.h"

// // // // // // // // // // // // // // // // // // // // 
// 
// Class Tree: super simple tree class with both labeled nodes and edges.
// There is NOT a parent pointer in this simple tree.
// 
// version: 2.0  Apr 30, 2019
//

// construct a node without a name
Tree::Tree() : name("")
{
}


// construct a node with a name
Tree::Tree(std::string n) : name(n)
{
}

// construct a node with a name and a named edge to a child that is a tree
Tree::Tree(std::string n, std::string ln, Tree *c) : name(n)
{
    addChild(ln, c);
}


// construct a node with two named children and named edges
// useful for binary and operator trees
Tree::Tree(std::string n, std::string leftn, Tree *leftc, std::string rightn, Tree *rightc) : name(n)
{
    addChild(leftn, leftc, rightn, rightc);
}

// destructor (none at the moment)
Tree::~Tree()
{
}

// add a node supplied as a pointer as a child and name the edge
Tree *Tree::addChild(std::string ln, Tree *c)
{
    if (c==NULL) {
        if (ln=="") printf("ERROR(Tree::add): you cannot add a NULL pointer as a child.\n");
        else printf("ERROR(Tree::add): you cannot add a NULL pointer as a child under the edgename \"%s\".\n",
                    ln.c_str());
        exit(1);
    }
    children.push_back(c);
    edgeNames.push_back(ln);

    return this;
}


// add two node2 supplied as pointers as a children and name the edges
// NOTE: this really doesn't limit it to just two children but this
// is a convenience function if you are making a binary tree
Tree *Tree::addChild(std::string leftn, Tree *leftc, std::string rightn, Tree *rightc)
{
    if (leftc==NULL) {
        if (leftn=="") printf("ERROR(Tree::add): you cannot add a NULL pointer as a left child.\n");
        else printf("ERROR(Tree::add): you cannot add a NULL pointer as a right child under the edgename \"%s\".\n", leftn.c_str());
        exit(1);
    }
    if (rightc==NULL) {
        if (rightn=="") printf("ERROR(Tree::add): you cannot add a NULL pointer as a right child.\n");
        else printf("ERROR(Tree::add): you cannot add a NULL pointer as a right child under the edgename \"%s\".\n", rightn.c_str());
        exit(1);
    }
    children.push_back(leftc);
    edgeNames.push_back(leftn);
    children.push_back(rightc);
    edgeNames.push_back(rightn);

    return this;
}

// print out just the contents of a node and spot check its contents
void Tree::printNode()
{
    printf("Node: %s\n", name.c_str());
    if (children.size() != edgeNames.size()) {
        printf("WARNING: number of children of node is %d but the number of edgeNames (empty or not) is %d\n",
               int(children.size()),
               int(edgeNames.size()));
    }
    for (int i=0; i<children.size(); i++) {
        printf("%2d %s %s\n",
               i,
               (children[i] == NULL ? "NULL" : "Child"),
               edgeNames[i].c_str());
    }
}

// the print helper function.   Not public
// if edge is true then the edge names are printed one for each child
// if edge is false then the edge names are ignored
void Tree::printAux(bool edge, int level)
{
    if (children.size()==0) {
        for (int j=0; j<level; j++) printf("   ");
        printf("%s\n", (name=="" ? "NONAME" : name.c_str()));
    }
    else {
        if (! edge) {
            for (int j=0; j<level; j++) printf("   ");
            printf("%s\n", (name=="" ? "NONAME" : name.c_str()));
        }

        for (unsigned int i=0; i<children.size(); i++) {
            for (int j=0; j<level; j++) printf("   ");
            if (edge) printf("%s=%s\n",
                   (name=="" ? "NONAME" : name.c_str()),
                   (edgeNames[i]=="" ? "NONAME" : edgeNames[i].c_str()));
            if (children[i]) children[i]->printAux(edge, level+1);
        }
    }
}


// find child whose edge name is given.   Return NULL if not found
Tree *Tree::getChild(const std::string edgeName)
{
    for (int i=0; i<edgeNames.size(); i++) {
        if (edgeName == edgeNames[i]) return children[i];
    }

    return NULL;
}


/*
// // // // // // // // // // // // // // // // 
//
// main for testing
//
int main()
{
    Tree n((char *)"DOG");
    Tree *m;

    n.print();
    n.printWithEdges();
    n.printNode();

    printf("NC: %d\n", n.numChildren());
    n.addChild((char *)"woof", new Tree());
    n.addChild((char *)"bark", new Tree((char *)"CAT"));
    printf("NC: %d\n", n.numChildren());
    n.print();
    printf("\n");
    n.printWithEdges();
    n.printNode();
    m = new Tree((char *)"BAT", (char *)"squeak", new Tree((char *)"GUANO"));
    m->print();
    printf("\n");

    Tree m2((char *)"Ochotona",
            (char *)"Y", new Tree((char *)"Pica"),
            (char *)"N", new Tree((char *)"Leprodae", (char *)"Y", new Tree((char *)"rabbit"),
                                (char *)"N", new Tree((char *)"hare")));
    m2.print();
    printf("\n");
    printf("\n");
    m2.printWithEdges();
    m2.printNode();
    
    m = new Tree((char *)"A", (char *)"a1", new Tree((char *)"1"),
                 (char *)"a2", new Tree((char *)"2"));
    m->print();
    printf("\n");
    m->printWithEdges();

    return 0;
}

*/
