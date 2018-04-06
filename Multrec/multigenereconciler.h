#ifndef MULTIGENERECONCILER_H
#define MULTIGENERECONCILER_H

#include <iostream>

#include <map>
#include "div/util.h"
#include "trees/newicklex.h"
#include "trees/node.h"
#include "trees/genespeciestreeutil.h"
#include "trees/treeiterator.h"

using namespace std;


/**
 * @brief The MultiGeneReconcilerInfo class is a basic structure to hold
 * various variables related to a partial mapping.  It is mainly used to pass all these
 * variables around into a single structure.
 */
class MultiGeneReconcilerInfo
{
public:
    unordered_map<Node*, Node*> partialMapping;
    int nbLosses;
    int dupHeightSum;
    bool isBad;

    MultiGeneReconcilerInfo()
    {
        isBad = false;
        nbLosses = 0;
        dupHeightSum = 0;
    }

    double GetCost(double dupcost, double losscost)
    {
        return dupcost * (double)dupHeightSum + losscost * (double)nbLosses;
    }

};



class MultiGeneReconciler
{
public:

    /**
     * @brief MultiGeneReconciler
     * @param geneTrees The set of gene trees contained in the forest.
     * @param speciesTree The species tree.
     * @param geneSpeciesMapping A mapping from the leaves of the gene trees to the leaves of the species tree.
     * @param dupcost The cost for one level of duplication.
     * @param losscost The cost for each loss.
     * @param maxDupHeight The maximum allowable duplication height.
     */
    MultiGeneReconciler(vector<Node*> &geneTrees, Node* speciesTree, unordered_map<Node*, Node*> &geneSpeciesMapping, double dupcost, double losscost, int maxDupHeight);

    /**
     * @brief Reconcile
     * Performs the reconciliation.  The return value contains the mapping, the sum of duplication heights and number of losses.
     * If isBad is true in the returned info, then it means there exists no solution.
     * @return
     */
    MultiGeneReconcilerInfo Reconcile();

    /**
     * @brief GetMappingCost
     * @param fullMapping A mapping of each node of each gene tree to the species tree.  We assume this mapping is valid without checking.
     * @return The total segmental dup + loss cost.
     */
    double GetMappingCost(unordered_map<Node*, Node*> &fullMapping);


    /**
     * @brief IsDuplication Returns true iff g is a duplication under partialMapping
     * @param g Internal node from some gene tree
     * @param partialMapping The current partial mapping.  Can actually be compelte.
     * @return true or false
     */
    bool IsDuplication(Node* g, unordered_map<Node*, Node*> &partialMapping);

private:

    vector<Node*> geneTrees;
    Node* speciesTree;
    unordered_map<Node*, Node*> geneSpeciesMapping;
    unordered_map<Node*, Node*> lcaMapping;
    double dupcost;
    double losscost;
    int maxDupHeight;

    //Main recursive function for the computation of a mapping.  Takes the partial mapping in info and tries to map additional
    //nodes.  The given mapping must be clean.  Returns a MultiGeneReconcilerInfo containing a complete mapping if one can be
    //found.  If not, the returned info object will have isBad = true.  The variable duplicationHeights has one entry for each species.
    MultiGeneReconcilerInfo ReconcileRecursive(MultiGeneReconcilerInfo &info, unordered_map<Node*, int> &duplicationHeights);

    //holds the current best solution, so that we can do some branch-and-bound early stop if we know we acnnot beat this in a recursion
    MultiGeneReconcilerInfo currentBestInfo;

    //key1 = species 1, key2 = species2, int = dist.  Not sure how well this performs, but let's try
    unordered_map< Node*, unordered_map<Node*, int> > speciesTreeDistances;

    //fills up the lcaMapping variables
    void ComputeLCAMapping();

    //true iff g is a key in partialMapping
    bool IsMapped(Node* g, unordered_map<Node*, Node*> &partialMapping);

    //returns the lsit of unmapped nodes whose two children are mapped
    vector<Node*> GetMinimalUnmappedNodes(unordered_map<Node*, Node*> &partialMapping);


    //returns the lowest node of the species tree on which g can be mapped, ie the lca of the mappings of the 2 children of g
    Node* GetLowestPossibleMapping(Node* g, unordered_map<Node*, Node*> &partialMapping);

    //false iff mapping g anywhere valid makes it a duplication
    bool IsRequiredDuplication(Node* g, unordered_map<Node*, Node*> &partialMapping);

    //true iff g is unmapped but its children are
    bool IsMinimalUnmapped(Node* g, unordered_map<Node*, Node*> &partialMapping);

    //true iff mapping g to its lowest possible place does not incerase duplication heights
    bool IsEasyDuplication(Node* g, unordered_map<Node*, Node*> &partialMapping, unordered_map<Node*, int> &duplicationHeights);

    //returns the duplication height at species of the subtree rooted at g
    int GetDuplicationHeightUnder(Node* g, Node* species, unordered_map<Node*, Node*> &partialMapping);



    //returns a node in the minimal nodes whose species is the lowest possible (multiple choices are possible, returns the first)
    Node* GetLowestMinimalNode(vector<Node*> &minimalNodes, unordered_map<Node*, Node*> &partialMapping);

    //returns the list of nodes on which minimalNode can be mapped to
    vector<Node*> GetPossibleSpeciesMapping(Node* minimalNode, unordered_map<Node*, Node*> &partialMapping);


    //Returns the number of edges between x and y in the species tree.  x and y must be comparable.
    int GetSpeciesTreeDistance(Node* x, Node* y);

    //Applies the cleaning phase on a partialMapping by mapping easy nodes until none are left.
    //This mapping can undergo modifications.  Only the minimal nodes and their ancestors can be modified.
    int CleanupPartialMapping(unordered_map<Node *, Node *> &partialMapping, unordered_map<Node*, int> &duplicationheights, vector<Node*> &minimalNodes);

};

#endif // MULTIGENERECONCILER_H
