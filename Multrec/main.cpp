#include <iostream>

#include <map>
#include "div/util.h"
#include "trees/newicklex.h"
#include "trees/node.h"
#include "trees/genespeciestreeutil.h"
#include "trees/treeiterator.h"
#include "multigenereconciler.h"

using namespace std;


int verbose = 0;


//labels the gene trees to prepare them for output.  Adds the species mapping to the out,
//plus _Spec or _Dup_nbx, where x is a dup id.  Also returns a map of dups per species,
//since we're computing it in this function anyway.  The value is a vector of int/Node pairs,
//where the int is the gene tree index and the node is a dup node in this tree.
//I agree that this function does more than just labelling with species mapping, but hey, life is tough.
map<Node*, vector< pair<int, Node*> > > LabelGeneTreesWithSpeciesMapping(vector<Node*> geneTrees, Node* speciesTree, MultiGeneReconciler &reconciler, MultiGeneReconcilerInfo &info, bool resetLabels = true)
{
    map<Node*, vector< pair<int, Node*> > > dups_per_species;
    int dup_counter = 1;
    for (int i = 0; i < geneTrees.size(); i++)
    {
        Node* genetree = geneTrees[i];

        TreeIterator* it = genetree->GetPostOrderIterator();
        while (Node* g = it->next())
        {
            if (!g->IsLeaf())
            {
                string lbl = g->GetLabel();
                if (lbl != "")
                    lbl += "_";

                if (resetLabels)
                    lbl = "";

                lbl += info.partialMapping[g]->GetLabel();

                if (reconciler.IsDuplication(g, info.partialMapping))
                {
                    lbl += "_Dup_nb" + Util::ToString(dup_counter);
                    dup_counter++;

                    vector< pair<int, Node*> > dups_for_s;
                    if (dups_per_species.find( info.partialMapping[g] ) != dups_per_species.end())
                        dups_for_s = dups_per_species[info.partialMapping[g]];
                    pair<int, Node*> p = make_pair(i + 1, g);
                    dups_for_s.push_back( p );
                    dups_per_species[info.partialMapping[g]] = dups_for_s;
                }
                else
                    lbl += "_Spec";
                g->SetLabel(lbl);
            }
        }
        genetree->CloseIterator(it);
    }

    return dups_per_species;
}

unordered_map<Node*, Node*> GetGeneSpeciesMapping(vector<Node*> geneTrees, Node* speciesTree, string species_separator, int species_index)
{
    //compute gene leaf to species leaf mapping.  Could be faster by using a map for the species leaves by name...
    unordered_map<Node*, Node*> geneSpeciesMapping;

    for (int i = 0; i < geneTrees.size(); i++)
    {
        unordered_map<Node*, Node*> tmpmap = GeneSpeciesTreeUtil::Instance()->GetGeneSpeciesMappingByLabel(geneTrees[i], speciesTree, species_separator, species_index);

        for (unordered_map<Node*, Node*>::iterator it = tmpmap.begin(); it != tmpmap.end(); ++it)
        {
            geneSpeciesMapping[ it->first ] = it->second;
        }
    }

    return geneSpeciesMapping;
}



void PrintHelp()
{
    cout    <<"------------------------------------------------------------------"<<endl
            <<"MULTREC - Multi-reconciliation program "<<endl
            <<"------------------------------------------------------------------"<<endl
            <<"Multrec takes as input a species tree S, a set of gene trees, a duplication cost, a loss cost and a parameter duplication height h.  The output is a mapping of the gene tree nodes to S that minimizes the segmental reconciliation cost, assuming that there exists such a mapping that has duplications sum-of-heights at most h.  If loss cost >= dup cost, the LCA mapping is returned."<<endl
            <<"The leaves of the gene trees must map to the leaves of S.  The gene tree leaves are assumed to have the format [species_name]__[gene_name], for example HUMAN_BRCA2 indicates that the gene is mapped to the leaf of S names HUMAN.  The gene/species separator can be changed with the -spsep argument, and the position of the species name in the gene name with the -spindex argument, indexed at 0.  "<<endl
            <<"If your genes are name e.g. GENENAME_SPECIESNAME_OTHERSTUFF, you can set -spsep \"_\" -spindex 1"<<endl
            <<endl
            <<"The format of the output is 5 lines, as follows"<<endl
            <<"COST=[total cost of mapping]"<<endl
            <<"DUPHEIGHT=[sum of duplication heights]"<<endl
            <<"NBLOSSES=[number of losses]"<<endl
            <<"SPECIESTREE=[species tree newick, with internal nodes labeled"<<endl
            <<"GENETREES=[all gene tree newick separated by ;.  Internal nodes are labeled by mapping]"<<endl
            <<endl
            <<"If no solution is found (when h is too small), then the output is simply"<<endl
            <<"NO SOLUTION FOUND"<<endl
            <<endl
            <<"Sample command line:"<<endl
            <<"./Multrec -d 10 -l 3 -gf ./sample_data/geneTrees.txt -sf ./sample_data/speciesTree.txt"
            <<endl
            <<"Required arguments:"<<endl
            <<"At least one of -g or -gf must be specified, and at least one of -s or -sf must be specified."<<endl
            <<"-g   [g1;g2;...;gk]   Here g1,g2,...,gk are gene trees"<<endl
            <<"                      represented in Newick format.  "<<endl
            <<"                      The gene trees are separated by the ; symbol.	"<<endl
            <<"-gf  [file]           file is the name of a file containing the list "<<endl
            <<"                      of gene trees, all in Newick format and separated "<<endl
            <<"                      by a ; symbol in the file."<<endl
            <<"-s   [newick]         The species tree in Newick format."<<endl
            <<"-sf  [file]           Name of the file containing species tree Newick."<<endl
            <<""<<endl
            <<"Optional arguments:"<<endl
            <<"--help                Print this help message."<<endl
            <<"-d   [double]         The cost for one height of duplication.  Default=3"<<endl
            <<"-l   [double]         The cost for one loss.  Default=1"<<endl
            <<"-h   [int]            Maximum allowed duplication sum-of-heights.  Default=20"<<endl
            <<"-o   [file]           Output file.  Default=output to console"<<endl
            <<"-spsep   [string]     Gene/species separator in the gene names.  Default=__"<<endl
            <<"-spindex [int]        Position of the species in the gene names, after "<<endl
            <<"                      being split by the gene/species separator.  Default=0"<<endl
            <<"--test                Launches a series of unit tests.  This includes small fixed "<<endl
            <<"                      examples with known outputs to expect, and larger random trees "<<endl
            <<"                      to see if the program terminates in an OK status on more complicated"<<endl
            <<"                      datasets.  "<<endl;
}



MultiGeneReconcilerInfo Execute(map<string, string> args)
{
    MultiGeneReconcilerInfo info;
    info.isBad = true;

    vector<Node*> geneTrees;
    Node* speciesTree = NULL;

    string species_separator = "__";
    int species_index = 0;
    double dupcost = 2;
    double losscost = 1;
    int maxDupheight = 20;

    //parse dup loss cost and max dup height
    if (args.find("d") != args.end())
    {
        dupcost = Util::ToDouble(args["d"]);
    }
    if (args.find("l") != args.end())
    {
        losscost = Util::ToDouble(args["l"]);
    }
    if (args.find("h") != args.end())
    {
        maxDupheight = Util::ToInt(args["h"]);
    }

    string outfile = "";
    if (args.find("o") != args.end())
    {
        outfile = args["o"];
    }
    //parse gene trees, either from command line or from file
    if (args.find("g") != args.end())
    {
        vector<string> gstrs = Util::Split( Util::ReplaceAll(args["g"], "\n", ""), ";", false);

        for (int i = 0; i < gstrs.size(); i++)
        {
            string str = gstrs[i];
            Node* tree = NewickLex::ParseNewickString(str, false);

            if (!tree)
            {
                cout<<"Error: there is a problem with input gene tree "<<str<<endl;
                return info;
            }

            geneTrees.push_back(tree);
        }
    }
    else if (args.find("gf") != args.end())
    {
        string gcontent = Util::GetFileContent(args["gf"]);

        vector<string> lines = Util::Split( gcontent, "\n");
        vector<string> gstrs;
        for (int l = 0; l < lines.size(); l++)
        {
            vector<string> trees_on_line = Util::Split( lines[l], ";", false);
            for (int t = 0; t < trees_on_line.size(); t++)
            {
                if (trees_on_line[t] != "")
                    gstrs.push_back(trees_on_line[t]);
            }
        }


        for (int i = 0; i < gstrs.size(); i++)
        {
            string str = gstrs[i];
            Node* tree = NewickLex::ParseNewickString(str, false);

            if (!tree)
            {
                cout<<"Error: there is a problem with input gene tree "<<str<<endl;
                return info;
            }

            geneTrees.push_back(tree);
        }
    }


    //parse species trees, either from command line or from file
    if (args.find("s") != args.end())
    {
        speciesTree = NewickLex::ParseNewickString(args["s"], false);

        if (!speciesTree)
        {
            cout<<"Error: there is a problem with the species tree."<<endl;
            return info;
        }
    }
    else if (args.find("sf") != args.end())
    {
        string scontent = Util::GetFileContent(args["sf"]);
        speciesTree = NewickLex::ParseNewickString(scontent, false);

        if (!speciesTree)
        {
            cout<<"Error: there is a problem with the species tree."<<endl;
            return info;
        }
    }




    //parse species separator and index
    if (args.find("spsep") != args.end())
    {
        species_separator = args["spsep"];
    }

    if (args.find("spindex") != args.end())
    {
        species_index = Util::ToInt(args["spindex"]);
    }




    if (geneTrees.size() == 0)
    {
        cout<<"No gene tree given.  Program will exit."<<endl;
        PrintHelp();
        return info;
    }
    if (!speciesTree)
    {
        cout<<"No species tree given.  Program will exit."<<endl;
        PrintHelp();
        return info;
    }

    if (dupcost < 0 || losscost <= 0)
    {
        cout<<"dupcost < 0 or losscost <= 0 are prohibited.  Program will exit."<<endl;
        return info;
    }


    /*if (losscost >= dupcost)
    {
        //do usual lca mapping stuff
    }
    else*/
    {

        if (dupcost/losscost > 20)
        {
            cout<<"WARNING: dupcost/losscost > 20 or losscost = 0.  Unless your trees are small, the program may not finish before the sun has grown large enough to gobble the earth."<<endl;
        }

        GeneSpeciesTreeUtil::Instance()->LabelInternalNodesUniquely(speciesTree);


        unordered_map<Node*, Node*> geneSpeciesMapping = GetGeneSpeciesMapping(geneTrees, speciesTree, species_separator, species_index);

        MultiGeneReconciler reconciler(geneTrees, speciesTree, geneSpeciesMapping, dupcost, losscost, maxDupheight);

        info = reconciler.Reconcile();

        string output = "";
        if (info.isBad)
        {
            output = "NO SOLUTION FOUND";
        }
        else
        {
            output += "<COST>\n" + Util::ToString(info.GetCost(dupcost, losscost)) + "\n</COST>\n";
            output += "<DUPHEIGHT>\n" + Util::ToString(info.dupHeightSum) + "\n</DUPHEIGHT>\n";
            output += "<NBLOSSES>\n" + Util::ToString(info.nbLosses) + "\n</NBLOSSES>\n";
            output += "<SPECIESTREE>\n" + NewickLex::ToNewickString(speciesTree) + "\n</SPECIESTREE>\n";

            map<Node*, vector< pair<int, Node*> > > dups_per_species = LabelGeneTreesWithSpeciesMapping(geneTrees, speciesTree, reconciler, info, false);

            output += "<GENETREES>\n";
            for (int t = 0; t < geneTrees.size(); t++)
            {
                output += NewickLex::ToNewickString(geneTrees[t]) + "\n";
            }
            output += "</GENETREES>\n";

            output += "<DUPS_PER_SPECIES>\n";
            TreeIterator* itsp = speciesTree->GetPostOrderIterator();
            while (Node* s = itsp->next())
            {
                if (dups_per_species.find(s) != dups_per_species.end())
                {
                    output += "[" + s->GetLabel() + "] ";
                    vector< pair<int, Node*> > dups_for_s = dups_per_species[s];

                    for (int d = 0; d < dups_for_s.size(); d++)
                    {
                        pair<int, Node*> p = dups_for_s[d];
                        string lbl = p.second->GetLabel();
                        lbl = Util::GetSubstringAfter(lbl, "_");

                        output += lbl + " (G" + Util::ToString(p.first) + ") ";

                    }
                    output += "\n";
                }
            }
            speciesTree->CloseIterator(itsp);
            output += "</DUPS_PER_SPECIES>\n";
        }

        if (outfile == "")
            cout<<output;
        else
            Util::WriteFileContent(outfile, output);

    }






    for (int i = 0; i < geneTrees.size(); i++)
    {
        delete geneTrees[i];
    }

    delete speciesTree;


    return info;
}




bool RunTest(vector<Node*> geneTrees, Node* speciesTree, unordered_map<Node*, Node*> geneSpeciesMapping,
             double dupcost, double losscost, int maxDupHeight,
             int expectedDupHeightSum, int expectedNbLosses, bool isExpectedBad, bool detailed = false)
{
    bool ok = true;

    MultiGeneReconciler reconciler(geneTrees, speciesTree, geneSpeciesMapping, dupcost, losscost, maxDupHeight);
    MultiGeneReconcilerInfo info = reconciler.Reconcile();



    if (detailed && !info.isBad)
    {
        LabelGeneTreesWithSpeciesMapping(geneTrees, speciesTree, reconciler, info);
        cout<<NewickLex::ToNewickString(speciesTree)<<endl;
        cout<<NewickLex::ToNewickString(geneTrees[0])<<endl;
        cout<<NewickLex::ToNewickString(geneTrees[1])<<endl;
        cout<<info.dupHeightSum<<" dups + "<<info.nbLosses<<" losses"<<endl;
    }

    if (!isExpectedBad)
    {
        if (info.isBad)
        {
            ok = false;
            cout<<"FAILED: info is bad and I don't know why."<<endl;
        }
        if (info.dupHeightSum != expectedDupHeightSum)
        {
            ok = false;
            cout<<"FAILED: dup height sum should be "<<expectedDupHeightSum<<" (not "<<info.dupHeightSum<<")"<<endl;
        }
        if (info.nbLosses != expectedNbLosses)
        {
            ok = false;
            cout<<"FAILED: losses should be "<<expectedNbLosses<<" (not "<<info.nbLosses<<")"<<endl;
        }
        double cost = reconciler.GetMappingCost(info.partialMapping);
        double cost2 = info.dupHeightSum * dupcost + info.nbLosses * losscost;
        if (cost - cost2 > 0.0000001)
        {
            ok = false;
            cout<<"FAILED: reconciler score does not match computed score ("
                <<cost<<" vs "<<cost2<<")"<<endl;
        }
    }
    else
    {
        if (!info.isBad)
        {
            ok = false;
            cout<<"FAILED: info is not bad but it should be..."<<endl;
        }
    }

    return ok;
}




bool TestRandomTrees()
{
    int nbTests = 1;
    int nbOK = 0;

    cout<<endl<<"*** TestRandomTrees ***"<<endl;

    cout<<endl<<"(testing "<<nbTests<<" random instance(s) - this might take a few minutes)"<<endl;

    for (int i = 0; i < nbTests; i++)
    {
        bool ok = true;
        //generate random sptree
        Node* sptree = new Node(false);
        int sleaves = rand() % 15 + 10;
        for (int l = 0; l < sleaves; l++)
        {
            Node* s = sptree->AddChild();
            s->SetLabel(Util::ToString(l));
        }
        sptree->BinarizeRandomly();

        //generate random gene trees
        int nbgeneTrees = rand() % 10 + 10;
        vector<Node*> geneTrees;
        for (int t = 0; t < nbgeneTrees; t++)
        {
            Node* genetree = new Node(false);
            int gleaves = rand() % (int)(sleaves * 2.5) + 5;
            for (int j = 0; j < gleaves; j++)
            {
                Node* g = genetree->AddChild();
                int spindex = rand() % sleaves;
                g->SetLabel(Util::ToString(spindex) + "__" + Util::ToString(j));
            }
            genetree->BinarizeRandomly();
            geneTrees.push_back(genetree);
        }

        GeneSpeciesTreeUtil::Instance()->LabelInternalNodesUniquely(sptree);

        unordered_map<Node*, Node*> gsMapping = GetGeneSpeciesMapping(geneTrees, sptree, "__", 0);

        MultiGeneReconciler reconciler(geneTrees, sptree, gsMapping, 1, 1, 1000);
        MultiGeneReconcilerInfo info = reconciler.Reconcile();

        cout<<"TEST "<<i + 1<<" : random trees: nbSpecies="<<sleaves<<" nbGeneTrees="<<geneTrees.size()<<endl;

        if (info.isBad)
        {
            cout<<"FAILED: lca mapping is bad"<<endl;
            ok = false;
        }
        else
        {
            cout<<"PASSED LCA MAPPING"<<endl;

            cout<<"Testing DUP2 with maxheight="<<min(30, info.dupHeightSum)<<endl;

            MultiGeneReconciler reconciler_2(geneTrees, sptree, gsMapping, 2, 1, min(30, info.dupHeightSum));
            MultiGeneReconcilerInfo info_2 = reconciler_2.Reconcile();

            if (info_2.isBad)
            {
                cout<<"FAILED: dupcost 2 is bad"<<endl;
                ok = false;
            }
            else if (info.dupHeightSum < 30 && info_2.dupHeightSum > info.dupHeightSum)
            {
                cout<<"FAILED: dupHeightSum 2 > dupHeightSum 1"<<endl;
                ok = false;
            }
            else
            {
                cout<<"PASSED DUPS2 TEST"<<endl;
            }

            cout<<"Testing DUP5 with maxheight="<<min(10, info.dupHeightSum)<<endl;

            MultiGeneReconciler reconciler_3(geneTrees, sptree, gsMapping, 5, 1, min(10, info.dupHeightSum));
            MultiGeneReconcilerInfo info_3 = reconciler_3.Reconcile();

            cout<<"PASSED DUPS5 TEST (hey, it terminated!)"<<endl;


        }



        delete sptree;
        for (int j = 0; j < geneTrees.size(); j++)
        {
            delete geneTrees[j];
        }

        if (ok)
            nbOK++;


    }

    cout<<"TOTAL = "<<nbOK<<"/"<<nbTests<<endl;
}



void TestBasicInstance()
{
    cout<<endl<<"*** TestBasicInstance ***"<<endl;

    map<string, string> args;

    string g1 = "((A__1, C__1),B__1);";
    string g2 = "((A__2, B__2),B__3);";
    string snewick = "((A,B),(C,D));";

    vector<Node*> geneTrees;
    geneTrees.push_back( NewickLex::ParseNewickString( g1 ) );
    geneTrees.push_back( NewickLex::ParseNewickString( g2 ) );

    Node* speciesTree = NewickLex::ParseNewickString(snewick);
    GeneSpeciesTreeUtil::Instance()->LabelInternalNodesUniquely(speciesTree);

    unordered_map<Node*, Node*> gsMapping = GetGeneSpeciesMapping(geneTrees, speciesTree, "__", 0);


    int nbOK = 0;
    int nbTests = 0;

    cout<<"TEST 1: delta = 0.2, lambda = 10 (LCA mapping)"<<endl;
    nbTests++;
    bool ok = RunTest(geneTrees, speciesTree, gsMapping, .2, 10, 20, 2, 5, false);
    if (ok){ nbOK++;  cout<<"PASSED!"<<endl; }


    cout<<"TEST 2: delta = 1 (LCA mapping)"<<endl;
    nbTests++;
    ok = RunTest(geneTrees, speciesTree, gsMapping, 1, 1, 20, 2, 5, false);
    if (ok){ nbOK++;  cout<<"PASSED!"<<endl; }

    cout<<"TEST 3: delta = 2.0001"<<endl;
    nbTests++;
    ok = RunTest(geneTrees, speciesTree, gsMapping, 2.0001, 1, 2, 1, 7, false);
    if (ok){ nbOK++;  cout<<"PASSED!"<<endl; }

    cout<<"TOTAL = "<<nbOK<<"/"<<nbTests<<endl;

    for (int i = 0; i < geneTrees.size(); i++)
    {
        delete geneTrees[i];
    }

    delete speciesTree;
}



void TestCaterpillarSpeciesTree()
{
    cout<<endl<<"*** TestCaterpillarSpeciesTree ***"<<endl<<endl;

    map<string, string> args;

    vector<string> slbls;
    vector<string> g1_labels;
    vector<string> g2_labels;
    vector<string> g3_labels;

    for (int i = 1; i <= 10; i++)
    {
        slbls.push_back(Util::ToString(i));
    }


    g1_labels.push_back("1");
    for (int i = 1; i <= 8; i++)
    {
        g1_labels.push_back("6");
    }


    g2_labels.push_back("1");
    for (int i = 2; i <= 5; i++)
    {
        g2_labels.push_back(Util::ToString(i));
        g2_labels.push_back(Util::ToString(i));
    }


    string snewick = NewickLex::GetCaterpillarNewick(slbls);
    string g1_newick = NewickLex::GetCaterpillarNewick(g1_labels);
    string g2_newick = NewickLex::GetCaterpillarNewick(g2_labels);


    vector<Node*> geneTrees;
    geneTrees.push_back( NewickLex::ParseNewickString( g1_newick ) );
    geneTrees.push_back( NewickLex::ParseNewickString( g2_newick ) );

    Node* speciesTree = NewickLex::ParseNewickString(snewick);
    GeneSpeciesTreeUtil::Instance()->LabelInternalNodesUniquely(speciesTree);

    unordered_map<Node*, Node*> gsMapping = GetGeneSpeciesMapping(geneTrees, speciesTree, "__", 0);


    int nbOK = 0;
    int nbTests = 0;

    bool ok = true;





    cout<<"Test 1: delta = 1.999"<<endl;
    nbTests++;
    ok = RunTest(geneTrees, speciesTree, gsMapping, 1.999, 1, 20, 11, 15, false);
    if (ok) {nbOK++; cout<<"PASSED!"<<endl;}

    cout<<"Test 2: delta = 2.0001"<<endl;
    nbTests++;
    ok = RunTest(geneTrees, speciesTree, gsMapping, 2.0001, 1, 20, 10, 17, false);
    if (ok) {nbOK++; cout<<"PASSED!"<<endl;}

    cout<<"Test 3: delta = 3.9999"<<endl;
    nbTests++;
    ok = RunTest(geneTrees, speciesTree, gsMapping, 3.9999, 1, 20, 10, 17, false);
    if (ok) {nbOK++; cout<<"PASSED!"<<endl;}

    cout<<"Test 4: delta = 5.0001"<<endl;
    nbTests++;
    ok = RunTest(geneTrees, speciesTree, gsMapping, 5.0001, 1, 20, 9, 22, false);
    if (ok) {nbOK++; cout<<"PASSED!"<<endl;}


    cout<<"Test 5: delta = 23/2 + 0.00001"<<endl;
    nbTests++;
    ok = RunTest(geneTrees, speciesTree, gsMapping, 23/2 + 0.00001, 1, 20, 7, 38, false);
    if (ok) {nbOK++; cout<<"PASSED!"<<endl;}


    cout<<"Test 6: delta = 1.00001 but maxdupheight = 7"<<endl;
    nbTests++;
    ok = RunTest(geneTrees, speciesTree, gsMapping, 1.00001, 1, 7, 7, 38, true);
    if (ok) {nbOK++; cout<<"PASSED!"<<endl;}


    cout<<"Test 7: maxdupheight = 6, should be bad"<<endl;
    nbTests++;
    ok = RunTest(geneTrees, speciesTree, gsMapping, 20, 1, 6, 7, 38, true);
    if (ok) {nbOK++; cout<<"PASSED!"<<endl;}



    cout<<"Test 8: delta = 0.2, lambda = 10"<<endl;
    nbTests++;
    ok = RunTest(geneTrees, speciesTree, gsMapping, .2, 10, 20, 11, 15, false);
    if (ok) {nbOK++; cout<<"PASSED!"<<endl;}



    cout<<"Test 9: delta = 5.0001 * 10, losses = 1 * 10"<<endl;
    nbTests++;
    ok = RunTest(geneTrees, speciesTree, gsMapping, 5.0001 * 10, 1 * 10, 20, 9, 22, false);
    if (ok) {nbOK++; cout<<"PASSED!"<<endl;}


    cout<<"TOTAL = "<<nbOK<<"/"<<nbTests<<endl;


    for (int i = 0; i < geneTrees.size(); i++)
    {
        delete geneTrees[i];
    }

    delete speciesTree;
}









int main(int argc, char *argv[])
{

    map<string, string> args;

    bool hasHelp = false;
    bool hasTest = false;

    //BUILD DICTIONARY OF ARGS
    string prevArg = "";
    for (int i = 0; i < argc; i++)
    {
        if (string(argv[i]) == "-v")
        {
            verbose = 1;
            prevArg = "";
        }
        else if (string(argv[i]) == "--help")
        {
            hasHelp = true;
        }
        else if (string(argv[i]) == "--test")
        {
            hasTest = true;
        }
        else
        {
            if (prevArg != "" && prevArg[0] == '-')
            {
                args[Util::ReplaceAll(prevArg, "-", "")] = string(argv[i]);
            }

            prevArg = string(argv[i]);
        }
    }

    //args["test"] = 1;

    if (args.find("test") != args.end() || hasTest)
    {
        TestBasicInstance();
        TestCaterpillarSpeciesTree();
        TestRandomTrees();

        return 0;
    }
    else if (args.find("help") != args.end() || hasHelp)
    {
        PrintHelp();
    }
    else
    {
        //ML's ad-hoc testing stuff for Windows.
        /*args["d"] = "10";
        args["l"] = "0.1";
        args["gf"] = "W:/Users/Manuel/Documents/GitHub/Multrec/Multrec/Multrec/sample_data/geneTrees.txt";
        args["sf"] = "W:/Users/Manuel/Documents/GitHub/Multrec/Multrec/Multrec/sample_data/speciesTree.txt";
        args["o"] = "W:/Users/Manuel/Desktop/tmp/out.txt";*/


        Execute(args);
        return 0;
    }
}


