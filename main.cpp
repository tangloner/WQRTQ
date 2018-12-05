#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include "QuadProg/QuadProg++.hh"
#include "RTree/RTree.h"
#include "minheap.h"
#include "operator.h"
#include "MQP.h"
#include "MWK.h"
#include "MQWK.h"

using namespace std;

void testRTree()
{
    std::cout<<"test"<<std::endl;
    test();
    std::cout<<"test_memory"<<std::endl;
    MemoryTest();
}

void testQuadProg()
{
    std::cout<<"test_quadprog"<<std::endl;
    quadprogpp::test_quadprog();
}

void testMQP()
{
    int k = 3;
    
    cout<<"testMQP begin"<<endl;
    double p01[] = {2,1};
    double p02[] = {6,3};
    double p03[] = {1,9};
    double p04[] = {9,3};
    double p05[] = {7,5};
    double p06[] = {5,8};
    double p07[] = {3,7};
    double q0[] = {4,4};
    DataPoint p1(p01, 2, 1);
    DataPoint p2(p02, 2, 2);
    DataPoint p3(p03, 2, 3);
    DataPoint p4(p04, 2, 4);
    DataPoint p5(p05, 2, 5);
    DataPoint p6(p06, 2, 6);
    DataPoint p7(p07, 2, 7);
    DataPoint q(q0, 2, 0);
    //DATATYPE ELEMTYPE DIMNUM REALELEMTYPE
    RTree<int, double, 2, double > tree;

    tree.Insert(q.rect_min(), q.rect_max(),q.getIndex());
    tree.Insert(p1.rect_min(), p1.rect_max(),p1.getIndex());
    tree.Insert(p2.rect_min(), p2.rect_max(),p2.getIndex());
    tree.Insert(p3.rect_min(), p3.rect_max(),p3.getIndex());
    tree.Insert(p4.rect_min(), p4.rect_max(),p4.getIndex());
    tree.Insert(p5.rect_min(), p5.rect_max(),p5.getIndex());
    tree.Insert(p6.rect_min(), p6.rect_max(),p6.getIndex());
    tree.Insert(p7.rect_min(), p7.rect_max(),p7.getIndex());

    cout<<"RTree in main: "<<endl;
    int itIndex = 0;
    RTree<int, double, 2, double>::Iterator it;
    for( tree.GetFirst(it); !tree.IsNull(it); tree.GetNext(it) )
    {
        int value = tree.GetAt(it);
        double boundsMin[2] = {0,0};
        double boundsMax[2] = {0,0};
        it.GetBounds(boundsMin, boundsMax);
        printf("it[%d] %d = (%f,%f,%f,%f)\n", itIndex++, value, boundsMin[0], boundsMin[1], boundsMax[0], boundsMax[1]);
    }

    double w01[] = {0.1, 0.9};
    double w02[] = {0.3, 0.7};
    double w03[] = {0.5, 0.5};
    double w04[] = {0.9, 0.1};
    WeightPoint w1(w01, 2, 1);
    WeightPoint w2(w02, 2, 2);
    WeightPoint w3(w03, 2, 3);
    WeightPoint w4(w04, 2, 4);
    vector<WeightPoint> Wm(2);
    Wm[0] = w1;
    Wm[1] = w4;
    cout<<endl;
    cout<<"Wm in main: "<<endl;
    for(auto iter = Wm.begin(); iter != Wm.end(); iter++)
    {
        for(int i = 0; i < iter->size(); i++)
        {
            cout<<(*iter)[i]<<",";
        }
        cout<<endl;
    }

    DataPoint q_m = MQP(tree, q, k, Wm);
    cout<<endl;
    q_m.toString();
    cout<<endl;
    cout<<"testMQP finished"<<endl;
}

void testMWK()
{
    int k = 3;
    double alpha = 0.5;
    double beta = 0.5;
    double sampleSize = 20;
    
    cout<<"testMQP begin"<<endl;
    double p01[] = {2,1};
    double p02[] = {6,3};
    double p03[] = {1,9};
    double p04[] = {9,3};
    double p05[] = {7,5};
    double p06[] = {5,8};
    double p07[] = {3,7};
    double q0[] = {4,4};
    DataPoint p1(p01, 2, 1);
    DataPoint p2(p02, 2, 2);
    DataPoint p3(p03, 2, 3);
    DataPoint p4(p04, 2, 4);
    DataPoint p5(p05, 2, 5);
    DataPoint p6(p06, 2, 6);
    DataPoint p7(p07, 2, 7);
    DataPoint q(q0, 2, 0);
    //DATATYPE ELEMTYPE DIMNUM REALELEMTYPE
    RTree<int, double, 2, double > tree;
    
    tree.Insert(q.rect_min(), q.rect_max(),q.getIndex());
    tree.Insert(p1.rect_min(), p1.rect_max(),p1.getIndex());
    tree.Insert(p2.rect_min(), p2.rect_max(),p2.getIndex());
    tree.Insert(p3.rect_min(), p3.rect_max(),p3.getIndex());
    tree.Insert(p4.rect_min(), p4.rect_max(),p4.getIndex());
    tree.Insert(p5.rect_min(), p5.rect_max(),p5.getIndex());
    tree.Insert(p6.rect_min(), p6.rect_max(),p6.getIndex());
    tree.Insert(p7.rect_min(), p7.rect_max(),p7.getIndex());
    
    cout<<"RTree in main: "<<endl;
    int itIndex = 0;
    RTree<int, double, 2, double>::Iterator it;
    for( tree.GetFirst(it); !tree.IsNull(it); tree.GetNext(it) )
    {
        int value = tree.GetAt(it);
        double boundsMin[2] = {0,0};
        double boundsMax[2] = {0,0};
        it.GetBounds(boundsMin, boundsMax);
        printf("it[%d] %d = (%f,%f,%f,%f)\n", itIndex++, value, boundsMin[0], boundsMin[1], boundsMax[0], boundsMax[1]);
    }
    
    double w01[] = {0.1, 0.9};
    double w02[] = {0.3, 0.7};
    double w03[] = {0.5, 0.5};
    double w04[] = {0.9, 0.1};
    WeightPoint w1(w01, 2, 1);
    WeightPoint w2(w02, 2, 2);
    WeightPoint w3(w03, 2, 3);
    WeightPoint w4(w04, 2, 4);
    vector<WeightPoint> Wm(2);
    Wm[0] = w1;
    Wm[1] = w4;
    
    vector<WeightPoint> Wmapos;
    int kapos;
        tie(Wmapos, kapos) = MWK(tree, q, k, Wm, sampleSize, alpha, beta);
    cout<<endl;
    cout<<"W'm = "<<endl;
    for(auto w = Wmapos.begin(); w != Wmapos.end(); w++)
    {
        w->toString();
    }
    cout<<endl;
    cout<<"k' = "<<kapos<<endl;
}

void testMQWK()
{
    int k = 3;
    double alpha = 0.5;
    double beta = 0.5;
    double gamma = 0.5;
    double lambda = 0.5;
    int dataSampleSize = 20;
    int weightSampleSize = 20;
    
    
    cout<<"testMQP begin"<<endl;
    double p01[] = {2,1};
    double p02[] = {6,3};
    double p03[] = {1,9};
    double p04[] = {9,3};
    double p05[] = {7,5};
    double p06[] = {5,8};
    double p07[] = {3,7};
    double q0[] = {4,4};
    DataPoint p1(p01, 2, 1);
    DataPoint p2(p02, 2, 2);
    DataPoint p3(p03, 2, 3);
    DataPoint p4(p04, 2, 4);
    DataPoint p5(p05, 2, 5);
    DataPoint p6(p06, 2, 6);
    DataPoint p7(p07, 2, 7);
    DataPoint q(q0, 2, 0);
    //DATATYPE ELEMTYPE DIMNUM REALELEMTYPE
    RTree<int, double, 2, double > tree;
    
    tree.Insert(q.rect_min(), q.rect_max(),q.getIndex());
    tree.Insert(p1.rect_min(), p1.rect_max(),p1.getIndex());
    tree.Insert(p2.rect_min(), p2.rect_max(),p2.getIndex());
    tree.Insert(p3.rect_min(), p3.rect_max(),p3.getIndex());
    tree.Insert(p4.rect_min(), p4.rect_max(),p4.getIndex());
    tree.Insert(p5.rect_min(), p5.rect_max(),p5.getIndex());
    tree.Insert(p6.rect_min(), p6.rect_max(),p6.getIndex());
    tree.Insert(p7.rect_min(), p7.rect_max(),p7.getIndex());
    
    cout<<"RTree in main: "<<endl;
    int itIndex = 0;
    RTree<int, double, 2, double>::Iterator it;
    for( tree.GetFirst(it); !tree.IsNull(it); tree.GetNext(it) )
    {
        int value = tree.GetAt(it);
        double boundsMin[2] = {0,0};
        double boundsMax[2] = {0,0};
        it.GetBounds(boundsMin, boundsMax);
        printf("it[%d] %d = (%f,%f,%f,%f)\n", itIndex++, value, boundsMin[0], boundsMin[1], boundsMax[0], boundsMax[1]);
    }
    
    double w01[] = {0.1, 0.9};
    double w02[] = {0.3, 0.7};
    double w03[] = {0.5, 0.5};
    double w04[] = {0.9, 0.1};
    WeightPoint w1(w01, 2, 1);
    WeightPoint w2(w02, 2, 2);
    WeightPoint w3(w03, 2, 3);
    WeightPoint w4(w04, 2, 4);
    vector<WeightPoint> Wm(2);
    Wm[0] = w1;
    Wm[1] = w4;
    
    DataPoint qapos;
    vector<WeightPoint> Wmapos;
    int kapos;
    tie(qapos, Wmapos, kapos) = MQWK(tree, q, k, Wm, dataSampleSize, weightSampleSize, alpha, beta, gamma, lambda);
    cout<<endl;
    cout<<"q = "<<endl;
    q.toString();
    cout<<endl;
    cout<<"W'm = "<<endl;
    for(auto w = Wmapos.begin(); w != Wmapos.end(); w++)
    {
        w->toString();
    }
    cout<<endl;
    cout<<"k' = "<<kapos<<endl;
}



int main(int argc, char* argv[]) {
//    testRTree();
//    testQuadProg();
//    testMinHeap();
//    learnVector();
//    testMQP();
    cout<<"hello world"<<endl;
    return 0;
}
