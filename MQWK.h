//
// Created by Guangnian_DONG on 2018/10/14 014.
//
#ifndef WQRTQ_MQWK_H
#define WQRTQ_MQWK_H

#include <iostream>
#include <vector>
#include <climits>
#include <cfloat>
#include <tuple>
#include "operator.h"
#include "MQP.h"
#include "MWK.h"
#include "minheap.h"
#include "QuadProg/QuadProg++.hh"
#include "RTree/RTree.h"

#define RTREE_TEMPLATE template<class DATATYPE, class ELEMTYPE, int NUMDIMS, class ELEMTYPEREAL, int TMAXNODES, int TMINNODES>
#define RTREE_QUAL RTree<DATATYPE, ELEMTYPE, NUMDIMS, ELEMTYPEREAL, TMAXNODES, TMINNODES>

RTREE_TEMPLATE
std::tuple<DataPoint, std::vector<WeightPoint>, int> MQWK(RTREE_QUAL P, DataPoint q,unsigned int k, std::vector<WeightPoint> Wm, unsigned int dataSampleSize, unsigned int weightSampleSize, double alpha, double beta, double gamma, double lambda);
double Penalty(DataPoint q, DataPoint qapos, std::vector<WeightPoint> Wm, std::vector<WeightPoint> Wmapos, int k, int kapos, double alpha, double beta, double gamma, double lambda);


/**input (1)dataset P on R-tree
 *       (2)query point q
 *       (3)k
 *       (4)why-not weighting std::vector set Wm
 *       (5)sample size |S|
 *       (6)sample size |Q|
 * output(1)modified W'm
 *       (2)modified k'
 *       (3)modified q'
 */
RTREE_TEMPLATE
std::tuple<DataPoint, std::vector<WeightPoint>, int> MQWK(RTREE_QUAL P, DataPoint q,unsigned int k, std::vector<WeightPoint> Wm, unsigned int dataSampleSize, unsigned int weightSampleSize, double alpha, double beta, double gamma, double lambda)
{
    std::vector<WeightPoint> Wmapos;
    int kapos;
    DataPoint qapos;
    
    std::vector<DataPoint> Q(dataSampleSize);
    DataPoint qmin = MQP(P, q, k, Wm);
    
    //sample |Q| query points from the space determined by qmin and q, and added them to Q
    srand(time(0));
    for(int i = 0; i < weightSampleSize; i++)
    {
        const int dim = q.getDim();
        double point[dim];
        for(int j = 0; j < dim; j++)
        {
            point[j] = ((double)rand()/RAND_MAX)*(q[j] - qmin[j]) + qmin[j];///此处可能有问题，随机数可能并不随机
        }
        DataPoint p(point, dim, 0);
        Q[i] = p;
    }
    
    double MinPenalty = DBL_MAX;
    for(auto qi = Q.begin(); qi != Q.end(); qi++)
    {
        std::vector<WeightPoint> Wmaa;
        int kaa;
        std::tie(Wmaa, kaa) = MWK(P, *qi, k, Wm, dataSampleSize, alpha, beta);
        double newPenalty = Penalty(q, *qi, Wm, Wmaa, k, kaa, alpha, beta, gamma, lambda);
        if(newPenalty < MinPenalty)
        {
            qapos = *qi;
            kapos = kaa;
            Wmapos = Wmaa;
        }
    }
    return std::make_tuple(qapos, Wmapos, kapos);
}


/**
 * @param q
 * @param qapos
 * @param Wm
 * @param Wmapos
 *
 * @param k
 * @param kapos
 * @param alpha
 * @param beta
 * @param gamma
 * @param lambda
 * @return
 */
double Penalty(DataPoint q, DataPoint qapos, std::vector<WeightPoint> Wm, std::vector<WeightPoint> Wmapos, int k, int kapos, double alpha, double beta, double gamma, double lambda)
{
    double peank, peanWm, peanq;
    //计算q
    double molecule_q = 0;
    double denominator_q = 0;
    for(int i = 0; i < q.getDim(); i++)
    {
        molecule_q += pow(q[i] - qapos[i], 2);
        denominator_q += pow(q[i], 2);
    }
    peanq = sqrt(molecule_q)/sqrt(denominator_q);
    //计算k
    double molecule_k = kapos>k?kapos-k:0;
    double denomonator_k = Wmapos[0].getRank();
    for(auto w = Wmapos.begin(); w != Wmapos.end(); w++)
    {
        if(denomonator_k < w->getRank())
        {
            denomonator_k = w->getRank();
        }
    }
    denomonator_k -= k;
    peank = molecule_k/denomonator_k;
    //计算Wm
    double molecule_Wm = 0;
    double denominator_Wm = 0;
    auto w = Wm.begin();
    auto wapos = Wmapos.begin();
    while(w != Wm.end() && wapos != Wmapos.end())
    {
        int d = w->getDim();
        double mole = 0;
        double deno = 0;
        for(int i = 0; i < d; i++)
        {
            mole += pow(w->getData()[i] - wapos->getData()[i], 2);
            deno += pow(wapos->getData()[i], 2);
        }
        molecule_Wm += sqrt(mole);
        denominator_Wm += sqrt(1+deno);
        w++;
        wapos++;
    }
    peanWm = molecule_Wm/denominator_Wm;
    
    return gamma*peanq + lambda*(alpha*peank + beta*peanWm);
}

#endif