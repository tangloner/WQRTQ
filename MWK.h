//
// Created by Guangnian_DONG on 2018/10/14 014.
//
#ifndef WQRTQ_MWK_H
#define WQRTQ_MWK_H

#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <ctime>
#include <tuple>
#include "operator.h"
#include "minheap.h"
#include "QuadProg/QuadProg++.hh"
#include "RTree/RTree.h"

#define RTREE_TEMPLATE template<class DATATYPE, class ELEMTYPE, int NUMDIMS, class ELEMTYPEREAL, int TMAXNODES, int TMINNODES>
#define RTREE_QUAL RTree<DATATYPE, ELEMTYPE, NUMDIMS, ELEMTYPEREAL, TMAXNODES, TMINNODES>

RTREE_TEMPLATE
std::vector<WeightPoint> MWK(RTREE_QUAL P, DataPoint q, int k, std::vector<WeightPoint> Wm, int S);
bool dominate(DataPoint q1, DataPoint q2);
bool lesserSort(DataPoint q1, DataPoint q2);
double Peanlty(std::vector<WeightPoint> Wapos, int kapos, std::vector<WeightPoint> Wm, int k, double alpha, double beta);


/** input(1)dataset P                           R-tree
 *       (2)query point q                       std::vector<float>
 *       (3)k                                   int
 *       (4)why-not weighting std::vector set Wm     std::vector set
 *       (5)sample size |S|                     int
 * output(1)modified W'm                        std::vector set
 *       (2)modified k'                         int
 */
RTREE_TEMPLATE
std::tuple<std::vector<WeightPoint>, int> MWK(RTREE_QUAL P, DataPoint q, int k, std::vector<WeightPoint> Wm, int sampleSize, double alpha, double beta){
        //D, I 分别存储
    std::vector<DataPoint> D, I;
    //FindIncom 找到与q incomparable和dominate q的节点
    //首先判断RTree中仍有数据
    typename RTREE_QUAL::Iterator it;
    for(P.GetFirst(it); it.IsNotNull(); ++it)
    {
        ELEMTYPE boundsMin[NUMDIMS];
        ELEMTYPE boundsMax[NUMDIMS];
        it.GetBounds(boundsMin, boundsMax);
        ELEMTYPE e[NUMDIMS];
        for(int i = 0; i < NUMDIMS; i++)
        {
            e[i] = (boundsMin[i] + boundsMax[i]) / 2;
        }
        DataPoint p(e, NUMDIMS, 0);
        //如果p dominate q
        if(dominate(p, q))
        {
            D.push_back(p);
        }
        //如果q不dominate p
        else if(!dominate(q, p))
        {
            I.push_back(p);
        }
    }
    
    int kmax = 0;
    std::vector<WeightPoint> S(sampleSize);
    std::vector<int> Srank(sampleSize);
    ///sample |S| weighting std::vector from the hyperplanes
    //即在I中所有点与q的连线和∑w=1的交面中随机抽取点放在S中
    //抽取向量总数为sample个
    for(int i = 0; i < sampleSize; i++)
    {
        //首先随机选择一个超平面
        srand(time(0));
        WeightPoint hyperplane0 = I[rand()%I.size()];///此处有问题，rand()产生的随机数上界仅为32767可能远远小于数据集长度
        std::vector<double> cof(hyperplane0.getDim());
        
        
        WeightPoint w;
        //再计算随机选择的超平面与∑w=1相交而成的超直线
        //总共d维
        for(int j = 0; j < q.getDim(); j++)
        {
        
        }
        for(int k = 0; k < q.getDim(); k++)
        {
            if(w[k] < 0)
            {
                printf("there is a 0 elem at w[%d]", k);
                exit(-1);
            }
        }
        S[i] = w;
    }
    
    
    
    for(auto s = S.begin(); s != S.end(); s++)
    {
        //q最高排名|D|+1
        int rank = D.size() + 1;
        q.setScore(*s);
        //compute the ranking of q based on D and I
        //q的排名取决于I中的数据
        for(auto it  = I.begin(); it != I.end(); it++)
        {
            it->setScore(*s);
            if(it->getScore() < q.getScore())
            {
                rank += 1;
            }
        }
        Srank[s - S.begin()] = rank;
    }
    
    //对S中的抽样点按照排名进行排序
    sort(S.begin(), S.end(), lesserSort);
    
    //计算Wm中每个用户对q的排位，并设置最大排位
    for(auto w = Wm.begin(); w != Wm.end(); w++)
    {
        //compute the rank of q in every w in Wm
        //q最高排名|D|+1
        int rank = D.size() + 1;
        q.setScore(*w);
        //compute the ranking of q based on D and I
        //q的排名取决于I中的数据
        for(auto it  = I.begin(); it != I.end(); it++)
        {
            it->setScore(*w);
            if(it->getScore() < q.getScore())
            {
                rank += 1;
            }
        }
        if(kmax < rank)
        {
            kmax = rank;
        }
    }
    
    //初始化CW为S的第一个元素
    std::vector<WeightPoint> CW(Wm.size(), S[0]);
    std::vector<WeightPoint> Wapos = Wm;
    int kapos = kmax;
    double Pmin = Peanlty(Wapos, kapos, Wm, kmax, alpha, beta);
    
    for(auto s = S.begin()+1; s != S.end(); s++)
    {
        //q在每个
        int rs = Srank[s - S.begin()];
        if(kmax < rs)
        {
            break;
        }
        
        bool CWUpdaated = false;
        auto cw = CW.begin();
        auto w = Wm.begin();
        while(cw != CW.end() && w != Wm.end())
        {
            //计算cw 和 s 对于所有 w 的差值并替换cw
            if(s->distance(*w) < cw->distance(*w))
            {
                *cw = *s;
                CWUpdaated = true;
            }
            cw++;
            w++;
        }
        if(CWUpdaated)
        {
            if(Peanlty(CW, rs, Wm, kmax, alpha, beta) < Pmin)
            {
                Wapos = CW;
                kapos = k>rs?k:rs;
                Pmin = Peanlty(Wapos, rs, Wm, kmax, alpha, beta);
            }
        }
    }
    return std::make_tuple(Wm, kapos);
}


bool dominate(DataPoint q1, DataPoint q2){
    float sum = 0;
    //如果维度不同，无法比较
    if(q1.size() != q2.size()){
        std::cout<<"q1, q2 has different dimensions, can not be compared"<<std::endl;
        return false;
    }
    for(int i = 0; i < q1.size(); i++){
        //如果存在一个维度q1比q2大，则不成立
        if(q1[i] > q2[i]){
            return false;
        }else{
            sum += q1[i] - q2[i];
        }
    }
    //q1存在一个维度比q2小则成立
    return sum < 0;
}

bool lesserSort(DataPoint q1, DataPoint q2)
{
    return q1.getRank()<q2.getRank();
}

double Peanlty(std::vector<WeightPoint> Wapos, int kapos, std::vector<WeightPoint> Wm, int k, double alpha, double beta)
{
    double peank, peanWm;
    //计算k
    double molecule_k = kapos-k>0?kapos-k:0;
    double denominator_k = Wapos[0].getRank();
    for(auto w = Wapos.begin(); w != Wapos.end(); w++)
    {
        if(denominator_k < w->getRank())
        {
            denominator_k = w->getRank();
        }
    }
    denominator_k -= k;
    peank = molecule_k/denominator_k;
    //计算Wm
    double molecule_Wm = 0;
    double denominator_Wm = 0;
    auto w = Wm.begin();
    auto wapos = Wapos.begin();
    while(w != Wm.end() && wapos != Wapos.end())
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
    
    return alpha*peank + beta*peanWm;
}

#endif