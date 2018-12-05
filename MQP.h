//
// Created by Guangnian_DONG on 2018/10/14 014.
//
#ifndef WQRTQ_MQP_H
#define WQRTQ_MQP_H

#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include "operator.h"
#include "minheap.h"
#include "QuadProg/QuadProg++.hh"
#include "RTree/RTree.h"

#define RTREE_TEMPLATE template<class DATATYPE, class ELEMTYPE, int NUMDIMS, class ELEMTYPEREAL, int TMAXNODES, int TMINNODES>
#define RTREE_QUAL RTree<DATATYPE, ELEMTYPE, NUMDIMS, ELEMTYPEREAL, TMAXNODES, TMINNODES>

RTREE_TEMPLATE
DataPoint MQP(RTREE_QUAL P, DataPoint q, int k, std::vector<WeightPoint> Wm);

/**input(1)dataset P                         R-Tree
 *      (2)query point q                     vector<double>
 *      (3)parameter k                       int
 *      (4)why-not weigting vector set Wm    vector set
 *output(1)q'                                vector
 */
RTREE_TEMPLATE
DataPoint MQP(RTREE_QUAL P, DataPoint q, int k, std::vector<WeightPoint> Wm)
{
    //first find every k-th point for each elem in Wm
    //最小堆存储RTree的节点迭代器，存取rectangle的话无法判断是否是叶子节点
    //最小堆, 用于计算每个missing vector的k-th point
    MinHeap<DataPoint, double> HP(20);
    std::vector<DataPoint > lambda(Wm.size()); //存储Wm集合中对应的排名k的数据点
    //首先层次遍历R树，计算每个missing vctor 的k-th 向量
    for(auto w= Wm.begin(); w != Wm.end(); w++)
    {
        ///该部分构造最小堆方法有问题，当数据集非常大的时候，需要遍历所有的数据，带来极大的计算开支
        // 按照datavector在w下的得分构造最小堆
        typename RTREE_QUAL::Iterator it;
        for(P.GetFirst(it); it.IsNotNull(); ++it)
        {
            //R树中的数据和Datapoint做转换
            ELEMTYPE boundsMin[NUMDIMS];
            ELEMTYPE boundsMax[NUMDIMS];
            it.GetBounds(boundsMin, boundsMax);
            ELEMTYPE e[NUMDIMS];
            for(int i = 0; i < NUMDIMS; i++)
            {
                e[i] = (boundsMin[i] + boundsMax[i]) / 2;
            }
            DataPoint p(e, NUMDIMS, 100);
            p.setScore(*w);
            HP.insert(p);
        }
        int count = 0;
        while(HP.size() > 0)
        {
            DataPoint e = HP.deheap();
            count++;
            if(count == k)
            {
                lambda[w-Wm.begin()] = e;
                break;
            }
        }
    }

    WeightPoint q_m; //修改后的q
    //约束条件Ax<=b
    //目标方程运算量
    quadprogpp::Matrix<double> G, CE, CI;
    //目标函数c
    //约束条件b
    //约束条件lb<=x<=ub
    quadprogpp::Vector<double> g0, ce0, ci0, x;
    //二次规划的变量个数
    int d, m, n;
    d = q.size();   //变量个数
    m = 0;          //约束中等式约束的个数
    n = Wm.size();  //约束中不等式约束的个数
    //而后根据每个missing vector 的k-th 向量形成二次规划的参数
    //初始化H = diag(2,2,...,2), c = (-2q[1],-2q[2],...,-2q[d])
    G.resize(d, d);
    g0.resize(d);
    for(int i = 0; i < d; i++)
    {
        for(int j = 0; j < d; j++)
        {
            if(i == j)  {G[i][j] = 2;}
            else        {G[i][j] = 0;}
        }
        g0[i] = -2*q[i];
    }
    //CI^T x + ci0 >= 0
    //CI为(n+2d)*n的矩阵，CI^T为n*(n+2d)的矩阵
    CI.resize(d, n+2*d);
    ci0.resize(n+2*d);
    for(int j = 0; j < n+2*d; j++)
    {
        //添加top-k约束
        if(j < n)
        {
            double sum = 0;
            for(int i = 0; i < d; i++)
            {
                //-Ax+b>=0
                CI[i][j] = (-1) * Wm[j][i];
                sum += Wm[j][i] * lambda[j][i];
            }
            ci0[j] = sum;
        }
        //添加下界约束
        else if(j >= n && j < n+d)
        {
            int idx = j - n;
            for(int i = 0; i < d; i++)
            {
                if(idx == i)    {CI[i][j] = 1;}
                else            {CI[i][j] = 0;}
            }
            ci0[j] = 0;
        }
        //添加上界约束
        else if(j >= n+d && j<n+2*d)
        {
            int idx = j - n - d;
            for(int i = 0; i < d; i++)
            {
                if(idx == i)    {CI[i][j] = -1;}
                else            {CI[i][j] = 0;}
            }
            ci0[j] = q[idx];
        }
    }
    CE.resize(d,0);
    ce0.resize(0);
    x.resize(d);
    
    std::cout<<std::endl;
    std::cout<<"Wm = "<<std::endl;
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < d; j++)
        {
            std::cout<<Wm[i][j]<<" ";
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
    std::cout<<"lambda = "<<std::endl;
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < d; j++)
        {
            std::cout<<lambda[i][j]<<" ";
        }
        std::cout<<std::endl;
    }

    std::cout<<std::endl;
    std::cout<<"G = "<<std::endl;
    for(int i = 0; i < d; i++)
    {
        for(int j = 0; j < d; j++)
        {
            std::cout<<G[i][j]<<" ";
        }
        std::cout<<std::endl;
    }

    std::cout<<std::endl;
    std::cout<<"g0 = "<<std::endl;
    for(int i = 0; i < d; i++)
    {
        std::cout<<g0[i]<<" ";
    }

    std::cout<<std::endl;
    std::cout<<"CI = "<<std::endl;
    for(int i = 0; i < d; i++)
    {
        for(int j = 0; j < n+2*d; j++)
        {
            std::cout<<CI[i][j]<<" ";
        }
        std::cout<<std::endl;
    }

    std::cout<<std::endl;
    std::cout<<"ci0 = "<<std::endl;
    for(int i = 0; i < n+2*d; i++)
    {
        std::cout<<ci0[i]<<" ";
    }
    std::cout<<std::endl;


    //QuadProg 通过二次规划解决问题
    double cost = quadprogpp::solve_quadprog(G, g0, CE, ce0, CI, ci0, x);
    double* qm = new double[x.size()];
    std::cout<<std::endl;
    std::cout<<"cost: "<<cost<<std::endl;
    for(int i = 0; i < x.size(); i++)
    {
        qm[i] = x[i];
        std::cout<<i<<": "<<qm[i]<<" ";
    }
    q_m.setData(qm, x.size());
    delete[] qm;
    return q_m;
}

#undef RTREE_TEMPLATE
#undef RTREE_QUAL

#endif