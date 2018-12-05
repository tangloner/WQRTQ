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
    //��С�Ѵ洢RTree�Ľڵ����������ȡrectangle�Ļ��޷��ж��Ƿ���Ҷ�ӽڵ�
    //��С��, ���ڼ���ÿ��missing vector��k-th point
    MinHeap<DataPoint, double> HP(20);
    std::vector<DataPoint > lambda(Wm.size()); //�洢Wm�����ж�Ӧ������k�����ݵ�
    //���Ȳ�α���R��������ÿ��missing vctor ��k-th ����
    for(auto w= Wm.begin(); w != Wm.end(); w++)
    {
        ///�ò��ֹ�����С�ѷ��������⣬�����ݼ��ǳ����ʱ����Ҫ�������е����ݣ���������ļ��㿪֧
        // ����datavector��w�µĵ÷ֹ�����С��
        typename RTREE_QUAL::Iterator it;
        for(P.GetFirst(it); it.IsNotNull(); ++it)
        {
            //R���е����ݺ�Datapoint��ת��
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

    WeightPoint q_m; //�޸ĺ��q
    //Լ������Ax<=b
    //Ŀ�귽��������
    quadprogpp::Matrix<double> G, CE, CI;
    //Ŀ�꺯��c
    //Լ������b
    //Լ������lb<=x<=ub
    quadprogpp::Vector<double> g0, ce0, ci0, x;
    //���ι滮�ı�������
    int d, m, n;
    d = q.size();   //��������
    m = 0;          //Լ���е�ʽԼ���ĸ���
    n = Wm.size();  //Լ���в���ʽԼ���ĸ���
    //�������ÿ��missing vector ��k-th �����γɶ��ι滮�Ĳ���
    //��ʼ��H = diag(2,2,...,2), c = (-2q[1],-2q[2],...,-2q[d])
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
    //CIΪ(n+2d)*n�ľ���CI^TΪn*(n+2d)�ľ���
    CI.resize(d, n+2*d);
    ci0.resize(n+2*d);
    for(int j = 0; j < n+2*d; j++)
    {
        //���top-kԼ��
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
        //����½�Լ��
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
        //����Ͻ�Լ��
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


    //QuadProg ͨ�����ι滮�������
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