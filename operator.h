//
// Created by Guangnian_DONG on 2018/10/21 021.
//
#ifndef WQRTQ_OPERATOR_H
#define WQRTQ_OPERATOR_H


#include <vector>
#include <iostream>
#include <cstdlib>
#include <cmath>

#define RTREE_TEMPLATE template<class DATATYPE, class ELEMTYPE, int NUMDIMS, class ELEMTYPEREAL, int TMAXNODES, int TMINNODES>
#define RTREE_QUAL RTree<DATATYPE, ELEMTYPE, NUMDIMS, ELEMTYPEREAL, TMAXNODES, TMINNODES>

#define MINHEAP_TEMPLATE template<class DATATYPE, class ELEMTYPE>
#define MINHEAP_QUAL MinHeap<DATATYPE, ELEMTYPE>


class DataPoint
{
private:
    // 存储数据点
    double* __data;
    int __index, __dim, __rank;
    double __score;

public:

    DataPoint(){};

    DataPoint(const double* data, int d, int idx)
    {
        setData(data, d);
        setRank(-1);
        __index = idx;
        __score = 0;
    }

    ~DataPoint()
    {
//        std::cout<<"delete datapoint "<<__index<<std::endl;
        delete[] __data;
        __data = NULL;
    }
    
    bool operator> (const DataPoint& p) const
    {
        return this->__score > p.__score;
    }
    
    bool operator< (const DataPoint& p) const
    {
        return this->__score < p.__score;
    }
    
    double operator[] (int i) const
    {
        if(i > __dim-1)
        {
            std::cout<<"下表范围0-"<<__dim<<"下标越界"<<std::endl;
            exit(-1);
        }
        return __data[i];
    }

    void setScore(const DataPoint w);
    double getScore();
    void setData(const double* data, const int d);
    double* getData();
    void setRank(const int rank);
    double distance(DataPoint p);
    void toString();
    int size();
    const double* rect_min();
    const double* rect_max();
    int getIndex();
    int getRank();
    int getDim();
};

void DataPoint::setScore(const DataPoint w)
{
    double sum = 0;
    for(int i = 0; i < __dim; i++)
    {
        sum += __data[i] * w[i];
    }
    __score = sum;
}

void DataPoint::setData(const double* data, const int d)
{
    __dim = d;
    __data = new double[__dim];
    for(int i = 0; i < __dim; i++)
    {
        __data[i] = data[i];
    }
}

void DataPoint::setRank(const int rank = -1)
{
    __rank = rank;
}

double DataPoint::distance(DataPoint p)
{
    double sum = 0;
    for(int i = 0; i < __dim; i++)
    {
        sum += pow(this->__data[i] - p.getData()[i], 2);
    }
    return sqrt(sum);
}

void DataPoint::toString()
{
    std::cout<<"index: "<<__index<<std::endl;
    for(int i = 0; i < __dim; i++)
    {
        std::cout<<"dim"<<i<<": "<<__data[i]<<" ";
    }
    std::cout<<std::endl;
    std::cout<<"score: "<<__score<<std::endl;
}

double* DataPoint::getData()
{
    return __data;
}

int DataPoint::size()
{
    return __dim;
}

const double* DataPoint::rect_min()
{
    return __data;
}

const double* DataPoint::rect_max()
{
    return __data;
}

int DataPoint::getIndex()
{
    return __index;
}

double DataPoint::getScore()
{
    return __score;
}

int DataPoint::getRank()
{
    return __rank;
}

int DataPoint::getDim()
{
    return __dim;
}


typedef DataPoint WeightPoint;


#endif //WQRTQ_OPERATOR_H
