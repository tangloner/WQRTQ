//
// Created by Guangnian_DONG on 2018/10/29 029.
//
#ifndef WQRTQ__minheap_H
#define WQRTQ__minheap_H

#include <iostream>

#define MINHEAP_TEMPLATE template<class DATATYPE, class ELEMTYPE>
#define MINHEAP_QUAL MinHeap<DATATYPE, ELEMTYPE>

template<class T>
void swap(T& a, T& b)
{
    T c = a;
    a = b;
    b = c;
}


///通过DATATYPE数据类型的__minheap一维数组来存储最小堆
template<class DATATYPE, class ELEMTYPE>
class MinHeap
{
private:
    DATATYPE* __minheap;
    int __size, __maxSize;

public:
    MinHeap()
    {
        __maxSize = 100;
        __minheap = new DATATYPE[__maxSize];
        __size = -1;
    };

    MinHeap(int maxSize)
    {
        __maxSize = maxSize;
//        __minheap = (DATATYPE*)malloc(__maxSize* sizeof(DATATYPE));
        __minheap = new DATATYPE[__maxSize];
        __size = -1;
    }

    ~MinHeap()
    {
        std::cout<<"正在删除最小堆"<<std::endl;
        delete[] __minheap;
        __minheap = NULL;
    }

    // 返回最小堆
    DATATYPE* getMinHeap();
    void adjustDown(int idx);
    void adjustUp(int idx);
    bool insert(DATATYPE x);
    DATATYPE deheap();

    void toString();
    bool isEmpty();
    bool isFull();
    int size();
};

MINHEAP_TEMPLATE
void MINHEAP_QUAL::adjustDown(int idx)
{
    if(isEmpty() || idx > __size)
    {
        std::cout<<"##################################################"<<std::endl;
        std::cout<<"empty heap or index larger than heap size..."<<std::endl;
        return;
    }
    //有孩子
    while (2 * idx + 1 <= __size)
    {
        int leftChild = 2 * idx + 1;
        int rightChild = 2 * idx + 2;
        //如果左孩子是最后一个节点
        if(leftChild == __size && __minheap[idx] > __minheap[leftChild])
        {
            swap(__minheap[idx], __minheap[leftChild]);
            return;
        }
        //有左孩子也有右孩子
        //左孩子比较大
        if(__minheap[leftChild] > __minheap[rightChild])
        {
            // 如果值比孩子大
            if(__minheap[idx] > __minheap[rightChild])
            {
                swap(__minheap[idx], __minheap[rightChild]);
                idx = rightChild;
            }
            // 如果数值已经小于孩子
            else
            {
                idx = __size;
            }
        }
        //右孩子比较大
        else
        {
            // 如果值比孩子大
            if(__minheap[idx] > __minheap[leftChild])
            {
                swap(__minheap[idx], __minheap[leftChild]);
                idx = leftChild;
            }
                // 如果数值已经小于孩子
            else
            {
                idx = __size;
            }
        }
    }
}


///从idx节点开始向上调整堆
MINHEAP_TEMPLATE
void MINHEAP_QUAL::adjustUp(int idx)
{
    while(idx > -1)
    {
        int father = (idx - 1) / 2;
        //如果存在父亲节点开始调整
        if(father > -1)
        {
            if(__minheap[idx] < __minheap[father])
            {
                swap(__minheap[idx], __minheap[father]);
                idx = father;
            }
            else
            {
                idx = -1;
            }
        }
        //如果不存在父亲节点，调整结束
        else
        {
            idx = -1;
        }
    }
}

MINHEAP_TEMPLATE
bool MINHEAP_QUAL::insert(DATATYPE x)
{
    if(isFull())
    {
        return false;
    }
    __minheap[++__size] = x;
    adjustUp(__size);
    return true;
}

MINHEAP_TEMPLATE
DATATYPE* MINHEAP_QUAL::getMinHeap()
{
    return __minheap;
}

MINHEAP_TEMPLATE
DATATYPE MINHEAP_QUAL::deheap()
{
    if(__size == -1)
    {
        std::cout<<"error: empty heap"<<std::endl;
        exit(-1);
    }
    DATATYPE item = __minheap[0];
    __minheap[0] = __minheap[__size--];
    adjustDown(0);
    return item;
}

MINHEAP_TEMPLATE
bool MINHEAP_QUAL::isEmpty()
{
    return __size == -1;
}

MINHEAP_TEMPLATE
bool MINHEAP_QUAL::isFull()
{
    return __size == __maxSize;
}

MINHEAP_TEMPLATE
void MINHEAP_QUAL::toString()
{
    for(int i = 0; i <=__size; i++)
    {
        __minheap[i].toString();
    }
}

MINHEAP_TEMPLATE
int MINHEAP_QUAL::size()
{
    return this->__size+1;
}

#undef MINHEAP_TEMPLATE
#undef MINHEAP_QUAL

#endif //WQRTQ__minheap_H
