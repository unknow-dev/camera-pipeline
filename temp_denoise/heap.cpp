//
// Created by Prikshet on 8/16/2019.
//

#include <stdio.h>
#include <string>
#include <iostream>
#include <random>
#include <math.h>
#include <algorithm>
#include <vector>
#include <list>


int i_parent(int i) {return (i-1) / 2;}
int i_left_child(int i) {return 2 * i + 1;}
int i_right_child(int i) {return 2 * i + 2;}

void swap_element(std::vector<std::vector<double> > &v, int index1, int index2) {
  std::vector<double> temp = v[index1];
  v[index1] = v[index2];
  v[index2] = temp;
}

void sift_down(std::vector<std::vector<double> > v, int start, int end)
{
  int root = start;

  while(i_left_child(root) <= end) {
    int child = i_left_child(root);
    int swap = root;

    if(v[swap].front() < v[child].front()) {
      swap = child;
    }

    if(child + 1 <= end && v[swap].front() < v[child + 1].front()) {
      swap = child + 1;
    }

    if(swap == root) {
      return;
    } else {
      swap_element(v, root, swap);
      root = swap;
    }
  }
}

void heapify(std::vector<std::vector<double> > v)
{
  int size = v.size();
  int start = i_parent(size - 1);

  while(start >= 0) {
    sift_down(v, start, size - 1);
    start -= 1;
  }
}

void sort_neighbors(std::vector<std::vector<double> > v)
{
  heapify(v);
  int end = v.size() - 1;

  while(end > 0) {
    swap_element(v, end, 0);
    end -= 1;
    sift_down(v, 0, end);
  }

  v.erase(v.begin() + 10, v.end());
}

//void push_in_heap(std::vector<std::vector<short>*>* my_std::vector, std::vector<short>* element)
//{
//    std::vector<std::vector<short>*>::iterator it = my_std::vector->begin();
//    if(my_std::vector->size() == 0)
//        my_std::vector->insert(my_std::vector->end(), element);
//    else if(my_std::vector->size() == 1 && my_std::vector->at(0)->at(2) > element->at(2))
//        my_std::vector->insert(it, element);
//    else if(my_std::vector->size() == 1 && my_std::vector->at(0)->at(2) < element->at(2))
//        my_std::vector->insert(my_std::vector->end(), element);
//    else {
//        for(uint i = 0; i < my_std::vector->size() - 1; ++i)
//            if(my_std::vector->size() > 1 &&
//                    (*(it + i + 1))->at(2) >= element->at(2)) {
//                my_std::vector->insert(it + i, element);
//                maintain_size(my_std::vector);
//                return;
//            }
//    }
//}
//
//void maintain_size(std::vector<std::vector<short>*>* my_std::vector)
//{
//    //cout<<"coming here --------------------------------------------------------------------"<<endl;
//    if(my_std::vector->size() > K) {
//       my_std::vector->erase(my_std::vector->begin() + my_std::vector->size() - 1);
//    }
//}
void print_v_i(std::vector<double> v_i);

void print_heap(std::vector<std::vector<double> > heap)
{
  std::cout<<"Heap"<<std::endl;
  
  for(uint i = 0; i < heap.size(); i++) {
    print_v_i(heap[i]);

    std::cout<<std::endl;
  }

  std::cout<<std::endl;
}

void print_v_i(std::vector<double> v_i)
{
  std::cout<<"Printing v_i"<<std::endl;
  for(int i = 0; i < 3; i++) {
    switch(i) {

    case 0:
      std::cout<<"x_i: ";
      break;

    case 1:
      std::cout<<"y_i: ";
      break;

    case 2:
      std::cout<<"D(P(q), P(q_i)): ";
      break;

    }

    std::cout<<v_i.at(i)<<std::endl;
  }
}

void print_coord(double* coord)
{
  std::cout<<"print_coord()"<<std::endl;
  std::cout<<"x "<<coord[0]<<std::endl;
  std::cout<<"y "<<coord[1]<<std::endl;
}


int main() {
  printf("Testing heap");

  std::vector<std::vector<double> > v;
  
  std::vector<double> v1;
  v1.push_back(-1); v1.push_back(-1); v1.push_back(-10);
  std::vector<double> v2;
  v2.push_back(4); v2.push_back(-43); v2.push_back(-2);
  std::vector<double> v3;
  v3.push_back(1); v3.push_back(23); v3.push_back(2);
  std::vector<double> v4;
  v4.push_back(0); v4.push_back(-3); v4.push_back(-133);
  std::vector<double> v5;
  v5.push_back(1); v5.push_back(4); v5.push_back(7);
  std::vector<double> v6;
  v6.push_back(-13); v6.push_back(-11); v6.push_back(-21); 

  v.push_back(v1); v.push_back(v2); v.push_back(v3); v.push_back(v4); v.push_back(v5); v.push_back(v6);

  heapify(v);
  print_heap(v);
}
