#include <stdio.h>
#include <string>
#include <iostream>
#include <random>
#include <math.h>
#include <algorithm>
#include <vector>
#include <list>
#include "temporal.h"
#include<time.h>


using namespace std;


int i_parent(int i) {return (i-1) / 2;}
int i_left_child(int i) {return 2 * i + 1;}
int i_right_child(int i) {return 2 * i + 2;}

void swap_element(vector<vector<int16_t>*>* a, int b, int c) 
{
	vector<int16_t>* temp = a->at(b);
	a->at(b) = a->at(c);
	a->at(c) = temp;
}

void sift_down(vector<vector<short>*> *a, int start, int end) 
{
	int root = start;

	while(i_left_child(root) <= end) {
		int child = i_left_child(root);
		int swap = root;

		if(a->at(swap)->front() < a->at(child)->front()) {
			swap = child;
		}

		if(child + 1 <= end && a->at(swap)->front() < a->at(child + 1)->front()) {
			swap = child + 1;
		}

		if(swap == root) {
			return;
		} else {
			swap_element(a, root, swap);
			root = swap;
		}
	}
}

void heapify(vector<vector<int16_t>*>* a) 
{
	int count = a->size();
	int start = i_parent(a->size() - 1);

	while(start >= 0) {
		sift_down(a, start, count - 1);
		start -= 1;
	}
}

void sort_neighbors(vector<vector<short>*>* a) 
{
	heapify(a);
	int end = a->size() - 1;

	while(end > 0) {
		swap_element(a, end, 0);
		end -= 1;
		sift_down(a, 0, end);
	}

    a->erase(a->begin() + 10, a->end());
}

//void push_in_heap(vector<vector<short>*>* my_vector, vector<short>* element)
//{
//    vector<vector<short>*>::iterator it = my_vector->begin();
//    if(my_vector->size() == 0)
//        my_vector->insert(my_vector->end(), element);
//    else if(my_vector->size() == 1 && my_vector->at(0)->at(2) > element->at(2))
//        my_vector->insert(it, element);
//    else if(my_vector->size() == 1 && my_vector->at(0)->at(2) < element->at(2))
//        my_vector->insert(my_vector->end(), element);
//    else {
//        for(uint i = 0; i < my_vector->size() - 1; ++i) 
//            if(my_vector->size() > 1 && 
//                    (*(it + i + 1))->at(2) >= element->at(2)) {    
//                my_vector->insert(it + i, element);
//                maintain_size(my_vector);
//                return;
//            }
//    } 
//}
//
//void maintain_size(vector<vector<short>*>* my_vector)
//{
//    //cout<<"coming here --------------------------------------------------------------------"<<endl;
//    if(my_vector->size() > K) {
//       my_vector->erase(my_vector->begin() + my_vector->size() - 1);
//    }    
//}

void print_heap(vector<vector<int16_t>> heap, int x, int y) 
{
	cout<<"Heap"<<endl;
	cout<<"y "<<y<<" x "<<x<<endl;
	for(uint i = 0; i < heap.size(); i++) {
		print_v_i(heap.at(i));

		cout<<endl;
	}

	cout<<endl;
}

void print_v_i(vector<short> v_i) 
{
	cout<<"Printing v_i"<<endl;
	for(int i = 0; i < 3; i++) {
		switch(i) {

				case 0:
					cout<<"x_i: ";
					break;

				case 1:
					cout<<"y_i: ";
					break;

				case 2:
					cout<<"D(P(q), P(q_i)): ";
					break;
			
			}

		cout<<v_i.at(i)<<endl;
	}
}

void print_coord(short* coord) 
{
	cout<<"print_coord()"<<endl;
	cout<<"x "<<coord[0]<<endl;
	cout<<"y "<<coord[1]<<endl;
}
