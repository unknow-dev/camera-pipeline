#include<stdio.h>
#include<vector>
#include<iostream>
using namespace std;

int main(int argc, char **argv) {
  vector<int> num;
  for (int i = 0; i <= 100; i++) {
    int a = i;
    num.push_back(a);
  }
  for (int i = 0; i <= 100; i++) {
    cout << num[i];
  }
}
