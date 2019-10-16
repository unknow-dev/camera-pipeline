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
#include "Halide.h"
#include "halide_image_io.h"
#include "temporal.h"
#include <unordered_map>

using namespace std;
using namespace Halide;
using namespace Halide::ConciseCasts;


#include "x_l_in.h"
#include "x_r_in.h"
#include "y_u_in.h"
#include "y_d_in.h"
#include "x_l_out.h"
#include "x_r_out.h"
#include "y_u_out.h"
#include "y_d_out.h"

const int n_channels = 3;

void get_frames(Halide::Runtime::Buffer<uint8_t>* &frame_b, string frames_dir)
{
    string path
    for(int i = 0; i < n_frames; i++) {
        path =  frames_dir + to_string(i + 1) + ".png";
        frame_b[i] = Tools::load_image(path);
    }
}

int main() {
    int n_frames;
    Halide::Runtime::Buffer<uint8_t> frame_b[n_frames];
    string frames_dir = "./frames/";
    get_frames(frame_b, frames_dir);

}




void generate_random_offsets_and_ssds(short current_frame, vector<short>* point, vector<vector<short>*>* neighbors)
{
for(ushort i = 0; i < K; i++) {
vector<short>* ssd_and_offset = get_neighbor_ssd(current_frame, point);
neighbors->push_back(ssd_and_offset);
}
sort_neighbors(neighbors);
//if(x == 33 && y == 33)
//    cout<<"neighbors of 33, 33 when generating randoms"<<endl;
//    print_neighbors(neighbors);
}

void initiate_neighbors(short current_frame, vector<vector<short>*>* neigh_heap[height][width])
{

srand(17);

for(short y = 0; y < height; y++)
for(short x = 0; x < width; x++) {

    vector<vector<short>*>* neighbors;
    neighbors = new vector<vector<short>*>;
    vector<short>* point = new vector<short>;

    point->push_back(x);
    point->push_back(y);

    generate_random_offsets_and_ssds(current_frame, point, neighbors); // 0 seconds

    neigh_heap[y][x] = neighbors;

} // 2.42 seconds


//cout<<"initial heap size "<<neighbors_h[0][0]->size()<<endl ;
}



vector<short>* get_neighbor_ssd(short current_frame, vector<short>* point)
{
    vector<short>* coord = new vector<short>;
    Buffer<short> pix(1, 1, 1, 1, 3);
    short offset_x = get_random_x();
    short offset_y = get_random_y();
    short x = point->at(0);
    short y = point->at(1);
    pix.set_min(x, y, x + offset_x, y + offset_y);
    D[current_frame].realize(pix);
    short ssd = pix(x, y, x + offset_x, y + offset_y, 0);
    coord->push_back(ssd);
    //cout<<"ssd "<<ssd;
    coord->push_back(offset_x);
    coord->push_back(offset_y);
    return coord;
}



void print_offset_and_ssd(vector<short>* offset_ssd)
{
    cout<<"ssd "<<offset_ssd->at(0)<<" x_i "<<offset_ssd->at(1)
        <<" y_i "<<offset_ssd->at(2)<<endl;
}

void print_neighbors(vector<vector<short>*>* heap)
{
for(uint i = 0; i < heap->size(); i++) {
print_offset_and_ssd(heap->at(i));
}
}

//vector<short>* get_random_coord()
//{
//    vector<short>* coord = new vector<short>;
//    coord->push_back(get_random_x());
//    coord->push_back(get_random_y());
//    return coord;
//}

short get_random_x()
{
    float random_x;
    for(;;) {
        random_x = width / 3 * box_muller_trans((float) rand() / RAND_MAX);
        if(random_x < width)
            return (short) random_x;
    }
}

short get_random_y()
{
    float random_y;
    for(;;) {
        random_y = height / 3 * box_muller_trans((float) rand() / RAND_MAX);
        if(random_y < height)
            return (short) random_y;
    }
}

// converts a uniform random variable into a standard normal variable
float box_muller_trans(float x)
{
    return sqrt(-2 * log(x)) * cos(2 * M_PI * x);;
}




void propagate_neighbors(short frame, vector<vector<short>*>* neigh_heap[height][width])
{

for(short y = 1; y < height; ++y)
for(short x = 1; x < width; ++x) {

propagate_scanline(frame, x, y, neigh_heap);

} // 1.6-4.1 seconds

for(short y = height - 2; y >= 0; --y)
for(short x = width - 2; x >= 0; --x) {
propagate_reverse_scanline(frame, x, y, neigh_heap);
}
}


void propagate_scanline(short frame, short x, short y, vector<vector<short>*>* neigh_heap[height][width])
{
short offset_x, offset_y, offset_ssd;
clock_t t = clock();


for(vector<vector<short>*>::iterator it = neigh_heap[y][x - 1]->begin();
it != neigh_heap[y][x - 1]->end(); ++it) {
vector<short>* new_neighbor = new vector<short>;
offset_ssd = (*it)->at(0);
offset_x = (*it)->at(1);
offset_y = (*it)->at(2);
new_neighbor->push_back(calculate_new_ssd(frame, x, y, offset_x, offset_y, offset_ssd, 'r'));
new_neighbor->push_back((short) offset_x + 1);
new_neighbor->push_back(offset_y);
neigh_heap[y][x]->push_back(new_neighbor);
}
t = clock() - t;
cout<<"time tester "<<(float)t/CLOCKS_PER_SEC<<" seconds"<<endl;

for(vector<vector<short>*>::iterator it = neigh_heap[y - 1][x]->begin();
it != neigh_heap[y - 1][x]->end(); ++it) {
vector<short>* new_neighbor = new vector<short>;
offset_ssd = (*it)->at(0);
offset_x = (*it)->at(1);
offset_y = (*it)->at(2);
new_neighbor->push_back(calculate_new_ssd(frame, x, y, offset_x, offset_y, offset_ssd, 'd'));
new_neighbor->push_back(offset_x);
new_neighbor->push_back((short) offset_y + 1);
neigh_heap[y][x]->push_back(new_neighbor);
}
sort_neighbors(neigh_heap[y][x]);
}

void propagate_reverse_scanline(short frame, short x, short y, vector<vector<short>*>* neigh_heap[height][width])
{
short offset_x, offset_y, offset_ssd;
for(vector<vector<short>*>::iterator it = neigh_heap[y][x + 1]->begin();
it != neigh_heap[y][x + 1]->end(); ++it) {
vector<short>* new_neighbor = new vector<short>;
offset_ssd = (*it)->at(0);
offset_x = (*it)->at(1);
offset_y = (*it)->at(2);
new_neighbor->push_back(calculate_new_ssd(frame, x, y, offset_x, offset_y, offset_ssd, 'l'));
new_neighbor->push_back((short) offset_x - 1);
new_neighbor->push_back(offset_y);
}

for(vector<vector<short>*>::iterator it = neigh_heap[y + 1][x]->begin();
it != neigh_heap[y + 1][x]->end(); ++it) {
vector<short>* new_neighbor = new vector<short>;
offset_ssd = (*it)->at(0);
offset_x = (*it)->at(1);
offset_y = (*it)->at(2);
new_neighbor->push_back(calculate_new_ssd(frame, x, y, offset_x, offset_y, offset_ssd, 'u'));
new_neighbor->push_back(offset_x);
new_neighbor->push_back((short)offset_y - 1);
}
sort_neighbors(neigh_heap[y][x]);
}

short calculate_new_ssd(short frame, short x, short y, short offset_x,
                        short offset_y, short offset_ssd, char direction)
{
    Buffer<short> add_b(1, 1, 1, 1, 3);
    Buffer<short> subtract_b(1, 1, 1, 1, 3);

    add_b.set_min(x, y, x + offset_x, y + offset_y, 0);
    subtract_b.set_min(x, y, x + offset_x, y + offset_y, 0);

    switch(direction) {

        case 'r':
            Dx_right_out[frame].realize(add_b);
            Dx_left_in[frame].realize(subtract_b);
            break;

        case 'd':
            Dy_down_out[frame].realize(add_b);
            Dy_up_in[frame].realize(subtract_b);
            break;

        case 'l':
            Dx_left_out[frame].realize(add_b);
            Dx_right_in[frame].realize(subtract_b);
            break;

        case 'u':
            Dy_up_out[frame].realize(add_b);
            Dy_down_in[frame].realize(subtract_b);
            break;

        default:
            cerr<<"wrong direction character"<<endl;
    }
    return offset_ssd - subtract_b(x, y, x + offset_x, y + offset_y, 0)
           + add_b(x, y, x + offset_x, y + offset_y, 0);
}


void random_search(short frame, vector<vector<short>*>* neighbors_h[height][width])
{
int M = min(log(width / 3), (double)K);
//cout<<"M "<<M<<endl;
for(int y = 0; y < height; ++y) {
for(int x = 0; x < width; ++x) {
for(int i = 0; i < M; ++i) {
    vector<short>* random_guess = new vector<short>;

    short offset_x = get_random_x() * pow(0.5, i);
    short offset_y = get_random_y() * pow(0.5, i);

    Buffer<short> pix(1, 1, 1, 1, 3);

    pix.set_min(x, y, x + offset_x, y + offset_y, 0);
    D[frame].realize(pix);

    short ssd = pix(x, y, x + offset_x, y + offset_y, 0);

    random_guess->push_back(ssd);
    random_guess->push_back(offset_x);
    random_guess->push_back(offset_y);

    //cout<<"ssd "<<ssd<<" offset_x "<<offset_x<<" offset_y "<<offset_y<<endl;

    neighbors_h[y][x]->push_back(random_guess);
    }

    sort_neighbors(neighbors_h[y][x]);
    //for(vector<vector<short>*>::iterator it = neighbors_h[y][x]->begin();
    //    it != neighbors_h[y][x]->end(); ++it) {
    //    cout<<"ssd "<<(*it)->at(0)<<" offset_x "<<(*it)->at(1)<<" offset_y "
    //        <<(*it)->at(2)<<endl;
    //}
}
}
}

vector<vector<short>*>* get_neighbors(short x, short y, vector<vector<short>*>* neighbors_h[height][width])
{
return neighbors_h[x][y];
}


float calc_weighted_ssd(short current_frame, short other_frame, vector<short>* patch_coord_current_frame,
                        vector<short>* patch_coord_other_frame, ushort neighbor)
{
    load_halide_functions_nlm(current_frame, other_frame,
                              patch_coord_current_frame, patch_coord_other_frame);

    Buffer<int16_t> weighted_ssd_buff(1, 1, 1, 1, 3);

    short x = patch_coord_current_frame->at(0);
    short y = patch_coord_current_frame->at(1);

    short x_ = patch_coord_other_frame->at(0);
    short y_ = patch_coord_other_frame->at(1);


    weighted_ssd_buff.set_min(x, y, x_, y_, 0);
    weighted_ssd[current_frame][other_frame + H][neighbor].realize(weighted_ssd_buff);
    return weighted_ssd_buff(x, y, x_, y_, 0);
}

short weighted_ssd_neighbor_counter = 0;
void load_halide_functions_nlm(short current_frame, short other_frame, vector<short>* patch_coord_current_frame, vector<short>* patch_coord_other_frame)
{
    Var x, y, x_i, y_i, c;
    RDom u(-s, s, -s, s);
    for(ushort neighbor = 0; neighbor < K; neighbor++) {
        weighted_ssd[current_frame][other_frame + H][neighbor](x, y, x_i, y_i, c) = i16(sum(pow((I[current_frame](x + u.x, y + u.y, c)
                                                                                                 - I[current_frame + other_frame](x_i + u.x, y_i + u.y, c)), 2)
                                                                                            * exp(-(pow(u.x, 2) + pow(u.y, 2)) /
                                                                                                  (float)(2 * pow(s / 2, 2)))));
    }
}