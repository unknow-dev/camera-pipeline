//
// Created by Prikshet on 8/4/2019.
//
//
// g++ energy*generate.cpp -g -std=c++11 -I ../tools -I ../include -L ../bin -lHalide -lpthread -ldl -o energy_generate && ./energy_generate

// g++ energy*generate.cpp -g -I ../include -I ../tools -L ../bin -lHalide `libpng-config --cflags --ldflags` -ljpeg -lpthread -ldl -o energy_generate -std=c++11 && ./energy_generate
// g++ sobel*generate.cpp -g -std=c++11 -I ../include -L ../bin -lHalide -lpthread -ldl -o sobel_generate && ./sobel_generate
// g++ support*generate.cpp -g -std=c++11 -I ../include -L ../bin -lHalide -lpthread -ldl -o support_generate && ./support_generate

#include<stdio.h>
#include<string>
#include <math.h>
#include "clock.h"
#include <vector>
#include "sobel_x_out.h"
#include "sobel_y_out.h"
#include "support.h"
#include "HalideBuffer.h"
#include "halide_image_io.h"
//#include <fstream>
#include "triangle.h"
#include "energy.h"
#include "Plane.h"
#include "energy_grad_des_gen.h"


#define MAXIMUM_FLOAT
#define MINIMUM_FLOAT

float gamma = 1;
float sigma = 1;
float mu = 1;
float beta = 1;
float tau = 0.9;

vector<vector<Expr>> get_triangle_points(vector<vector<Expr>> support_points, vector<Expr> triangle) {
    vector<vector<Expr>> support_points_triangle;
    Expr x, y, z;
    vector<Expr> point;
    for(int i = 0; i < support_points.size(); i++) {

        x = select(is_inside(support_points[i], triangle), support_points[i][0], -1.f);
        y = select(is_inside(support_points[i], triangle), support_points[i][1], -1.f);
        x = select(is_inside(support_points[i], triangle), support_points[i][2], -1.f);
        point.push_back(x);
        point.push_back(y);
        point.push_back(z);
        support_points_triangle.push_back(point);
//        if(is_inside(support_points[i], triangle)) {
//            support_points_triangle.push_back(support_points[i]);
//        }
    }

    return support_points_triangle;
}

// This function finds the solution of system of
// linear equations using cramer's rule
vector<double> find_solution(double coeff[3][4])
{
    // Matrix d using coeff as given in cramer's rule
    double d[3][3] = {
            { coeff[0][0], coeff[0][1], coeff[0][2] },
            { coeff[1][0], coeff[1][1], coeff[1][2] },
            { coeff[2][0], coeff[2][1], coeff[2][2] },
    };
    // Matrix d1 using coeff as given in cramer's rule
    double d1[3][3] = {
            { coeff[0][3], coeff[0][1], coeff[0][2] },
            { coeff[1][3], coeff[1][1], coeff[1][2] },
            { coeff[2][3], coeff[2][1], coeff[2][2] },
    };
    // Matrix d2 using coeff as given in cramer's rule
    double d2[3][3] = {
            { coeff[0][0], coeff[0][3], coeff[0][2] },
            { coeff[1][0], coeff[1][3], coeff[1][2] },
            { coeff[2][0], coeff[2][3], coeff[2][2] },
    };
    // Matrix d3 using coeff as given in cramer's rule
    double d3[3][3] = {
            { coeff[0][0], coeff[0][1], coeff[0][3] },
            { coeff[1][0], coeff[1][1], coeff[1][3] },
            { coeff[2][0], coeff[2][1], coeff[2][3] },
    };

    // Calculating Determinant of Matrices d, d1, d2, d3
    double D = determinantOfMatrix(d);
    double D1 = determinantOfMatrix(d1);
    double D2 = determinantOfMatrix(d2);
    double D3 = determinantOfMatrix(d3);
    printf("D is : %lf \n", (double)D);
    printf("D1 is : %lf \n", (double)D1);
    printf("D2 is : %lf \n", (double)D2);
    printf("D3 is : %lf \n", (double)D3);

//    vector<double> result(3, 0);
//
//    double x = select(D != 0, D1 / D, (D1 == 0 && D2 == 0 && D3 == 0), -1.f, -2.f);
//    double y = select(D != 0, D2 / D, (D1 == 0 && D2 == 0 && D3 == 0), -1.f, -2.f);
//    double z = select(D != 0, D3 / D, (D1 == 0 && D2 == 0 && D3 == 0), -1.f, -2.f);
//    result[0] = select(D != 0, x, -1.f);
//    result[1] = select(D != 0, y, -1.f);
//    result[2] = select(D != 0, z, -1.f);

    // Case 1
    if (D != 0) {
        // Coeff have a unique solution. Apply Cramer's Rule
        double x = D1 / D;
        double y = D2 / D;
        double z = D3 / D; // calculating z using cramer's rule
        printf("Value of x is : %lf\n", (double)x);
        printf("Value of y is : %lf\n", (double)y);
        printf("Value of z is : %lf\n", (double)z);

        result[0] = x;
        result[1] = y;
        result[2] = z;

    }
//     Case 2
    else {
        if (D1 == 0 && D2 == 0 && D3 == 0)
            printf("Infinite solutions\n");
        else if (D1 != 0 || D2 != 0 || D3 != 0)
            printf("No solutions\n");

    }

    return result;


}

void interpolate(Halide::Runtime::Buffer<double> &a, Halide::Runtime::Buffer<double>  &b,\
Halide::Runtime::Buffer<double>  &c, double** triangulation) {
    for(int i = 0; i < triangulation.size() - 2; i += 3) {
        double x = 0, y = 0, z = 0, xx = 0, xy = 0, yy = 0, xz = 0, yz = 0;

        for (int i = 0; i < triangulation.size(); i++) {
            x += triangulation[i][0];
            y += triangulation[i][1];
            z += triangulation[i][2];
            xx += triangulation[i][0] * triangulation[i][0];
            xy += triangulation[i][0] * triangulation[i][1];
            yy += triangulation[i][1] * triangulation[i][1];
            xz += triangulation[i][0] * triangulation[i][2];
            yz += triangulation[i][1] * triangulation[i][2];
        }
        double coeffs[3][4] = {{xx, xy, x, xz}, {xy, yy, y, yz}, {x, y, (double)triangulation.size(), z}};
        vector<double> result;
        result = find_solution(coeffs);

        a[i] = result[0];
        b[i] = result[1];
        c[i] = result[2];
    }
}

vector<vector<Expr>> convert_to_2d_vector(ImageParam points, Param<int> support_point_size) {
    vector<vector<Expr>> vecs;
    vector<Expr> point;
    for(int p = 0; p < support_point_size - 1; p += 2) {
        for(int i = 0; i < 2; i++) {
            point[0] = points(p);
            point[1] = points(p + 1);
            vecs.push_back(point);
        }
    }
    return vecs;
}

class Observations {
    map<pair, float> observations;

    void add(float x, float y, float d) {
        std::pair p = {x, y};
        observations.insert(p, d);
    }
};

class Observation {
    float u_n;
    float v_n;
    float f_n;

    Observation(u_n, v_n, f_n) {
        this.u_n = u_n;
        this.v_n = v_n
        this.f_n = f_n;
    }


    Observation() {
        u_n = 0;
        v_n = 0;
        f_n = 0;
    }


    Observation f(float d) {
        Observation o = new Observation(l.u_n - d, v_n, observations(l.u_n, v_n));
        return o;
    }
};

double calc_prior(d_n, S, o_n) {

    if(abs(d_n - mu) < 3 * sigma || inN_s(d_n)) {
        return gamma + exp(-pow((d_n - mean_func(S, o_n)), 2) / (2 * sigma * sigma))
    }
    return 0;
}



bool inN_s(point) {
    for(int i = 0; i < support_points.size(); i += 2) {
        if(point[0] == support_points[i] && point[0] == support_points[i + 1]) {
            return true;
        }
    }

    return false;
}


float image_likelihood(Observation r, Observation l, float d_n) {
    if(l.u_n == r.u_n + d_n && l.u_n == r.v_n) {
        return exp(-beta * one_norm(l.f_n - r.f_n));
    }
    return 0;
}


Observations get_obs_points(support_points, triangle) {
    Observations o = new Observations();
    for(int i = 0; i < support_points.size() - 1; i += 2) {
        if(is_inside(triangle[0].x, triangle[0].y, triangle[1].x, triangle[1].y, triangle[2].z, triangle[2].z, support_points[i], support_points[i+1])) {
            o.add(support_points[0], support_points[1], NULL);
        }
    }
}

bool is_best(float curr, float &best, float &second_best) {
    if (curr > best) {
        best = curr;
    }
    if(curr < best && curr > second_best) {
        second_best = curr;
    }
}

int k = 4;

void remove(auto & tuple) {
    for (int i = 0; i < k * 3; i++) {
        tuple[i] = 0.f;
    }
}
void remove_spuriours_matches(Halide::Runtime::Buffer<float> & delaunay[12]) {

    for(int y = 0; y = left_buffer.height(); y++) {
        for(int x = 0; x = left_buffer.width(); x++) {
            float best = MAXIMUM_FLOAT;
            float second_best = MINIMUM_FLOAT;
            for(int i = 0; i < 3 * k; i += 3) {
                sort_best(delaunay[i], best, second_best, vector<int> remove_indices);
            }
            if(best/second_best > tau) {
                remove(delaunay(x, y));
            }
        }
    }


}

double[] to_array(triangleio* support_points) {
    return support_points->pointlist;
}


float gaussian_factor(dist) {
    float term = 1/10;
    return exp(-dist* term);
}

int num_features = 6;
Halide::Runtime::Buffer<float> *generate_features(Halide::Runtime::Buffer<uint8_t> buffer, Halide::Runtime::Buffer<float> sobel) {
    Halide::Runtime::Buffer<float> features[num_features];
    features(x, y)[0] = sobel(x, y)[0];
    features(x, y)[1] = buffer(x, y)[1];
    RDom area(-10, 20, -10, 20);
    Func gaussian;
    gaussian(x, y) = sum(buffer(x - area.x, y - area.y) * gaussian_factor(sqrt(area.x * area.x + area.y * area.y)));

    features(x, y)[2] = gaussian_factor(x, y) / 100;

    Func spatial;
    spatial(x, y) =sum(buffer(x, y) - buffer(x - area.x, y - area.y));

    features(x, y)[3] = spatial(x, y) / 100;

    Func range;
    range(x, y) = maximum(buffer(x - area.x, y - area.y)) - minimum(buffer(x - area.x, y - area.y));

    features(x, y)[4] = range(x, y) / 100;

    Func color_constancy;
    color_constancy(x, y) = sum(buffer(x - area.x, y - area.y, 0) - buffer(x - area.x, y - area.y, 1) + buffer(x - area.x, y - area.y, 1) - buffer(x - area.x, y - area.y, 2) + buffer(x - area.x, y - area.y, 2) - buffer(x - area.x, y - area.y, 0));

    features(x, y)[5] = color_constancy(x, y) / 100;\

}

int main(int argc, char **argv) {

    double t1 = current_time();
    Halide::Runtime::Buffer<uint8_t> left_buffer = Halide::Tools::load_image("images/stereo/bike_smallest.jpg");
    Halide::Runtime::Buffer<uint8_t> right_buffer = Halide::Tools::load_image("images/stereo/bike_smallest.jpg");
    Halide::Runtime::Buffer<float> disparities(left_buffer.width(), left_buffer.height());
    Halide::Runtime::Buffer<float> depth(left_buffer.width(), left_buffer.height());
    Halide::Runtime::Buffer<float> output_x_l(left_buffer.width(), left_buffer.height());
    Halide::Runtime::Buffer<float> output_y_l(left_buffer.width(), left_buffer.height());
    Halide::Runtime::Buffer<float> output_l(left_buffer.width(), left_buffer.height());
    Halide::Runtime::Buffer<float> output_x_r(right_buffer.width(), right_buffer.height());
    Halide::Runtime::Buffer<float> output_y_r(right_buffer.width(), right_buffer.height());
    Halide::Runtime::Buffer<float> output_r(right_buffer.width(), right_buffer.height());

    Halide::Runtime::Buffer<float> delaunay_l[12];
    Halide::Runtime::Buffer<float> delaunay_r[12];

    for (int i = 0; i < 12; i++) {
        Halide::Runtime::Buffer<float> b(left_buffer.width(), left_buffer.height());
        delaunay_l[i] = b;
    }
    for (int i = 0; i < 12; i++) {
        Halide::Runtime::Buffer<float> b(right_buffer.width(), right_buffer.height());
        output_r[i] = b;
    }

    for (int y = 0; y < left_buffer.height(); y++) {
        for (int x = 0; x < left_buffer.width(); x++) {
            disparities(x, y) = f32(0);
        }
    }


    int error_x_l = sobel_x_out(left_buffer, output_x_l);
    int error_y_l = sobel_y_out(left_buffer, output_y_l);
    int error_l = sobel_out(left_buffer, output_l);

    if (error_x_l || error_y_l || error_l) {
        printf("pipeline returned error\n");
    }

    int error_support_l = argmin(output_x_l, output_y_l, left_buffer.width(), left_buffer.height(),
            delaunay_l[0],\
            delaunay_l[1],\
            delaunay_l[2],\
            delaunay_l[3],\
            delaunay_l[4],\
            delaunay_l[5],\
            delaunay_l[6],\
            delaunay_l[7],\
            delaunay_l[8],\
            delaunay_l[9],\
            delaunay_l[10],\
            delaunay_l[11]);

    if (error_support_l) {
        printf("error in support points\n");
    }

    float point_list_l[left_buffer.width() * left_buffer.height() * 4];
    int count = 0;
    for(int xi = 0; xi < left_buffer.width(); xi++) {
        for(int yi = 0; yi < left_buffer.height(); yi++) {
            printf("support points for point: (%d, %d)\n", xi, yi);
            for(int i = 0; i < 4; i++) {
                printf("argmin %d: %f\n", 3 * i, delaunay_l[3 * i](xi, yi));
                printf("argmin %d: %f\n", 3 * i + 1, delaunay_l[3 * i + 1](xi, yi));
                printf("argmin %d: %f\n", 3 * i + 2, delaunay_l[3 * i + 2](xi, yi));
                point_list_l[count] = delaunay_l[3 * i + 1](xi, yi);
                point_list_l[count + 1] = delaunay_l[3 * i + 2](xi, yi);
                count += 2;
            }
        }
    }

    int error_x_r = sobel_x_out(right_buffer, output_x_r);
    int error_y_r = sobel_y_out(right_buffer, output_y_r);
    int error_r = sobel_out(right_buffer, output_r);

    if (error_x_r || error_y_r || error_r) {
        printf("pipeline returned error\n");
    }

    int error_support_r = argmin(output_x_r, output_y_r, left_buffer.width(), left_buffer.height(),
                                 delaunay_r[0],\
            delaunay_r[1],\
            delaunay_r[2],\
            delaunay_r[3],\
            delaunay_r[4],\
            delaunay_r[5],\
            delaunay_r[6],\
            delaunay_r[7],\
            delaunay_r[8],\
            delaunay_r[9],\
            delaunay_r[10],\
            delaunay_r[11]);

    remove_suprious_matches(delaunay_l);

    remove_suprious_matches(delaunay_r);

    if (error_support_r) {
        printf("error in support points r\n");
    }

    float point_list_r[right_buffer.width() * right_buffer.height() * 4];
    count = 0;
    for(int xi = 0; xi < right_buffer.width(); xi++) {
        for(int yi = 0; yi < right_buffer.height(); yi++) {
            printf("support points for point: (%d, %d)\n", xi, yi);
            for(int i = 0; i < 4; i++) {
                printf("argmin %d: %f\n", 3 * i, delaunay_r[3 * i](xi, yi));
                printf("argmin %d: %f\n", 3 * i + 1, delaunay_r[3 * i + 1](xi, yi));
                printf("argmin %d: %f\n", 3 * i + 2, delaunay_r[3 * i + 2](xi, yi));
                point_list_r[count] = delaunay_r[3 * i + 1](xi, yi);
                point_list_r[count + 1] = delaunay_r[3 * i + 2](xi, yi);
                count += 2;
            }
        }
    }

    triangulateio* support_points_l;
    support_points_l->pointlist = point_list_l;
    triangulateio* tri_result_l;
    triangulate((char*) "z", support_points_l, tri_result_l, NULL);


    triangulateio* support_points_r;
    support_points_r->pointlist = point_list_r;
    triangulateio* tri_result_r;
    triangulate((char*) "z", support_points_r, tri_result_r, NULL);

    int num_triangulations = tri_result_l.size();
    Halide::Runtime::Buffer<double> a(num_triangulations);
    Halide::Runtime::Buffer<double> b(num_triangulations);
    Halide::Runtime::Buffer<double> c(num_triangulations);
    interpolate(a, b, c, to_array(tri_result_l));

    
    int err_energy = p(disparities, a, b, c, generate_features(left_buffer, output_l), generate_features(right_buffer, output_r), to_array(support_points_l), to_array(support_points_r));

    if(err_energy) {
        printf("energy pipeline returned error\n");
    }
    int err_depth_map = depth_map(disparities, depth);

    if(err_depth_map) {
        printf("depth map pipeline error\n");
    }
    save_image(depth, "images/stereo/result.png");
}