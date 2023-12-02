//#include<iostream>
//#include <random>
//
//using namespace std;
//
//#define NMAX 512
//#define MAXN 100
//
//
//const double EPS = 1E-8;
//
//
//struct Point {
//    double x, y;
//
//    Point() {}
//
//    Point(double x, double y) : x(x), y(y) {}
//};
//// Hàm sinh ngẫu nhiên Point với random seed
//void generateRandomPoints(Point points[], int count, unsigned int seed) {
//    std::mt19937 rng(seed); // Mersenne Twister 19937 generator
//    std::uniform_real_distribution<double> dist(-50.0, 50.0); // Phạm vi giá trị từ -10.0 đến 10.0
//
//    for (int i = 0; i < count; ++i) {
//        points[i].x = dist(rng);
//        points[i].y = dist(rng);
//    }
//}
//void swap1(Point* a, Point* b) {
//    Point temp;
//    temp.x = a->x;
//    temp.y = a->y;
//
//    a->x = b->x;
//    a->y = b->y;
//
//    b->x = temp.x;
//    b->y = temp.y;
//}
//double cross(Point o, Point a, Point b) {
//    return (a.x - o.x) * (b.y - o.y) - (b.x - o.x) * (a.y - o.y);
//}
//double dis(Point a, Point b) {
//    return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
//}
//void Jarvis(Point* in_poly, int& n_poly) {
//    Point p_max, p_k;
//    int max_index, k_index;
//    int Stack[NMAX] = {}, top1, top2;
//    double sign;
//    Point right_point[10], left_point[10];
//
//    for (int i = 0; i < n_poly; i++) {
//        if (in_poly[i].y < in_poly[0].y ||
//            in_poly[i].y == in_poly[0].y && in_poly[i].x < in_poly[0].x) {
//            Point* j = &(in_poly[0]);
//            Point* k = &(in_poly[i]);
//            swap1(j, k);
//        }
//        if (i == 0) {
//            p_max = in_poly[0];
//            max_index = 0;
//        }
//        if (in_poly[i].y > p_max.y ||
//            in_poly[i].y == p_max.y && in_poly[i].x > p_max.x) {
//            p_max = in_poly[i];
//            max_index = i;
//        }
//    }
//
//    if (max_index == 0) {
//        max_index = 1;
//        p_max = in_poly[max_index];
//    }
//
//    k_index = 0, Stack[0] = 0, top1 = 0;
//    while (k_index != max_index) {
//        p_k = p_max;
//        k_index = max_index;
//        for (int i = 1; i < n_poly; i++) {
//            sign = cross(in_poly[Stack[top1]], in_poly[i], p_k);
//            if ((sign > 0) || ((sign == 0) && (dis(in_poly[Stack[top1]], in_poly[i]) >
//                dis(in_poly[Stack[top1]], p_k)))) {
//                p_k = in_poly[i];
//                k_index = i;
//            }
//        }
//        top1++;
//        Stack[top1] = k_index;
//    }
//    for (int i = 0; i <= top1; i++) right_point[i] = in_poly[Stack[i]];
//
//    k_index = 0, Stack[0] = 0, top2 = 0;
//
//    while (k_index != max_index) {
//        p_k = p_max;
//        k_index = max_index;
//        for (int i = 1; i < n_poly; i++) {
//            sign = cross(in_poly[Stack[top2]], in_poly[i], p_k);
//            if ((sign < 0) || (sign == 0) && (dis(in_poly[Stack[top2]], in_poly[i]) >
//                dis(in_poly[Stack[top2]], p_k))) {
//                p_k = in_poly[i];
//                k_index = i;
//            }
//        }
//        top2++;
//        Stack[top2] = k_index;
//    }
//    for (int i = top2 - 1; i >= 0; i--) left_point[i] = in_poly[Stack[i]];
//
//    for (int i = 0; i < top1 + top2; i++) {
//        if (i <= top1) {
//            in_poly[i] = right_point[i];
//        }
//        else {
//            in_poly[i] = left_point[top2 - (i - top1)];
//        }
//    }
//    n_poly = top1 + top2;
//}
//int main() {
//    // Đặt random seed
//    unsigned int seed = 42;
//
//    // Sinh ngẫu nhiên 20 điểm với random seed và lưu vào mảng Point[]
//    Point points[20];
//    generateRandomPoints(points, 20, seed);
//
//    // In kết quả
//    for (const auto& point : points) {
//        std::cout << "Point: (" << point.x << ", " << point.y << ")\n";
//    }
//    int n_poly = 20;
//    Jarvis(points, n_poly);
//    // In ra giá trị x, y của các phần tử trong mảng
//    std::cout << "======================================================" << std::endl;
//    std::cout << "tap hop bao loi:" << n_poly << " diem" << std::endl;
//    for (int i = 0; i < n_poly; i++) {
//        std::cout << "[" << points[i].x << ", " << points[i].y << "]" << std::endl;
//    }
//    std::cout << "====================================================" << std::endl;
//
//
//
//    return 0;
//}
//
