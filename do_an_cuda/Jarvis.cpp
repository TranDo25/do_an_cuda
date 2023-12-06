#include <iostream>
#include<cmath>
#include<random>
#define NMAX 512
#define MAXN 100

using namespace std;

const double EPS = 1E-8;
const double PI = 3.14159265358979323846;
int sig(double d) { return (d > EPS) - (d < -EPS); }



struct Point {
	double x, y;

	Point() {}

	Point(double x, double y) : x(x), y(y) {}
};
double cross(Point o, Point a, Point b) {
	return (a.x - o.x) * (b.y - o.y) - (b.x - o.x) * (a.y - o.y);
}
bool point_same(Point& a, Point& b) {
	return sig(a.x - b.x) == 0 && sig(a.y - b.y) == 0;
}
double dis(Point a, Point b) {
	return (a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y);
}
void swap1(Point* a, Point* b) {
	Point temp;
	temp.x = a->x;
	temp.y = a->y;

	a->x = b->x;
	a->y = b->y;

	b->x = temp.x;
	b->y = temp.y;
}
void generateRandomPoints(Point points[], int count, unsigned int seed) {
	std::mt19937 rng(seed); // Mersenne Twister 19937 generator
	std::uniform_real_distribution<double> dist(-50.0, 50.0); // Phạm vi giá trị từ -10.0 đến 10.0

	for (int i = 0; i < count; ++i) {
		points[i].x = dist(rng);
		points[i].y = dist(rng);
	}
}
void Jarvis_and_index(Point* in_poly, int& n_poly,
	int* points_to_convex_ind) {
	int n_input = n_poly;
	Point input_poly[20];
	for (int i = 0; i < n_input; i++) {
		input_poly[i].x = in_poly[i].x;
		input_poly[i].y = in_poly[i].y;
	}
	Point p_max, p_k;
	int max_index, k_index;
	int Stack[20], top1, top2;
	double sign;
	Point right_point[10], left_point[10];

	for (int i = 0; i < n_poly; i++) {
		if (in_poly[i].y < in_poly[0].y ||
			in_poly[i].y == in_poly[0].y && in_poly[i].x < in_poly[0].x) {
			Point* j = &(in_poly[0]);
			Point* k = &(in_poly[i]);
			swap1(j, k);
		}
		if (i == 0) {
			p_max = in_poly[0];
			max_index = 0;
		}
		if (in_poly[i].y > p_max.y ||
			in_poly[i].y == p_max.y && in_poly[i].x > p_max.x) {
			p_max = in_poly[i];
			max_index = i;
		}
	}
	if (max_index == 0) {
		max_index = 1;
		p_max = in_poly[max_index];
	}

	k_index = 0, Stack[0] = 0, top1 = 0;
	while (k_index != max_index) {
		p_k = p_max;
		k_index = max_index;
		for (int i = 1; i < n_poly; i++) {
			sign = cross(in_poly[Stack[top1]], in_poly[i], p_k);
			if ((sign > 0) || ((sign == 0) && (dis(in_poly[Stack[top1]], in_poly[i]) >
				dis(in_poly[Stack[top1]], p_k)))) {
				p_k = in_poly[i];
				k_index = i;
			}
		}
		top1++;
		Stack[top1] = k_index;
	}
	for (int i = 0; i <= top1; i++) {
		right_point[i] = in_poly[Stack[i]];
	}

	k_index = 0, Stack[0] = 0, top2 = 0;

	while (k_index != max_index) {
		p_k = p_max;
		k_index = max_index;
		for (int i = 1; i < n_poly; i++) {
			sign = cross(in_poly[Stack[top2]], in_poly[i], p_k);
			if ((sign < 0) || (sign == 0) && (dis(in_poly[Stack[top2]], in_poly[i]) >
				dis(in_poly[Stack[top2]], p_k))) {
				p_k = in_poly[i];
				k_index = i;
			}
		}
		top2++;
		Stack[top2] = k_index;
	}

	for (int i = top2 - 1; i >= 0; i--) {
		left_point[i] = in_poly[Stack[i]];
	}

	for (int i = 0; i < top1 + top2; i++) {
		if (i <= top1) {
			in_poly[i] = right_point[i];
		}
		else {
			in_poly[i] = left_point[top2 - (i - top1)];
		}
	}
	n_poly = top1 + top2;

	for (int i = 0; i < n_poly; i++) {
		for (int j = 0; j < n_input; j++) {
			if (point_same(in_poly[i], input_poly[j])) {
				points_to_convex_ind[i] = j;
				break;
			}
		}
	}
}
int main() {
	// Đặt random seed
	unsigned int seed = 111;

	// Sinh ngẫu nhiên 20 điểm với random seed và lưu vào mảng Point[]
	Point points[20];
	/*	= {
						{1, 4},
						{1, 2},
						{2, 1},
						{3, 2},
						{4, -1},
						{-2, -1},
						{-1, 1},
						{-3, 1},
						{-2 ,3},
						{-3, -1},
						{4,1}

	};*/
	generateRandomPoints(points, 20, seed);
	int count_tmp = 0;
	// In kết quả
	for (const auto& point : points) {
		std::cout << "Point " << count_tmp++ << ": (" << point.x << ", " << point.y << ")\n";
	}
	int n_poly = 20;
	int point_to_convex_indx[20] = { -1, -1, -1, -1, -1, -1, -1, -1, -1,
									-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
	Jarvis_and_index(points, n_poly, point_to_convex_indx);
	// In ra giá trị x, y của các phần tử trong mảng
	std::cout << "======================================================" << std::endl;
	std::cout << "tap hop bao loi:" << n_poly << " diem" << std::endl;
	for (int i = 0; i < n_poly; i++) {
		std::cout << "[" << points[i].x << ", " << points[i].y << "]" << std::endl;
	}
	std::cout << "====================================================" << std::endl;
	cout << "chi so cac phan tu trong mang la" << endl;
	for (int i = 0; i < 20; i++) {
		cout << point_to_convex_indx[i] << " ";
	}
	return 0;

}