//#include<iostream>
//#include<cmath>
////#include<random>
//#define NMAX 512
//#define MAXN 100
//#include <chrono>
//using namespace std;
//
//const double EPS = 1E-8;
//const double PI = 3.14159265358979323846;
//
//
//struct Point {
//	double x, y;
//
//	Point() {}
//
//	Point(double x, double y) : x(x), y(y) {}
//};
//
//int sig(double d) { return (d > EPS) - (d < -EPS); }
//
//bool point_same(Point& a, Point& b) {
//	return sig(a.x - b.x) == 0 && sig(a.y - b.y) == 0;
//}
//
////======================hàm custom ? du?i=========================
//
//
//double f_max(double x, double y) {
//	if (x > y) {
//		return x;
//	}
//	else {
//		return y;
//	}
//}
//
//
//
//void copy_points(Point* src, Point* dst, int& n_src, int& n_dst) {
//	for (int i = 0; i < n_src; i++) {
//		dst[i].x = src[i].x;
//		dst[i].y = src[i].y;
//
//	}
//	n_dst = n_src;
//}
//
//void find_all_point(Point* P, int n, Point p, int ptest_index[]) {
//	for (int i = 0; i < n; i++) {
//		ptest_index[i] = -1;
//	}
//
//	int count = 0;
//
//	for (int i = 0; i < n; i++) {
//		if (P[i].x == p.x && P[i].y == p.y) {
//			ptest_index[count++] = i;
//		}
//	}
//}
//
//void insert_point_to_index(Point* P, int& n, int index, Point p) {
//	if (index < 0 || index > n) {
//		return;
//	}
//
//	for (int i = n; i > index; i--) {
//		P[i] = P[i - 1];
//	}
//
//	P[index].x = p.x;
//	P[index].y = p.y;
//
//	n++;
//}
//
//void delete_point_by_index(Point* P, int& n, int index) {
//	if (index < 0 || index >= n) {
//		return;
//	}
//
//	for (int i = index; i < n - 1; i++) {
//		P[i] = P[i + 1];
//	}
//	P[n - 1].x = -9999;
//	P[n - 1].y = -9999;
//	n--;
//}
//
//void move_point_to_end(Point* P, int n, int index) {
//	Point temp = P[index];
//
//	for (int i = index; i < n - 1; i++) {
//		P[i] = P[i + 1];
//	}
//
//	P[n - 1] = temp;
//}
//
//bool allclose(const Point& p1, const Point& p2, double rtol = 1e-5, double atol = 1e-8) {
//	return fabs(p1.x - p2.x) / (atol + rtol * f_max(fabs(p1.x), fabs(p2.x))) <= 1.0 &&
//		fabs(p1.y - p2.y) / (atol + rtol * f_max(fabs(p1.y), fabs(p2.y))) <= 1.0;
//}
//
//void delete_point(Point* P, int& n, Point pdoubt) {
//	int index = -1;
//	for (int i = 0; i < n; i++) {
//		if (P[i].x == pdoubt.x && P[i].y == pdoubt.y) {
//			index = i;
//			break;
//		}
//	}
//
//	if (index != -1) {
//		for (int i = index + 1; i < n; i++) {
//			P[i - 1] = P[i];
//		}
//
//		P[n - 1].x = -9999;
//		P[n - 1].y = -9999;
//		n--;
//	}
//}
//
//int find_index(Point* P, int n, Point pdoubt) {
//	for (int i = 0; i < n; i++) {
//		if (P[i].x == pdoubt.x && P[i].y == pdoubt.y) {
//			return i;
//		}
//	}
//
//	return -1;
//}
//
//void add_two_matrix(double A[1][2], double B[1][2], int m, int n, double C[1][2]) {
//
//	for (int i = 0; i < m; i++) {
//		for (int j = 0; j < n; j++) {
//			C[i][j] = 0;
//		}
//	}
//	for (int i = 0; i < m; i++) {
//		for (int j = 0; j < n; j++) {
//			C[i][j] = A[i][j] + B[i][j];
//		}
//	}
//}
//
//void multiply_matrix_with_double(double matrix[1][2], int n, int m, double scalar, double result[1][2]) {
//
//	for (int i = 0; i < n; i++) {
//		for (int j = 0; j < m; j++) {
//			result[i][j] = matrix[i][j] * scalar;
//		}
//	}
//
//}
//
//void deletePoint(Point* arr, int& n, Point p) {
//	int index = 0;
//	for (int i = 0; i < n; i++) {
//		if (arr[i].x == p.x && arr[i].y == p.y) {
//			index = i;
//			break;
//		}
//	}
//
//	for (int i = index + 1; i < n; i++) {
//		arr[i - 1] = arr[i];
//	}
//
//	arr[n - 1].x = -9999;
//	arr[n - 1].y = -9999;
//
//	n--;
//}
//
//void convert_double_to_point(double matrix[1][2], Point& point) {
//	point.x = matrix[0][0];
//	point.y = matrix[0][1];
//}
//
//double dot_product(double matrix[1][2], Point point) {
//	double dot_product = 0;
//	dot_product += matrix[0][0] * point.x + matrix[0][1] * point.y;
//	return dot_product;
//}
//
//double get_max_value(double mul_dp_xtranspose[1][20], int rows, int n_poly) {
//	double max_value = mul_dp_xtranspose[0][0];
//
//
//	for (int i = 0; i < rows; i++) {
//		for (int j = 0; j < n_poly; j++) {
//			if (mul_dp_xtranspose[i][j] > max_value) {
//				max_value = mul_dp_xtranspose[i][j];
//			}
//		}
//	}
//
//	return max_value;
//}
//void transpose_matrix(double A[][2], int cols, int rows, double A_transpose[][20]) {
//
//	for (int i = 0; i < rows; i++) {
//		for (int j = 0; j < cols; j++) {
//			A_transpose[j][i] = A[i][j];
//		}
//	}
//}
//
//void transpose_dp(double matrix[2][1], double transposed_matrix[1][2]) {
//	transposed_matrix[0][0] = matrix[0][0];
//	transposed_matrix[0][1] = matrix[1][0];
//}
//
//
//void convert_point_to_matrix(Point* points, int n, double X[][2]) {
//	for (int i = 0; i < n; i++) {
//		X[i][0] = points[i].x;
//		X[i][1] = points[i].y;
//	}
//}
//
//void divide_matrix_by_double_and_return_new_matrix(double matrix[2][1], double dp[2][1], int rows, int cols, double d) {
//
//	for (int i = 0; i < rows; i++) {
//		for (int j = 0; j < cols; j++) {
//			dp[i][j] = matrix[i][j] / d;
//		}
//	}
//}
//
//double norm_2(Point p) {
//	double x2 = p.x * p.x;
//	double y2 = p.y * p.y;
//
//	return sqrt(x2 + y2);
//}
//
//void multiply_matrix(double A[2][2], double B[2][1], double C[2][1], int m, int n, int p) {
//
//	for (int i = 0; i < m; i++) {
//		for (int j = 0; j < p; j++) {
//			C[i][j] = 0;
//		}
//	}
//
//	for (int i = 0; i < m; i++) {
//		for (int j = 0; j < p; j++) {
//			for (int k = 0; k < n; k++) {
//				C[i][j] += A[i][k] * B[k][j];
//			}
//		}
//	}
//}
//
//void multiply_matrix(double A[1][2], double B[2][20], double C[1][20], int m, int n, int p) {
//
//	for (int i = 0; i < m; i++) {
//		for (int j = 0; j < p; j++) {
//			C[i][j] = 0;
//		}
//	}
//	for (int i = 0; i < m; i++) {
//		for (int j = 0; j < p; j++) {
//			for (int k = 0; k < n; k++) {
//				C[i][j] += A[i][k] * B[k][j];
//			}
//		}
//	}
//}
//void subtract_points(Point p1, Point p2, Point& p3) {
//	p3.x = p1.x - p2.x;
//	p3.y = p1.y - p2.y;
//}
//
//void convert_point_to_matrix(Point p, double matrix[2][1]) {
//	matrix[0][0] = p.x;
//	matrix[1][0] = p.y;
//}
//void convert_point_to_row_matrix(Point p, double matrix[1][2]) {
//	matrix[0][0] = p.x;
//	matrix[0][1] = p.y;
//}
//
//int find_point_index(Point* Ptest, int n_Ptest, Point pdoubt) {
//
//	int index = -1;
//	for (int i = 0; i < n_Ptest; i++) {
//		if (Ptest[i].x == pdoubt.x && Ptest[i].y == pdoubt.y) {
//			index = i;
//			break;
//		}
//	}
//
//	return index;
//}
//int findMinPointIndex(const Point points[], int n) {
//	if (n <= 0) {
//		return -1;
//	}
//
//	int minIndex = 0;
//
//	for (int i = 1; i < n; ++i) {
//		if (points[i].y < points[minIndex].y) {
//			minIndex = i;
//		}
//		else if (points[i].y == points[minIndex].y && points[i].x < points[minIndex].x) {
//			minIndex = i;
//		}
//	}
//
//	return minIndex;
//}
//
//void moveMinPointsToFront(Point points[], int n) {
//	int minIndex = findMinPointIndex(points, n);
//	Point points_premitive[20];
//	for (int i = 0; i < n; i++) {
//		points_premitive[i] = points[i];
//	}
//	if (minIndex != -1) {
//		Point temp[20];
//		int tempIndex = 0;
//
//		for (int i = minIndex; i < n; ++i) {
//			temp[tempIndex++] = points[i];
//		}
//
//		for (int i = 0; i < minIndex; ++i) {
//			points[i + n - minIndex] = points_premitive[i];
//		}
//
//		for (int i = 0; i < tempIndex; ++i) {
//			points[i] = temp[i];
//		}
//	}
//}
//
//Point* OuterConvexApproximation_and_index(Point* in_poly, int& n_poly, int* points_to_convex_ind) {
//	int n_input = n_poly;
//	Point input_poly[20];
//	for (int i = 0; i < n_input; i++) {
//		input_poly[i].x = in_poly[i].x;
//		input_poly[i].y = in_poly[i].y;
//	}
//
//	//======================code Jarvis thay tu ben duoi===========================
//	double delta = 0.0;
//	double alpha = -PI / 2;
//	double R[2][2];
//
//	R[0][0] = cos(alpha);
//	R[0][1] = sin(alpha);
//	R[1][0] = -sin(alpha);
//	R[1][1] = cos(alpha);
//
//
//	Point D[MAXN];
//
//	D[0] = Point(1.0, 0.0);
//	D[1] = Point(0.0, 1.0);
//	D[2] = Point(-1.0, 0.0);
//	D[3] = Point(0.0, -1.0);
//	int D_size = 4;
//
//
//	double min_x = in_poly[0].x;
//	double min_y = in_poly[0].y;
//	double max_x = in_poly[0].x;
//	double max_y = in_poly[0].y;
//
//	for (int i = 1; i < n_poly; i++) {
//		if (in_poly[i].x < min_x) {
//			min_x = in_poly[i].x;
//		}
//
//		if (in_poly[i].y < min_y) {
//			min_y = in_poly[i].y;
//		}
//
//		if (in_poly[i].x > max_x) {
//			max_x = in_poly[i].x;
//		}
//
//		if (in_poly[i].y > max_y) {
//			max_y = in_poly[i].y;
//		}
//	}
//
//	Point P[MAXN];
//
//	P[0] = Point(max_x, max_y);
//	P[1] = Point(min_x, max_y);
//	P[2] = Point(min_x, min_y);
//	P[3] = Point(max_x, min_y);
//
//	int size_P = 4;
//
//
//	Point Pdoubt[MAXN];
//
//	Pdoubt[0] = Point(max_x, max_y);
//	Pdoubt[1] = Point(min_x, max_y);
//	Pdoubt[2] = Point(min_x, min_y);
//	Pdoubt[3] = Point(max_x, min_y);
//	int size_Pdoubt = 4;
//
//	Point Ptest[MAXN];
//
//	Ptest[0] = Point(max_x, max_y);
//	Ptest[1] = Point(min_x, max_y);
//	Ptest[2] = Point(min_x, min_y);
//	Ptest[3] = Point(max_x, min_y);
//
//	int size_Ptest = 4;
//
//	while (size_Pdoubt > 0) {
//
//		Point pdoubt = Pdoubt[0];
//		int pdoubt_index_idx = find_point_index(Ptest, size_Ptest, pdoubt);
//		int pdoubt_minus_index = (pdoubt_index_idx + size_Ptest - 1) % size_Ptest;
//		int pdoubt_plus_index = (pdoubt_index_idx + 1) % size_Ptest;
//		Point pdoubt_minus = Ptest[pdoubt_minus_index];
//		Point pdoubt_plus = Ptest[pdoubt_plus_index];
//
//		Point result_sub_pminus_plus = Point(0, 0);
//		subtract_points(pdoubt_minus, pdoubt_plus, result_sub_pminus_plus);
//
//		double transposed_matrix[2][1];
//		convert_point_to_matrix(result_sub_pminus_plus, transposed_matrix);
//		double result_mul_matrix[2][1];
//		multiply_matrix(R, transposed_matrix, result_mul_matrix, 2, 2, 1);
//
//		double norm_sub_pminus_pplus = norm_2(result_sub_pminus_plus);
//
//		double dp[2][1];
//		divide_matrix_by_double_and_return_new_matrix(result_mul_matrix, dp, 2, 1, norm_sub_pminus_pplus);
//
//		double X[20][2];
//		convert_point_to_matrix(in_poly, n_poly, X);
//
//		double X_transpose[2][20];
//		transpose_matrix(X, 2, n_poly, X_transpose);
//
//		double dp_transpose[1][2];
//		transpose_dp(dp, dp_transpose);
//
//		double mul_dp_xtranspose[1][20];
//		multiply_matrix(dp_transpose, X_transpose, mul_dp_xtranspose, 1, 2, n_poly);
//
//		double beta_dp = get_max_value(mul_dp_xtranspose, 1, n_poly);
//
//		if (beta_dp == dot_product(dp_transpose, pdoubt_plus)) {
//
//			Point dp_transpose_point = Point(0, 0);
//			convert_double_to_point(dp_transpose, dp_transpose_point);
//			D[D_size] = dp_transpose_point;
//			D_size++;
//			deletePoint(P, size_P, pdoubt);
//			deletePoint(Pdoubt, size_Pdoubt, pdoubt);
//			deletePoint(Ptest, size_Ptest, pdoubt);
//		}
//
//		else if (dot_product(dp_transpose, pdoubt) - beta_dp > delta) {
//
//			double lambda_p = (beta_dp - dot_product(dp_transpose, pdoubt_minus))
//				/ (dot_product(dp_transpose, pdoubt) - dot_product(dp_transpose, pdoubt_minus));
//
//			double pdoubt_minus_convert_to_matrix[1][2];
//			convert_point_to_row_matrix(pdoubt_minus, pdoubt_minus_convert_to_matrix);
//
//			double pdoubt_convert_to_matrix[1][2];
//			convert_point_to_row_matrix(pdoubt, pdoubt_convert_to_matrix);
//
//			double A[1][2];
//			multiply_matrix_with_double(pdoubt_minus_convert_to_matrix, 1, 2, (1 - lambda_p), A);
//
//			double B[1][2];
//			multiply_matrix_with_double(pdoubt_convert_to_matrix, 1, 2, lambda_p, B);
//
//			double p_hat_minus[1][2];
//			add_two_matrix(A, B, 1, 2, p_hat_minus);
//
//			double pdoubt_plus_convert_to_matrix[1][2];
//			convert_point_to_row_matrix(pdoubt_plus, pdoubt_plus_convert_to_matrix);
//
//			double C[1][2];
//			multiply_matrix_with_double(pdoubt_plus_convert_to_matrix, 1, 2, (1 - lambda_p), C);
//
//			double p_hat_plus[1][2];
//			add_two_matrix(C, B, 1, 2, p_hat_plus);
//
//			Point dp_transpose_point = Point(0, 0);
//			convert_double_to_point(dp_transpose, dp_transpose_point);
//			D[D_size] = dp_transpose_point;
//			D_size++;
//
//			int pdoubt_indexp = find_index(P, size_P, pdoubt);
//
//			int pdoubt_index = find_index(Pdoubt, size_Pdoubt, pdoubt);
//
//			delete_point(Pdoubt, size_Pdoubt, pdoubt);
//			Point p_hat_minus_point = Point(0, 0);
//			Point p_hat_plus_point = Point(0, 0);
//			convert_double_to_point(p_hat_minus, p_hat_minus_point);
//			convert_double_to_point(p_hat_plus, p_hat_plus_point);
//			if (allclose(pdoubt, p_hat_minus_point) && allclose(pdoubt, p_hat_plus_point)) {
//
//				int ptest_index = find_point_index(Ptest, size_Ptest, pdoubt);
//				move_point_to_end(Ptest, size_Ptest, ptest_index);
//			}
//			else {
//				delete_point_by_index(P, size_P, pdoubt_indexp);
//				int ptest_index = find_point_index(Ptest, size_Ptest, pdoubt);
//				delete_point_by_index(Ptest, size_Ptest, ptest_index);
//				Point p_hat_plus_point(0, 0);
//				convert_double_to_point(p_hat_plus, p_hat_plus_point);
//				if (allclose(p_hat_plus_point, pdoubt_plus)) {
//				}
//				else {
//					Point p_hat_plus_point = Point(0, 0);
//					convert_double_to_point(p_hat_plus, p_hat_plus_point);
//					insert_point_to_index(P, size_P, pdoubt_indexp, p_hat_plus_point);
//					insert_point_to_index(Pdoubt, size_Pdoubt, pdoubt_index, p_hat_plus_point);
//					insert_point_to_index(Ptest, size_Ptest, ptest_index, p_hat_plus_point);
//				}
//				Point p_hat_minus_point = Point(0, 0);
//				convert_double_to_point(p_hat_minus, p_hat_minus_point);
//				if (allclose(p_hat_minus_point, pdoubt_minus)) {
//				}
//				else {
//					Point p_hat_minus_point = Point(0, 0);
//					convert_double_to_point(p_hat_minus, p_hat_minus_point);
//
//					insert_point_to_index(P, size_P, pdoubt_indexp, p_hat_minus_point);
//					insert_point_to_index(Pdoubt, size_Pdoubt, pdoubt_index, p_hat_minus_point);
//					insert_point_to_index(Ptest, size_Ptest, ptest_index, p_hat_minus_point);
//
//				}
//
//			}
//		}
//		else {
//
//			delete_point(Pdoubt, size_Pdoubt, pdoubt);
//			int ptest_index[MAXN];
//			find_all_point(Ptest, size_Ptest, pdoubt, ptest_index);
//			int first_index = ptest_index[0];
//			move_point_to_end(Ptest, size_Ptest, first_index);
//		}
//
//	}
//
//
//	moveMinPointsToFront(P, size_P);
//	copy_points(P, in_poly, size_P, n_poly);
//	//===============h?t ph?n code Jarvis=================
//
//
//
//
//	for (int i = 0; i < n_poly; i++) {
//		for (int j = 0; j < n_input; j++) {
//			if (point_same(in_poly[i], input_poly[j])) {
//				points_to_convex_ind[i] = j;
//				break;
//			}
//		}
//	}
//
//
//
//	//	===========không copy do?n này vào code chính========
//	return in_poly;
//}
//
////void generateRandomPoints(Point points[], int count, unsigned int seed) {
////	std::mt19937 rng(seed); // Mersenne Twister 19937 generator
////	std::uniform_real_distribution<double> dist(-50.0, 50.0); // Ph?m vi giá tr? t? -10.0 d?n 10.0
////
////	for (int i = 0; i < count; ++i) {
////		points[i].x = dist(rng);
////		points[i].y = dist(rng);
////	}
////}
//
//int main() {
//	int t = 10;
//	while (t-- > 0) {
//		// Ð?t random seed
////	unsigned int seed = 111;
//
//	// Sinh ng?u nhiên 20 di?m v?i random seed và luu vào m?ng Point[]
//		Point points[20] = {
//				Point(1, 4),
//				Point(1, 2),
//				Point(2, 1),
//				Point(3, 2),
//				Point(4, -1),
//				Point(-2, -1),
//				Point(-1, 1),
//				Point(-3, 1),
//				Point(-2 ,3),
//				Point(-3, -1),
//				Point(4,1)
//
//		};
//		//	Point points[20];
//		//	generateRandomPoints(points, 20, seed);
//		int count_tmp = 0;
//		// In k?t qu?
//		int n_poly = 11;
//
//		for (int i = 0; i < n_poly; i++) {
//			std::cout << "[" << points[i].x << ", " << points[i].y << "]\n";
//		}
//		int point_to_convex_indx[20] = { -1, -1, -1, -1, -1, -1, -1, -1, -1,
//										-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
//		OuterConvexApproximation_and_index(points, n_poly, point_to_convex_indx);
//		// In ra giá tr? x, y c?a các ph?n t? trong m?ng
//		std::cout << "======================================================" << std::endl;
//		std::cout << "tap hop bao loi:" << n_poly << " diem" << std::endl;
//		for (int i = 0; i < n_poly; i++) {
//			std::cout << "[" << points[i].x << ", " << points[i].y << "]" << std::endl;
//		}
//		std::cout << "====================================================" << std::endl;
//
//		std::cout << "chi so cac phan tu trong mang la:" << std::endl;
//
//
//		for (int i = 0; i < 20; i++) {
//			std::cout << point_to_convex_indx[i] << " ";
//		}
//		auto start_time = std::chrono::high_resolution_clock::now(); // Bắt đầu đo thời gian
//
//		// Các lệnh trong vòng lặp main
//
//		auto end_time = std::chrono::high_resolution_clock::now(); // Kết thúc đo thời gian
//		auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
//
//		std::cout << "Thời gian thực hiện lần " << 10 - t << ": " << duration.count() << " microseconds" << std::endl;
//
//		// Giả sử chi phí là 1 đơn vị cho mỗi vòng lặp
//		int cost_per_iteration = 1;
//		int total_cost = (10 - t) * cost_per_iteration;
//		std::cout << "Chi phí tích lũy: " << total_cost << std::endl;
//	}
//
//
//}