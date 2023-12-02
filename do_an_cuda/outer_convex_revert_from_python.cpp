//#include<iostream>
//#include <random>
//#include<cfloat>
//#include<cmath>
//using namespace std;
//
//#define NMAX 512
//#define MAXN 100
//
//
//const double EPS = 1E-8;
//
//const double PI = 3.14159265358979323846;
//
//struct Point {
//	double x, y;
//
//	Point() {}
//
//	Point(double x, double y) : x(x), y(y) {}
//};
////=================================================================
////
//// // Hàm tính giai thừa
////double factorial(int n) {
////	if (n == 0 || n == 1) {
////		return 1;
////	}
////	else {
////		return n * factorial(n - 1);
////	}
////}
////// Hàm tính lũy thừa
////double power(double base, int exponent) {
////	double result = 1.0;
////	for (int i = 0; i < exponent; ++i) {
////		result *= base;
////	}
////	return result;
////}
//// Hàm tính sin bằng phương pháp Taylor
////double sin(double x) {
////	int terms = 10;
////	double result = 0.0;
////	for (int n = 0; n < terms; ++n) {
////		// Sử dụng công thức Taylor để tính sin
////		result += power(-1, n) * power(x, 2 * n + 1) / factorial(2 * n + 1);
////	}
////	return result;
////}
////double sqrt(double x) {
////	if (x < 0) {
////		return -1;
////	}
////
////	double guess = x / 2;
////	double previous_guess = 0;
////
////	while (guess != previous_guess) {
////		previous_guess = guess;
////		guess = (guess + x / guess) / 2;
////	}
////
////	return guess;
////}
////// Hàm tính cos bằng phương pháp Taylor, sử dụng hàm mySin
////double cos(double x) {
////	return sqrt(1 - sin(x) * sin(x));
////}
//
//double f_max(double x, double y) {
//	if (x > y) {
//		return x;
//	}
//	else {
//		return y;
//	}
//}
//double fabs(double x) {
//	//nễu x > 0, trả về giá trị dương, nếu x < 0, trả về giá trị âm
//	if (x >= 0) {
//		return x;
//	}
//	else {
//		return -x;
//	}
//}
//
//
////copy phần tử của mảng này chuyển sang mảng khác
//void copy_points(Point* src, Point* dst, int& n_src, int& n_dst) {
//	// Sao chép từng phần tử của mảng
//	for (int i = 0; i < n_src; i++) {
//		if (src[i].x != DBL_MIN && src[i].y != DBL_MIN) {
//			dst[i].x = src[i].x;
//			dst[i].y = src[i].y;
//		}
//
//	}
//	//thay đổi lại số phần tử của mảng n_poly
//	n_dst = n_src;
//}
//
////tìm tất cả các index của Point trong Point*
//void find_all_point(Point* P, int n, Point p, int ptest_index[]) {
//	// Khởi tạo mảng kết quả
//	for (int i = 0; i < n; i++) {
//		ptest_index[i] = DBL_MIN;
//	}
//
//	// Khởi tạo biến đếm
//	int count = 0;
//
//	// Duyệt qua mảng P
//	for (int i = 0; i < n; i++) {
//		// Nếu phần tử P[i] khớp với phần tử cần tìm
//		if (P[i].x == p.x && P[i].y == p.y) {
//			// Thêm chỉ số của phần tử P[i] vào mảng kết quả
//			ptest_index[count++] = i;
//		}
//	}
//}
//
////chèn Point vào Point* theo index
//void insert_point_to_index(Point* P, int& n, int index, Point p) {
//	// Kiểm tra index hợp lệ
//	if (index < 0 || index > n) {
//		return;
//	}
//
//	// Di chuyển các phần tử sau vị trí cần chèn lên một vị trí
//	for (int i = n; i > index; i--) {
//		P[i] = P[i - 1];
//	}
//
//	// Chèn phần tử vào vị trí cần chèn
//	P[index].x = p.x;
//	P[index].y = p.y;
//
//	// Cập nhật kích thước của mảng
//	n++;
//}
//
////xoá Point trong Point* theo index
//void delete_point_by_index(Point* P, int& n, int index) {
//	// Kiểm tra index hợp lệ
//	if (index < 0 || index >= n) {
//		return;
//	}
//
//	// Di chuyển các phần tử sau phần tử cần xóa lên một vị trí
//	for (int i = index; i < n - 1; i++) {
//		P[i] = P[i + 1];
//	}
//	P[n - 1].x = DBL_MIN;
//	P[n - 1].y = DBL_MIN;
//	// Cập nhật kích thước của mảng
//	n--;
//}
//
////di chuyển 1 Point xuống cuối con trỏ Point*
//void move_point_to_end(Point* P, int n, int index) {
//	// Lưu trữ phần tử cần di chuyển
//	Point temp = P[index];
//
//	// Di chuyển các phần tử sau phần tử cần di chuyển xuống một vị trí
//	for (int i = index; i < n - 1; i++) {
//		P[i] = P[i + 1];
//	}
//
//	// Đặt phần tử cần di chuyển vào vị trí cuối cùng
//	P[n - 1] = temp;
//}
//
////kiểm tra allclose của 2 Point
//bool allclose(const Point& p1, const Point& p2, double rtol = 1e-5, double atol = 1e-8) {
//	// Tính độ lệch tương đối giữa hai phần tử x
//	double rel_diff_x = fabs(p1.x - p2.x) / (atol + rtol * f_max(fabs(p1.x), fabs(p2.x)));
//
//	// Tính độ lệch tương đối giữa hai phần tử y
//	double rel_diff_y = fabs(p1.y - p2.y) / (atol + rtol * f_max(fabs(p1.y), fabs(p2.y)));
//
//	// Kiểm tra xem cả hai độ lệch tương đối đều nhỏ hơn hoặc bằng 1.0
//	return rel_diff_x <= 1.0 && rel_diff_y <= 1.0;
//}
//
////xoá 1 Point trong Point* không dùng index
//void delete_point(Point* P, int& n, Point pdoubt) {
//	// Tìm chỉ số của phần tử cần xoá
//	int index = -1;
//	for (int i = 0; i < n; i++) {
//		if (P[i].x == pdoubt.x && P[i].y == pdoubt.y) {
//			index = i;
//			break;
//		}
//	}
//
//	// Nếu tìm thấy phần tử cần xoá
//	if (index != -1) {
//		// Di chuyển các phần tử sau phần tử cần xoá lên một vị trí
//		for (int i = index + 1; i < n; i++) {
//			P[i - 1] = P[i];
//		}
//
//		// Thay thế phần tử cuối cùng của Point* bằng nullptr
//		P[n - 1].x = DBL_MIN;
//		P[n - 1].y = DBL_MIN;
//		//giảm n đi 1 đơn vị
//		n--;
//	}
//}
//
////tìm index của Point trong Point*
//int find_index(Point* P, int n, Point pdoubt) {
//	// Duyệt qua tất cả các phần tử của P
//	for (int i = 0; i < n; i++) {
//		// Nếu phần tử thứ i khớp với pdoubt
//		if (P[i].x == pdoubt.x && P[i].y == pdoubt.y) {
//			// Trả về chỉ số của phần tử thứ i
//			return i;
//		}
//	}
//
//	// Phần tử pdoubt không tồn tại
//	return -1;
//}
//
////cộng hai ma trận cùng chiều với nhau
//void add_two_matrix(double A[1][2], double B[1][2], int m, int n, double C[1][2]) {
//
//	for (int i = 0; i < m; i++) {
//		for (int j = 0; j < n; j++) {
//			C[i][j] = 0;
//		}
//	}
//	for (int i = 0; i < m; i++) {
//		for (int j = 0; j < n; j++) {
//			// Cộng các phần tử tương ứng của hai ma trận
//			C[i][j] = A[i][j] + B[i][j];
//		}
//	}
//}
//
////nhân ma trận double** với 1 số double
//void multiply_matrix_with_double(double matrix[1][2], int n, int m, double scalar, double result[1][2]) {
//
//	// Lặp qua các phần tử của ma trận cần nhân và nhân với số double cần nhân
//	for (int i = 0; i < n; i++) {
//		for (int j = 0; j < m; j++) {
//			result[i][j] = matrix[i][j] * scalar;
//		}
//	}
//
//}
//
////xoá Point trong Point*
//void deletePoint(Point* arr, int& n, Point p) {
//	// Tìm vị trí của Point cần xoá
//	int index = 0;
//	for (int i = 0; i < n; i++) {
//		if (arr[i].x == p.x && arr[i].y == p.y) {
//			index = i;
//			break;
//		}
//	}
//
//	// Di chuyển các Point phía sau vị trí cần xoá lên một vị trí
//	for (int i = index + 1; i < n; i++) {
//		arr[i - 1] = arr[i];
//	}
//
//	// Gán giá trị null cho vị trí cần xoá
//	arr[n - 1].x = DBL_MIN;
//	arr[n - 1].y = DBL_MIN;
//	//giảm số phần tử của mảng đi
//	n--;
//}
//
////chuyển ma trận 2x1 về Point
//void convert_double_to_point(double matrix[1][2], Point& point) {
//	// Gán giá trị cho điểm
//	point.x = matrix[0][0];
//	point.y = matrix[0][1];
//}
//
////tích vo hướng của ma trận và point
//double dot_product(double matrix[1][2], Point point) {
//	// Khởi tạo tích vô hướng
//	double dot_product = 0;
//	dot_product += matrix[0][0] * point.x + matrix[0][1] * point.y;
//	return dot_product;
//}
//
////lấy ra phần tử lớn nhất trong ma trận
//double get_max_value(double mul_dp_xtranspose[1][20], int rows, int n_poly) {
//	// Khởi tạo giá trị max
//	double max_value = -DBL_MAX;
//
//
//	// Lặp qua từng phần tử của mảng
//	for (int i = 0; i < rows; i++) {
//		for (int j = 0; j < n_poly; j++) {
//			// So sánh giá trị hiện tại với giá trị max
//			if (mul_dp_xtranspose[i][j] > max_value) {
//				max_value = mul_dp_xtranspose[i][j];
//			}
//		}
//	}
//
//	// Trả về giá trị max và vị trí của nó
//	return max_value;
//}
//void transpose_matrix(double A[][2], int cols, int rows, double A_transpose[][20]) {
//
//	// Duyệt qua tất cả các phần tử của ma trận A
//	for (int i = 0; i < rows; i++) {
//		for (int j = 0; j < cols; j++) {
//			// Gán giá trị của phần tử A[i][j] cho vị trí B[j][i]
//			A_transpose[j][i] = A[i][j];
//		}
//	}
//}
//
////ma trận chuyển vị của con trỏ double**, chuyển dọc sang ngang
//void transpose_dp(double matrix[2][1], double transposed_matrix[1][2]) {
//	// Gán giá trị cho ma trận chuyển vị
//	transposed_matrix[0][0] = matrix[0][0];
//	transposed_matrix[0][1] = matrix[1][0];
//}
//
//
////chuyển Point sang double**, áp dụng cho chuyển đổi mảng in_poly thành X
//void convert_point_to_matrix(Point* points, int n, double X[][2]) {
//	// Lặp qua từng phần tử Point
//	for (int i = 0; i < n; i++) {
//		X[i][0] = points[i].x;
//		X[i][1] = points[i].y;
//	}
//}
//
//// chia ma trận cho một số trả về một ma trận mới
//void divide_matrix_by_double_and_return_new_matrix(double matrix[2][1], double dp[2][1], int rows, int cols, double d) {
//
//	// Chia từng phần tử của mảng
//	for (int i = 0; i < rows; i++) {
//		for (int j = 0; j < cols; j++) {
//			dp[i][j] = matrix[i][j] / d;
//		}
//	}
//}
//
////tính chuẩn euclid của Point
//double norm_2(Point p) {
//	// Tính bình phương của từng phần tử của Point
//	double x2 = p.x * p.x;
//	double y2 = p.y * p.y;
//
//	// Trả về căn bậc hai của tổng bình phương
//	return sqrt(x2 + y2);
//}
//
////nhân hai ma trận với nhau
//void multiply_matrix(double A[2][2], double B[2][1], double C[2][1], int m, int n, int p) {
//	// Khởi tạo ma trận kết quả
//
//	for (int i = 0; i < m; i++) {
//		for (int j = 0; j < p; j++) {
//			C[i][j] = 0;
//		}
//	}
//	//for (int i = 0; i < m; i++) {
//	//	for (int j = 0; j < p; j++) {
//	//		cout << C[i][j] << " ";
//	//	}
//	//	cout << endl;
//	//}
//	// Nhân hai ma trận
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
//	// Khởi tạo ma trận kết quả
//
//	for (int i = 0; i < m; i++) {
//		for (int j = 0; j < p; j++) {
//			C[i][j] = 0;
//		}
//	}
//	//for (int i = 0; i < m; i++) {
//	//	for (int j = 0; j < p; j++) {
//	//		cout << C[i][j] << " ";
//	//	}
//	//	cout << endl;
//	//}
//	// Nhân hai ma trận
//	for (int i = 0; i < m; i++) {
//		for (int j = 0; j < p; j++) {
//			for (int k = 0; k < n; k++) {
//				C[i][j] += A[i][k] * B[k][j];
//			}
//		}
//	}
//}
////trừ hai điểm cho nhau
//void subtract_points(Point p1, Point p2, Point& p3) {
//	// Hàm trừ 2 Point cho nhau
//	p3.x = p1.x - p2.x;
//	p3.y = p1.y - p2.y;
//}
//
//void convert_point_to_matrix(Point p, double matrix[2][1]) {
//	// Gán giá trị cho ma trận
//	matrix[0][0] = p.x;
//	matrix[1][0] = p.y;
//}
////chuyển Point sang ma trận cột chỉ dành cho ma trận trả về kích thước 1x2
//void convert_point_to_row_matrix(Point p, double matrix[1][2]) {
//	// Gán giá trị cho ma trận
//	matrix[0][0] = p.x;
//	matrix[0][1] = p.y;
//}
//
////tìm chỉ số của point trong Point*
//int find_point_index(Point* Ptest, Point pdoubt) {
//	// Hàm trả về chỉ số của Point trong mảng Pest
//
//	int index = -1;
//	for (int i = 0; i < MAXN; i++) {
//		if (Ptest[i].x == pdoubt.x && Ptest[i].y == pdoubt.y) {
//			index = i;
//			break;
//		}
//	}
//	if (index == -1) {
//		//cout << "khong tim duoc chi so cua phan tu trong mang!" << endl;
//	}
//	return index;
//}
//void swap1(Point* a, Point* b) {
//	Point temp;
//	temp.x = a->x;
//	temp.y = a->y;
//
//	a->x = b->x;
//	a->y = b->y;
//
//	b->x = temp.x;
//	b->y = temp.y;
//}
////=======================kết thúc sửa code======================
//Point* OuterConvexApproximation(Point* in_poly, int& n_poly) {
//	Point p_max, p_k;
//	int max_index, k_index;
//	int Stack[NMAX] = {}, top1, top2;
//	double sign;
//	Point right_point[10], left_point[10];
//
//	for (int i = 0; i < n_poly; i++) {
//		if (in_poly[i].y < in_poly[0].y ||
//			in_poly[i].y == in_poly[0].y && in_poly[i].x < in_poly[0].x) {
//			Point* j = &(in_poly[0]);
//			Point* k = &(in_poly[i]);
//			swap1(j, k);
//		}
//		if (i == 0) {
//			p_max = in_poly[0];
//			max_index = 0;
//		}
//		if (in_poly[i].y > p_max.y ||
//			in_poly[i].y == p_max.y && in_poly[i].x > p_max.x) {
//			p_max = in_poly[i];
//			max_index = i;
//		}
//	}
//	//global
//	double delta = 0.0;
//	//khởi tạo mảng xoay R, các hàm sin cos ở dưới cũng dùng thư viện cmath
//	double alpha = -PI / 2;
//	//khởi tạo R
//	double R[2][2];
//
//	R[0][0] = cos(alpha);
//	R[0][1] = sin(alpha);
//	R[1][0] = -sin(alpha);
//	R[1][1] = cos(alpha);
//	// Tạo một mảng các đối tượng Point
//
//	Point D[MAXN];
//	for (int i = 0; i < MAXN; i++) {
//		D[i].x = DBL_MIN;
//		D[i].y = DBL_MIN;
//	}
//	D[0] = Point(1.0, 0.0);
//	D[1] = Point(0.0, 1.0);
//	D[2] = Point(-1.0, 0.0);
//	D[3] = Point(0.0, -1.0);
//	int D_size = 4;
//
//
//	// Khởi tạo giá trị ban đầu
//	double min_x = DBL_MAX;
//	double min_y = DBL_MAX;
//	double max_x = -DBL_MAX;
//	double max_y = -DBL_MAX;
//
//	// Duyệt qua mảng
//	for (int i = 0; i < n_poly; i++) {
//		// Cập nhật giá trị nhỏ nhất
//		if (in_poly[i].x < min_x) {
//			min_x = in_poly[i].x;
//		}
//
//		// Cập nhật giá trị nhỏ nhất
//		if (in_poly[i].y < min_y) {
//			min_y = in_poly[i].y;
//		}
//
//		// Cập nhật giá trị lớn nhất
//		if (in_poly[i].x > max_x) {
//			max_x = in_poly[i].x;
//		}
//
//		// Cập nhật giá trị lớn nhất
//		if (in_poly[i].y > max_y) {
//			max_y = in_poly[i].y;
//		}
//	}
//
//	//khởi tạo P
//	Point P[MAXN];
//	for (int i = 0; i < MAXN; i++) {
//		P[i].x = DBL_MIN;
//		P[i].y = DBL_MIN;
//	}
//	P[0] = Point(max_x, max_y);
//	P[1] = Point(min_x, max_y);
//	P[2] = Point(min_x, min_y);
//	P[3] = Point(max_x, min_y);
//
//	int size_P = 4;
//
//
//	//khởi tạo Pdoubt
//	Point Pdoubt[MAXN];
//	for (int i = 0; i < MAXN; i++) {
//		Pdoubt[i].x = DBL_MIN;
//		Pdoubt[i].y = DBL_MIN;
//	}
//	Pdoubt[0] = Point(max_x, max_y);
//	Pdoubt[1] = Point(min_x, max_y);
//	Pdoubt[2] = Point(min_x, min_y);
//	Pdoubt[3] = Point(max_x, min_y);
//	int size_Pdoubt = 4;
//
//	//khởi tạo Ptest
//	Point Ptest[MAXN];
//	for (int i = 0; i < MAXN; i++) {
//		Ptest[i].x = DBL_MIN;
//		Ptest[i].y = DBL_MIN;
//	}
//	Ptest[0] = Point(max_x, max_y);
//	Ptest[1] = Point(min_x, max_y);
//	Ptest[2] = Point(min_x, min_y);
//	Ptest[3] = Point(max_x, min_y);
//
//	int size_Ptest = 4;
//
//	int count = 0, count4 = 0, count1 = 0, count2 = 0, count3 = 0;
//
//
//	while (size_Pdoubt > 0) {
//		count += 1;
//		//không được giải phóng pdoubt vì sẽ làm mất Pdoubt
//		Point pdoubt = Pdoubt[0];
//		//lấy được chỉ số của pdoubt trong Ptest
//		int pdoubt_index_idx = find_point_index(Ptest, pdoubt);
//		int pdoubt_minus_index = (pdoubt_index_idx + size_Ptest - 1) % size_Ptest;
//		// Lấy chỉ số của mảng liền sau
//		int pdoubt_plus_index = (pdoubt_index_idx + 1) % size_Ptest;
//		// Lấy phần tử tại vị trí liền trước
//		Point pdoubt_minus = Ptest[pdoubt_minus_index];
//		// Lấy phần tử tại vị trí liền sau
//		Point pdoubt_plus = Ptest[pdoubt_plus_index];
//
//		//bắt đầu tính dp
//		Point result_sub_pminus_plus = { 0, 0 };
//		subtract_points(pdoubt_minus, pdoubt_plus, result_sub_pminus_plus);
//
//		// Chuyển đổi Point p thành ma trận 2 hàng 1 cột
//		double transposed_matrix[2][1];
//		convert_point_to_matrix(result_sub_pminus_plus, transposed_matrix);
//		//thực hiện nhân ma trận chuyển vị trên với R
//		double result_mul_matrix[2][1];
//		multiply_matrix(R, transposed_matrix, result_mul_matrix, 2, 2, 1);
//
//		//tính chuẩn norm 2 của kết quả này
//		double norm_sub_pminus_pplus = norm_2(result_sub_pminus_plus);
//
//		//dp ở đây có kích thước 2x1, là ma trận cột
//		double dp[2][1];
//		divide_matrix_by_double_and_return_new_matrix(result_mul_matrix, dp, 2, 1, norm_sub_pminus_pplus);
//
//		//chuyển X về ma trận cho dễ làm việc
//		double X[20][2];
//		convert_point_to_matrix(in_poly, n_poly, X);
//
//		double X_transpose[2][20];
//		transpose_matrix(X, 2, n_poly, X_transpose);
//
//		//dp_tranpose 1x2
//		double dp_transpose[1][2];
//		transpose_dp(dp, dp_transpose);
//
//		double mul_dp_xtranspose[1][20];
//		multiply_matrix(dp_transpose, X_transpose, mul_dp_xtranspose, 1, 2, n_poly);
//
//		//tính Bdp
//		double beta_dp = get_max_value(mul_dp_xtranspose, 1, n_poly);
//
//		//cần tính tích vô hướng trong điều kiện dưới
//		if (beta_dp == dot_product(dp_transpose, pdoubt_plus)) {
//			count1 += 1;
//			//chuyển đổi dp_transpose về kiểu Point để nạp vào tập D
//			Point dp_transpose_point = { 0,0 };
//			convert_double_to_point(dp_transpose, dp_transpose_point);
//			D[D_size] = dp_transpose_point;
//			D_size++;
//			deletePoint(P, size_P, pdoubt);
//			deletePoint(Pdoubt, size_Pdoubt, pdoubt);
//			deletePoint(Ptest, size_Ptest, pdoubt);
//
//		}
//		//chuyển đổi pdoubt sang dạng vector<double<double>>
//
//		else if (dot_product(dp_transpose, pdoubt) - beta_dp > delta) {
//			count2++;
//			//tính lambda_p
//			double lambda_p = (beta_dp - dot_product(dp_transpose, pdoubt_minus))
//				/ (dot_product(dp_transpose, pdoubt) - dot_product(dp_transpose, pdoubt_minus));
//			//1x2 nhưng thực tế cần 2x1
//			//đây thực chất là 1x2 nhưng trên danh nghĩa là 2x1
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
//			Point dp_transpose_point = { 0,0 };;
//			convert_double_to_point(dp_transpose, dp_transpose_point);
//			D[D_size] = dp_transpose_point;
//			D_size++;
//
//			int pdoubt_indexp = find_index(P, size_P, pdoubt);
//
//			int pdoubt_index = find_index(Pdoubt, size_Pdoubt, pdoubt);
//
//			delete_point(Pdoubt, size_Pdoubt, pdoubt);
//			Point p_hat_minus_point = { 0,0 };
//			Point p_hat_plus_point = { 0,0 };
//			convert_double_to_point(p_hat_minus, p_hat_minus_point);
//			convert_double_to_point(p_hat_plus, p_hat_plus_point);
//			if (allclose(pdoubt, p_hat_minus_point) && allclose(pdoubt, p_hat_plus_point)) {
//				count4++;
//				int ptest_index = find_point_index(Ptest, pdoubt);
//				move_point_to_end(Ptest, size_Ptest, ptest_index);
//			}
//			else {
//				delete_point_by_index(P, size_P, pdoubt_indexp);
//				int ptest_index = find_point_index(Ptest, pdoubt);
//				delete_point_by_index(Ptest, size_Ptest, ptest_index);
//				Point p_hat_plus_point(0, 0);
//				convert_double_to_point(p_hat_plus, p_hat_plus_point);
//				if (allclose(p_hat_plus_point, pdoubt_plus)) {
//					//cout << "1" << endl;
//				}
//				else {
//					Point p_hat_plus_point = { 0,0 };
//					convert_double_to_point(p_hat_plus, p_hat_plus_point);
//					insert_point_to_index(P, size_P, pdoubt_indexp, p_hat_plus_point);
//					insert_point_to_index(Pdoubt, size_Pdoubt, pdoubt_index, p_hat_plus_point);
//					insert_point_to_index(Ptest, size_Ptest, ptest_index, p_hat_plus_point);
//				}
//				Point p_hat_minus_point = { 0,0 };
//				convert_double_to_point(p_hat_minus, p_hat_minus_point);
//				if (allclose(p_hat_minus_point, pdoubt_minus)) {
//					//cout << "2" << endl;
//				}
//				else {
//					Point p_hat_minus_point = { 0,0 };
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
//			count3++;
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
//	copy_points(P, in_poly, size_P, n_poly);
//
//
//	//=================không copy phần dưới===================
//	return in_poly;
//}
//
//// Hàm sinh ngẫu nhiên Point với random seed
//void generateRandomPoints(Point points[], int count, unsigned int seed) {
//	std::mt19937 rng(seed); // Mersenne Twister 19937 generator
//	std::uniform_real_distribution<double> dist(-50.0, 50.0); // Phạm vi giá trị từ -10.0 đến 10.0
//
//	for (int i = 0; i < count; ++i) {
//		points[i].x = dist(rng);
//		points[i].y = dist(rng);
//	}
//}
//int main() {
//		// Đặt random seed
//	unsigned int seed = 1000;
//
//	// Sinh ngẫu nhiên 20 điểm với random seed và lưu vào mảng Point[]
//	Point points[20];
//	generateRandomPoints(points, 20, seed);
//
//	// In kết quả
//	for (const auto& point : points) {
//		std::cout << "Point: (" << point.x << ", " << point.y << ")\n";
//	}
//	int n_poly = 20;
//	std::cout << "======================================================" << std::endl;
//	std::cout << "tap hop cac diem point:" << std::endl;
//	for (int i = 0; i < n_poly; i++) {
//		std::cout << "[" << points[i].x << ", " << points[i].y << "]" << std::endl;
//	}
//	std::cout << "====================================================" << std::endl;
//
//	Point* result = OuterConvexApproximation(points, n_poly);
//	// In ra giá trị x, y của các phần tử trong mảng
//	std::cout << "======================================================" << std::endl;
//	std::cout << "tap hop bao loi:" << n_poly << " diem" << std::endl;
//	for (int i = 0; i < n_poly; i++) {
//		std::cout << "[" << result[i].x << ", " << result[i].y << "]" << std::endl;
//	}
//	std::cout << "====================================================" << std::endl;
//
//	
//
//	return 0;
//}
//
//
