#include <iostream>

#define NMAX 512
#define NMIN -999999
#define DBL_MAX 1.7976931348623158e+308



const double EPS = 1E-8;
const double M_PI = 3.14159265358979323846;


struct Point {
	double x, y;

	Point() {}

	Point(double x, double y) : x(x), y(y) {}
};

int sig(double d) { return (d > EPS) - (d < -EPS); }

bool point_same(Point& a, Point& b) {
	return sig(a.x - b.x) == 0 && sig(a.y - b.y) == 0;
}

//======================hàm custom ở dưới=========================
// Hàm tính giai thừa
double factorial(int n) {
    if (n == 0 || n == 1) {
        return 1;
    } else {
        return n * factorial(n - 1);
    }
}
// Hàm tính lũy thừa
double power(double base, int exponent) {
    double result = 1.0;
    for (int i = 0; i < exponent; ++i) {
        result *= base;
    }
    return result;
}
// Hàm tính sin bằng phương pháp Taylor
double sin(double x) {
    int terms = 10;
    double result = 0.0;
    for (int n = 0; n < terms; ++n) {
        // Sử dụng công thức Taylor để tính sin
        result += power(-1, n) * power(x, 2 * n + 1) / factorial(2 * n + 1);
    }
    return result;
}
double sqrt(double x) {
    if (x < 0) {
        return -1;
    }

    double guess = x / 2;
    double previous_guess = 0;

    while (guess != previous_guess) {
        previous_guess = guess;
        guess = (guess + x / guess) / 2;
    }

    return guess;
}
// Hàm tính cos bằng phương pháp Taylor, sử dụng hàm mySin
double cos(double x) {
    // Sử dụng công thức cos(x) = sqrt(1 - sin^2(x))
    return sqrt(1 - sin(x) * sin(x));
}

double f_max(double x, double y) {
    if (x > y) {
        return x;
    }
    else {
        return y;
    }
}
double fabs(double x) {
    if (x >= 0) {
        return x;
    }
    else {
        return -x;
    }
}


//copy phần tử của mảng này chuyển sang mảng khác
void copy_points(Point* src, Point* dst, int& n_src, int& n_dst) {
	// Sao chép từng phần tử của mảng
	for (int i = 0; i < n_src; i++) {
		if (src[i].x != NMIN && src[i].y != NMIN) {
			dst[i].x = src[i].x;
			dst[i].y = src[i].y;
		}

	}
	//thay đổi lại số phần tử của mảng n_poly
	n_dst = n_src;
}

//tìm tất cả các index của Point trong Point*
int* find_all_point(Point* P, int n, Point p) {
	// Khởi tạo mảng kết quả
	int* result = new int[n];
	for (int i = 0; i < n; i++) {
		result[i] = NMIN;
	}

	// Khởi tạo biến đếm
	int count = 0;

	// Duyệt qua mảng P
	for (int i = 0; i < n; i++) {
		// Nếu phần tử P[i] khớp với phần tử cần tìm
		if (P[i].x == p.x && P[i].y == p.y) {
			// Thêm chỉ số của phần tử P[i] vào mảng kết quả
			result[count++] = i;
		}
	}

	// Trả về mảng kết quả
	return result;
}

//chèn Point vào Point* theo index
void insert_point_to_index(Point* P, int& n, int index, Point p) {
	// Kiểm tra index hợp lệ
	if (index < 0 || index > n) {
		return;
	}

	// Di chuyển các phần tử sau vị trí cần chèn lên một vị trí
	for (int i = n; i > index; i--) {
		P[i] = P[i - 1];
	}

	// Chèn phần tử vào vị trí cần chèn
	P[index].x = p.x;
	P[index].y = p.y;

	// Cập nhật kích thước của mảng
	n++;
}

//xoá Point trong Point* theo index
void delete_point_by_index(Point* P, int& n, int index) {
	// Kiểm tra index hợp lệ
	if (index < 0 || index >= n) {
		return;
	}

	// Di chuyển các phần tử sau phần tử cần xóa lên một vị trí
	for (int i = index; i < n - 1; i++) {
		P[i] = P[i + 1];
	}
	P[n - 1].x = NMIN;
	P[n - 1].y = NMIN;
	// Cập nhật kích thước của mảng
	n--;
}

//di chuyển 1 Point xuống cuối con trỏ Point*
void move_point_to_end(Point* P, int n, int index) {
	// Lưu trữ phần tử cần di chuyển
	Point temp = P[index];

	// Di chuyển các phần tử sau phần tử cần di chuyển xuống một vị trí
	for (int i = index; i < n - 1; i++) {
		P[i] = P[i + 1];
	}

	// Đặt phần tử cần di chuyển vào vị trí cuối cùng
	P[n - 1] = temp;
}

//kiểm tra allclose của 2 Point
bool allclose(const Point& p1, const Point& p2, double rtol = 1e-5, double atol = 1e-8) {
	// Tính độ lệch tương đối giữa hai phần tử x
	double rel_diff_x = fabs(p1.x - p2.x) / (atol + rtol * f_max(fabs(p1.x), fabs(p2.x)));

	// Tính độ lệch tương đối giữa hai phần tử y
	double rel_diff_y = fabs(p1.y - p2.y) / (atol + rtol * f_max(fabs(p1.y), fabs(p2.y)));

	// Kiểm tra xem cả hai độ lệch tương đối đều nhỏ hơn hoặc bằng 1.0
	return rel_diff_x <= 1.0 && rel_diff_y <= 1.0;
}

//xoá 1 Point trong Point* không dùng index
void delete_point(Point* P, int& n, Point pdoubt) {
	// Tìm chỉ số của phần tử cần xoá
	int index = -1;
	for (int i = 0; i < n; i++) {
		if (P[i].x == pdoubt.x && P[i].y == pdoubt.y) {
			index = i;
			break;
		}
	}

	// Nếu tìm thấy phần tử cần xoá
	if (index != -1) {
		// Di chuyển các phần tử sau phần tử cần xoá lên một vị trí
		for (int i = index + 1; i < n; i++) {
			P[i - 1] = P[i];
		}

		// Thay thế phần tử cuối cùng của Point* bằng nullptr
		P[n - 1].x = NMIN;
		P[n - 1].y = NMIN;
		//giảm n đi 1 đơn vị
		n--;
	}
}

//tìm index của Point trong Point*
int find_index(Point* P, int n, Point pdoubt) {
	// Duyệt qua tất cả các phần tử của P
	for (int i = 0; i < n; i++) {
		// Nếu phần tử thứ i khớp với pdoubt
		if (P[i].x == pdoubt.x && P[i].y == pdoubt.y) {
			// Trả về chỉ số của phần tử thứ i
			return i;
		}
	}

	// Phần tử pdoubt không tồn tại
	return -1;
}

//cộng hai ma trận cùng chiều với nhau
void add_two_matrix(double A[1][2], double B[1][2], int m, int n, double C[1][2]) {

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			C[i][j] = 0;
		}
	}
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			// Cộng các phần tử tương ứng của hai ma trận
			C[i][j] = A[i][j] + B[i][j];
		}
	}
}

//nhân ma trận double** với 1 số double
void multiply_matrix_with_double(double matrix[1][2], int n, int m, double scalar, double result[1][2]) {

	// Lặp qua các phần tử của ma trận cần nhân và nhân với số double cần nhân
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			result[i][j] = matrix[i][j] * scalar;
		}
	}

}

//xoá Point trong Point*
void deletePoint(Point* arr, int& n, Point p) {
	// Tìm vị trí của Point cần xoá
	int index = 0;
	for (int i = 0; i < n; i++) {
		if (arr[i].x == p.x && arr[i].y == p.y) {
			index = i;
			break;
		}
	}

	// Di chuyển các Point phía sau vị trí cần xoá lên một vị trí
	for (int i = index + 1; i < n; i++) {
		arr[i - 1] = arr[i];
	}

	// Gán giá trị null cho vị trí cần xoá
	arr[n - 1].x = NMIN;
	arr[n - 1].y = NMIN;
	//giảm số phần tử của mảng đi
	n--;
}

//chuyển ma trận 2x1 về Point
void convert_double_to_point(double matrix[1][2], Point& point) {
	// Gán giá trị cho điểm
	point.x = matrix[0][0];
	point.y = matrix[0][1];
}

//tích vo hướng của ma trận và point
double dot_product(double matrix[1][2], Point point) {
	// Khởi tạo tích vô hướng
	double dot_product = 0;
	dot_product += matrix[0][0] * point.x + matrix[0][1] * point.y;
	return dot_product;
}

//lấy ra phần tử lớn nhất trong ma trận
double get_max_value(double mul_dp_xtranspose[1][50], int rows, int n_poly) {
	// Khởi tạo giá trị max
	double max_value = -DBL_MAX;
	int max_index = -1;

	// Lặp qua từng phần tử của mảng
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < n_poly; j++) {
			// So sánh giá trị hiện tại với giá trị max
			if (mul_dp_xtranspose[i][j] > max_value) {
				max_value = mul_dp_xtranspose[i][j];
			}
		}
	}

	// Trả về giá trị max và vị trí của nó
	return max_value;
}
void transpose_matrix(double A[][2], int cols, int rows, double A_transpose[][50]) {

	// Duyệt qua tất cả các phần tử của ma trận A
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			// Gán giá trị của phần tử A[i][j] cho vị trí B[j][i]
			A_transpose[j][i] = A[i][j];
		}
	}
}

//ma trận chuyển vị của con trỏ double**, chuyển dọc sang ngang
void transpose_dp(double matrix[2][1], double transposed_matrix[1][2]) {
	// Gán giá trị cho ma trận chuyển vị
	transposed_matrix[0][0] = matrix[0][0];
	transposed_matrix[0][1] = matrix[1][0];
}


//chuyển Point sang double**, áp dụng cho chuyển đổi mảng in_poly thành X
void convert_point_to_matrix(Point* points, int n, double X[][2]) {
	// Lặp qua từng phần tử Point
	for (int i = 0; i < n; i++) {
		X[i][0] = points[i].x;
		X[i][1] = points[i].y;
	}
}

// chia ma trận cho một số trả về một ma trận mới
void divide_matrix_by_double_and_return_new_matrix(double matrix[2][1], double dp[2][1], int rows, int cols, double d) {

	// Chia từng phần tử của mảng
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			dp[i][j] = matrix[i][j] / d;
		}
	}
}

//tính chuẩn euclid của Point
double norm_2(Point p) {
	// Tính bình phương của từng phần tử của Point
	double x2 = p.x * p.x;
	double y2 = p.y * p.y;

	// Trả về căn bậc hai của tổng bình phương
	return sqrt(x2 + y2);
}

//nhân hai ma trận với nhau
void multiply_matrix(double A[2][2], double B[2][1], double C[2][1], int m, int n, int p) {
	// Khởi tạo ma trận kết quả

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < p; j++) {
			C[i][j] = 0;
		}
	}
	//for (int i = 0; i < m; i++) {
	//	for (int j = 0; j < p; j++) {
	//		cout << C[i][j] << " ";
	//	}
	//	cout << endl;
	//}
	// Nhân hai ma trận
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < p; j++) {
			for (int k = 0; k < n; k++) {
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
}

void multiply_matrix(double A[1][2], double B[2][50], double C[1][50], int m, int n, int p) {
	// Khởi tạo ma trận kết quả

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < p; j++) {
			C[i][j] = 0;
		}
	}
	//for (int i = 0; i < m; i++) {
	//	for (int j = 0; j < p; j++) {
	//		cout << C[i][j] << " ";
	//	}
	//	cout << endl;
	//}
	// Nhân hai ma trận
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < p; j++) {
			for (int k = 0; k < n; k++) {
				C[i][j] += A[i][k] * B[k][j];
			}
		}
	}
}
//trừ hai điểm cho nhau
void subtract_points(Point p1, Point p2, Point& p3) {
	// Hàm trừ 2 Point cho nhau
	p3.x = p1.x - p2.x;
	p3.y = p1.y - p2.y;
}

void convert_point_to_matrix(Point p, double matrix[2][1]) {
	// Gán giá trị cho ma trận
	matrix[0][0] = p.x;
	matrix[1][0] = p.y;
}
//chuyển Point sang ma trận cột chỉ dành cho ma trận trả về kích thước 1x2
void convert_point_to_row_matrix(Point p, double matrix[1][2]) {
	// Gán giá trị cho ma trận
	matrix[0][0] = p.x;
	matrix[0][1] = p.y;
}

//tìm chỉ số của point trong Point*
int find_point_index(Point* Ptest, Point pdoubt) {
	// Hàm trả về chỉ số của Point trong mảng Pest

	int index = -1;
	for (int i = 0; i < NMAX; i++) {
		if (Ptest[i].x == pdoubt.x && Ptest[i].y == pdoubt.y) {
			index = i;
			break;
		}
	}
	if (index == -1) {
		//cout << "khong tim duoc chi so cua phan tu trong mang!" << endl;
	}
	return index;
}

Point* OuterConvexApproximation_and_index(Point* in_poly, int& n_poly, int* points_to_convex_ind) {
	//global
	int n_input = n_poly;
	//global
	Point* input_poly = new Point[51];
	for (int i = 0; i < n_input; i++) {
		input_poly[i].x = i;
		input_poly[i].y = i;
	}
	for (int i = 0; i < n_input; i++) {
		input_poly[i].x = in_poly[i].x;
		input_poly[i].y = in_poly[i].y;
	}
	//global
	double δ = 0.0;
	//khởi tạo mảng xoay R, các hàm sin cos ở dưới cũng dùng thư viện cmath
	double alpha = -M_PI / 2;
	//khởi tạo R
	double R[2][2];

	R[0][0] = cos(alpha);
	R[0][1] = sin(alpha);
	R[1][0] = -sin(alpha);
	R[1][1] = cos(alpha);

	// Tạo một mảng các đối tượng Point
	//nhớ dùng xong phải gán lại về bằng -999 hết
	//Point *D = new Point[NMAX];
	Point D[NMAX];
	for (int i = 0; i < NMAX; i++) {
		D[i].x = NMIN;
		D[i].y = NMIN;
	}
	D[0] = Point(1.0, 0.0);
	D[1] = Point(0.0, 1.0);
	D[2] = Point(-1.0, 0.0);
	D[3] = Point(0.0, -1.0);
	int D_size = 4;


	// Khởi tạo giá trị ban đầu
	double min_x = DBL_MAX;
	double min_y = DBL_MAX;
	double max_x = -DBL_MAX;
	double max_y = -DBL_MAX;

	// Duyệt qua mảng
	for (int i = 0; i < n_poly; i++) {
		// Cập nhật giá trị nhỏ nhất
		if (in_poly[i].x < min_x) {
			min_x = in_poly[i].x;
		}

		// Cập nhật giá trị nhỏ nhất
		if (in_poly[i].y < min_y) {
			min_y = in_poly[i].y;
		}

		// Cập nhật giá trị lớn nhất
		if (in_poly[i].x > max_x) {
			max_x = in_poly[i].x;
		}

		// Cập nhật giá trị lớn nhất
		if (in_poly[i].y > max_y) {
			max_y = in_poly[i].y;
		}
	}

	//khởi tạo P
	Point P[NMAX];
	for (int i = 0; i < NMAX; i++) {
		P[i].x = NMIN;
		P[i].y = NMIN;
	}
	P[0] = Point(max_x, max_y);
	P[1] = Point(min_x, max_y);
	P[2] = Point(min_x, min_y);
	P[3] = Point(max_x, min_y);

	int size_P = 4;


	//khởi tạo Pdoubt
	Point Pdoubt[NMAX];
	for (int i = 0; i < NMAX; i++) {
		Pdoubt[i].x = NMIN;
		Pdoubt[i].y = NMIN;
	}
	Pdoubt[0] = Point(max_x, max_y);
	Pdoubt[1] = Point(min_x, max_y);
	Pdoubt[2] = Point(min_x, min_y);
	Pdoubt[3] = Point(max_x, min_y);
	int size_Pdoubt = 4;

	//khởi tạo Ptest
	Point Ptest[NMAX];
	for (int i = 0; i < NMAX; i++) {
		Ptest[i].x = NMIN;
		Ptest[i].y = NMIN;
	}
	Ptest[0] = Point(max_x, max_y);
	Ptest[1] = Point(min_x, max_y);
	Ptest[2] = Point(min_x, min_y);
	Ptest[3] = Point(max_x, min_y);

	int size_Ptest = 4;


	int count = 0, count4 = 0, count1 = 0, count2 = 0, count3 = 0;


	while (size_Pdoubt > 0) {
		count += 1;
		//không được giải phóng pdoubt vì sẽ làm mất Pdoubt
		Point pdoubt = Pdoubt[0];
		//lấy được chỉ số của pdoubt trong Ptest
		int pdoubt_index_idx = find_point_index(Ptest, pdoubt);
		int pdoubt_minus_index = (pdoubt_index_idx + size_Ptest - 1) % size_Ptest;
		// Lấy chỉ số của mảng liền sau
		int pdoubt_plus_index = (pdoubt_index_idx + 1) % size_Ptest;
		// Lấy phần tử tại vị trí liền trước
		Point pdoubt_minus = Ptest[pdoubt_minus_index];
		// Lấy phần tử tại vị trí liền sau
		Point pdoubt_plus = Ptest[pdoubt_plus_index];

		//bắt đầu tính dp
		Point result_sub_pminus_plus = { 0, 0 };
		subtract_points(pdoubt_minus, pdoubt_plus, result_sub_pminus_plus);

		// Chuyển đổi Point p thành ma trận 2 hàng 1 cột
		double transposed_matrix[2][1];
		convert_point_to_matrix(result_sub_pminus_plus, transposed_matrix);
		//thực hiện nhân ma trận chuyển vị trên với R
		double result_mul_matrix[2][1];
		multiply_matrix(R, transposed_matrix, result_mul_matrix, 2, 2, 1);

		//tính chuẩn norm 2 của kết quả này
		double norm_sub_pminus_pplus = norm_2(result_sub_pminus_plus);

		//dp ở đây có kích thước 2x1, là ma trận cột
		double dp[2][1];
		divide_matrix_by_double_and_return_new_matrix(result_mul_matrix, dp, 2, 1, norm_sub_pminus_pplus);

		//chuyển X về ma trận cho dễ làm việc
		double X[50][2];
		convert_point_to_matrix(in_poly, n_poly, X);

		double X_transpose[2][50];
		transpose_matrix(X, 2, n_poly, X_transpose);

		//dp_tranpose 1x2
		double dp_transpose[1][2];
		transpose_dp(dp, dp_transpose);

		double mul_dp_xtranspose[1][50];
		multiply_matrix(dp_transpose, X_transpose, mul_dp_xtranspose, 1, 2, n_poly);

		//tính Bdp
		double βdp = get_max_value(mul_dp_xtranspose, 1, n_poly);

		//cần tính tích vô hướng trong điều kiện dưới
		if (βdp == dot_product(dp_transpose, pdoubt_plus)) {
			count1 += 1;
			//chuyển đổi dp_transpose về kiểu Point để nạp vào tập D
			Point dp_transpose_point = { 0,0 };
			convert_double_to_point(dp_transpose, dp_transpose_point);
			D[D_size] = dp_transpose_point;
			D_size++;
			deletePoint(P, size_P, pdoubt);
			deletePoint(Pdoubt, size_Pdoubt, pdoubt);
			deletePoint(Ptest, size_Ptest, pdoubt);

		}
		//chuyển đổi pdoubt sang dạng vector<double<double>>

		else if (dot_product(dp_transpose, pdoubt) - βdp > δ) {
			count2++;
			//tính λp
			double λp = (βdp - dot_product(dp_transpose, pdoubt_minus))
				/ (dot_product(dp_transpose, pdoubt) - dot_product(dp_transpose, pdoubt_minus));
			//1x2 nhưng thực tế cần 2x1
			//đây thực chất là 1x2 nhưng trên danh nghĩa là 2x1

			double pdoubt_minus_convert_to_matrix[1][2];
			convert_point_to_row_matrix(pdoubt_minus, pdoubt_minus_convert_to_matrix);

			double pdoubt_convert_to_matrix[1][2];
			convert_point_to_row_matrix(pdoubt, pdoubt_convert_to_matrix);

			double A[1][2];
			multiply_matrix_with_double(pdoubt_minus_convert_to_matrix, 1, 2, (1 - λp), A);

			double B[1][2];
			multiply_matrix_with_double(pdoubt_convert_to_matrix, 1, 2, λp, B);

			double p_hat_minus[1][2];
			add_two_matrix(A, B, 1, 2, p_hat_minus);

			double pdoubt_plus_convert_to_matrix[1][2];
			convert_point_to_row_matrix(pdoubt_plus, pdoubt_plus_convert_to_matrix);

			double C[1][2];
			multiply_matrix_with_double(pdoubt_plus_convert_to_matrix, 1, 2, (1 - λp), C);

			double p_hat_plus[1][2];
			add_two_matrix(C, B, 1, 2, p_hat_plus);

			Point dp_transpose_point = { 0,0 };;
			convert_double_to_point(dp_transpose, dp_transpose_point);
			D[D_size] = dp_transpose_point;
			D_size++;

			int pdoubt_indexp = find_index(P, size_P, pdoubt);

			int pdoubt_index = find_index(Pdoubt, size_Pdoubt, pdoubt);

			delete_point(Pdoubt, size_Pdoubt, pdoubt);
			Point p_hat_minus_point = { 0,0 };
			Point p_hat_plus_point = { 0,0 };
			convert_double_to_point(p_hat_minus, p_hat_minus_point);
			convert_double_to_point(p_hat_plus, p_hat_plus_point);
			if (allclose(pdoubt, p_hat_minus_point) && allclose(pdoubt, p_hat_plus_point)) {
				count4++;
				int ptest_index = find_point_index(Ptest, pdoubt);
				move_point_to_end(Ptest, size_Ptest, ptest_index);
			}
			else {
				delete_point_by_index(P, size_P, pdoubt_indexp);
				int ptest_index = find_point_index(Ptest, pdoubt);
				delete_point_by_index(Ptest, size_Ptest, ptest_index);
				Point p_hat_plus_point = { 0,0 };
				convert_double_to_point(p_hat_plus, p_hat_plus_point);
				if (allclose(p_hat_plus_point, pdoubt_plus)) {
					//cout << "1" << endl;
				}
				else {
					Point p_hat_plus_point = { 0,0 };
					convert_double_to_point(p_hat_plus, p_hat_plus_point);
					insert_point_to_index(P, size_P, pdoubt_indexp, p_hat_plus_point);
					insert_point_to_index(Pdoubt, size_Pdoubt, pdoubt_index, p_hat_plus_point);
					insert_point_to_index(Ptest, size_Ptest, ptest_index, p_hat_plus_point);
				}
				Point p_hat_minus_point = { 0,0 };
				convert_double_to_point(p_hat_minus, p_hat_minus_point);
				if (allclose(p_hat_minus_point, pdoubt_minus)) {
					//cout << "2" << endl;
				}
				else {
					Point p_hat_minus_point = { 0,0 };
					convert_double_to_point(p_hat_minus, p_hat_minus_point);

					insert_point_to_index(P, size_P, pdoubt_indexp, p_hat_minus_point);
					insert_point_to_index(Pdoubt, size_Pdoubt, pdoubt_index, p_hat_minus_point);
					insert_point_to_index(Ptest, size_Ptest, ptest_index, p_hat_minus_point);

				}

			}
		}
		else {
			count3++;
			delete_point(Pdoubt, size_Pdoubt, pdoubt);
			int* ptest_index = find_all_point(Ptest, size_Ptest, pdoubt);
			int first_index = ptest_index[0];
			move_point_to_end(Ptest, size_Ptest, first_index);
		}

	}
	copy_points(P, in_poly, size_P, n_poly);






	for (int i = 0; i < n_poly; i++) {
		for (int j = 0; j < n_input; j++) {
			if (point_same(in_poly[i], input_poly[j])) {
				points_to_convex_ind[i] = j;
				break;
			}
		}
	}

	return in_poly;
}


int main() {

//	Point* in_poly = new Point[50];
//
//	in_poly[0] = { -36.8462, -4.13499 };
//	in_poly[1] = { -28.1041, 17.8865 };
//	in_poly[2] = { 43.4693, 1.94164 };
//	in_poly[3] = { -46.5428, 2.97002 };
//	in_poly[4] = { -49.2302, -43.3158 };
//	in_poly[5] = { 18.6773, 43.0436 };
//	in_poly[6] = { 2.69288, 15.3919 };
//	in_poly[7] = { 20.1191, 26.2198 };
//	in_poly[8] = { -45.2535, -17.1766 };



    Point in_poly[9] = { {0, 0}, {1, 1}, {2, 0}, {2, 2}, {1, 3}, {0, 2}, {1, 9}, {1, 8}, {2, 5} };
	//    points[9] = {25.641, -13.4661};
	//    points[10] = {48.255, 25.3356};
	//    points[11] = {-42.7314, 38.4707};
	//    points[12] = {-6.35886, -2.22682};
	//    points[13] = {-22.5093, -33.3493};
	//    points[14] = {39.7656, -43.9436};
	//    points[15] = {0.452289, -18.0967};
	//    points[16] = {-0.602331, -40.9267};
	//    points[17] = {-42.6251, -11.5858};
	//    points[18] = {41.3817, -3.55542};
	//    points[19] = {-44.9916, 27.0205};
	//    points[20] = {-37.4635, 18.8455};
	//    points[21] = {12.9543, 22.5412};
	//    points[22] = {38.8572, -19.3678};
	//    points[23] = {1.32737, 34.5982};
	//    points[24] = {34.1511, -8.46054};
	//    points[25] = {-3.20826, -32.1672};
	//    points[26] = {7.16548, -46.6946};
	//    points[27] = {-0.151988, 24.8293};
	//    points[28] = {39.0737, 34.204};
	//    points[29] = {-28.7248, -36.9573};
	//    points[30] = {-22.5412, -8.57067};
	//    points[31] = {20.982, -26.0089};
	//    points[32] = {-18.246, 15.2059};
	//    points[33] = {18.1346, -11.2275};
	//    points[34] = {-35.2467, 34.5576};
	//    points[35] = {45.5409, -35.1848};
	//    points[36] = {-9.12333, 6.48987};
	//    points[37] = {-1.14855, 46.1095};
	//    points[38] = {-30.0243, 12.9269};
	//    points[39] = {15.1254, 30.3073};
	//    points[40] = {-2.35682, -29.675};
	//    points[41] = {40.1673, -35.7979};
	//    points[42] = {-8.9687, 38.5648};
	//    points[43] = {-33.7801, -13.4661};
	//    points[44] = {-36.4891, -4.46927};
	//    points[45] = {-4.76998, 43.1674};
	//    points[46] = {-28.4752, 40.8922};
	//    points[47] = {36.086, 0.595588};
	//    points[48] = {31.7561, -3.7755};
	//    points[49] = {13.2739, 32.4697};


	int n_poly = 9;

	int points_to_convex_ind[9] = { -1, -1, -1, -1, -1, -1, -1, -1, -1 };
    std::cout << "======================================================" << std::endl;
    std::cout << "tap hop cac diem point:" << std::endl;
	for (int i = 0; i < n_poly; i++) {
        std::cout << "[" << in_poly[i].x << ", " << in_poly[i].y << "]" << std::endl;
	}
    std::cout << "====================================================" << std::endl;

	Point* result = OuterConvexApproximation_and_index(in_poly, n_poly, points_to_convex_ind);
	// In ra giá trị x, y của các phần tử trong mảng
	std::cout << "======================================================" << std::endl;
    std::cout << "tap hop bao loi:" << n_poly << " diem" << std::endl;
	for (int i = 0; i < n_poly; i++) {
        std::cout << "[" << result[i].x << ", " << result[i].y << "]" << std::endl;
	}
    std::cout << "====================================================" << std::endl;

	// Giải phóng bộ nhớ

    std::cout << "chi so cua cac diem convex so voi tap goc:" << std::endl;
	for (int i = 0; i < n_poly; i++) std::cout << points_to_convex_ind[i] << " ";
//	delete[] points;

	return 0;
}


