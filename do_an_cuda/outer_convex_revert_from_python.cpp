#define _USE_MATH_DEFINES
#include <utility>
#include <vector>
#include <iostream>
#include <vector>
#include <algorithm>
#include <tuple>
#include <cstdlib>  // Để sử dụng hàm rand()
#include <ctime>    // Để sử dụng hàm srand()
#include <iostream>
#include <random>
//const double PI = 3.14159265358979323846;
#define MAXN 100
#define NMAX 512
#define NMIN -9999
using namespace std;
struct Point {
	double x, y;
	Point() {}
	Point(double x, double y) : x(x), y(y) {}
};
Point convert_double_to_point(double** matrix) {
	// Khởi tạo điểm
	Point point = { 0, 0 };

	// Gán giá trị cho điểm
	point.x = matrix[0][0];
	point.y = matrix[1][0];

	// Trả về điểm
	return point;
}
double dot_product(double** matrix, Point point) {
	// Khởi tạo tích vô hướng
	double dot_product = 0;
	dot_product += matrix[0][0] * point.x + matrix[0][1] * point.y;
	return dot_product;
}
double get_max_value(double** mul_dp_xtranspose, int rows, int n_poly) {
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

double** transpose_dp(double** matrix) {
	// Khởi tạo ma trận chuyển vị
	double** transposed_matrix = new double* [1];
	transposed_matrix[0] = new double[2];

	// Gán giá trị cho ma trận chuyển vị
	transposed_matrix[0][0] = matrix[0][0];
	transposed_matrix[0][1] = matrix[1][0];

	// Trả về ma trận chuyển vị
	return transposed_matrix;
}
double** convert_point_to_matrix(Point* points, int n) {
	// Khởi tạo ma trận 2 chiều
	double** matrix = new double* [n];
	for (int i = 0; i < n; i++) {
		matrix[i] = new double[2];
	}

	// Lặp qua từng phần tử Point
	for (int i = 0; i < n; i++) {
		matrix[i][0] = points[i].x;
		matrix[i][1] = points[i].y;
	}

	return matrix;
}

double** divide_matrix_by_double_and_return_new_matrix(double** matrix, int rows, int cols, double d) {
	// Tạo mảng double** mới
	double** new_matrix = new double* [rows];
	for (int i = 0; i < rows; i++) {
		new_matrix[i] = new double[cols];
	}

	// Chia từng phần tử của mảng
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			new_matrix[i][j] = matrix[i][j] / d;
		}
	}

	// Trả về mảng double** mới
	return new_matrix;
}
double norm_2(Point p) {
	// Tính bình phương của từng phần tử của Point
	double x2 = p.x * p.x;
	double y2 = p.y * p.y;

	// Trả về căn bậc hai của tổng bình phương
	return sqrt(x2 + y2);
}

double** multiply_matrix(double** A, double** B, int m, int n, int p) {
	// Khởi tạo ma trận kết quả
	double** C = new double* [m];
	for (int i = 0; i < m; i++) {
		C[i] = new double[p];
	}
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

	return C;
}
Point subtract_points(Point p1, Point p2) {
	// Hàm trừ 2 Point cho nhau

	Point p3;
	p3.x = p1.x - p2.x;
	p3.y = p1.y - p2.y;

	return p3;
}
double** convert_point_to_matrix(Point p) {
	// Khởi tạo ma trận 2 hàng 1 cột

	double** matrix = new double* [2];
	matrix[0] = new double[1];
	matrix[1] = new double[1];

	// Gán giá trị cho ma trận
	matrix[0][0] = p.x;
	matrix[1][0] = p.y;

	return matrix;
}
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
		cout << "khong tim duoc chi so cua phan tu trong mang!" << endl;
	}
	return index;
}

int size_array_point(Point* D) {
	// Hàm trả về số phần tử trong mảng có giá trị x, y đều khác -9999

	int count = 0;
	for (int i = 0; i < NMAX; i++) {
		if (D[i].x != -9999 && D[i].y != -9999) {
			count++;
		}
	}

	return count;
}
void copy_vector_to_array(Point* dest, const vector<Point>& src, int size) {
	for (int i = 0; i < size; i++) {
		dest[i].x = src[i].x;
		dest[i].y = src[i].y;
	}
}
//void remove_extra_points(Point*& points, size_t size, const std::vector<Point>& reference_points) {
//	// Tìm chỉ số của phần tử đầu tiên không cần xóa
//	size_t new_size = std::remove_if(points, points + size, [&reference_points](Point* point) {
//		return std::find(reference_points.begin(), reference_points.end(), *point) == reference_points.end();
//		}) - points;
//
//	// Tạo một mảng tạm để lưu trữ các Point cần giữ lại
//	Point* output = new Point[new_size];
//
//	// Copy các Point cần giữ lại sang mảng tạm
//	for (size_t i = 0; i < new_size; i++) {
//		output[i] = points[i];
//	}
//
//	// Thay thế mảng con trỏ bằng mảng tạm
//	delete[] points;
//	points = output;
//}

std::vector<Point> to_point_vector(const std::vector<std::tuple<double, double>>& input) {
	std::vector<Point> output;
	for (const auto& tuple : input) {
		double x = std::get<0>(tuple);
		double y = std::get<1>(tuple);
		output.push_back(Point(x, y));
	}
	return output;
}
void PrintVector(vector<tuple<double, double>>& result) {
	cout << "[";
	for (int i = 0; i < result.size(); i++) {
		cout << "(" << get<0>(result[i]) << ", " << get<1>(result[i]) << ")";
		if (i < result.size() - 1) {
			cout << ", " << endl;
		}
	}
	cout << "]" << endl;
}

vector<tuple<double, double>> GetEnoughPoints(Point* in_poly, int n_poly) {
	// Khởi tạo mảng kết quả
	vector<tuple<double, double>> out_poly(n_poly);

	// Lặp qua mảng in_poly
	for (int i = 0; i < n_poly; i++) {
		// Lấy ra phần tử thứ i trong mảng in_poly
		out_poly[i] = make_tuple(in_poly[i].x, in_poly[i].y);
	}

	return out_poly;
}
vector<int> find_list_index(vector<tuple<double, double>> Ptest, tuple<double, double> pdoubt) {
	// Tạo mảng kết quả chứa các chỉ số tìm thấy
	vector<int> result;

	// Lặp qua từng phần tử trong mảng Ptest
	for (int i = 0; i < Ptest.size(); i++) {
		// Kiểm tra xem phần tử hiện tại có phải là pdoubt không
		if (Ptest[i] == pdoubt) {
			// Thêm chỉ số hiện tại vào mảng kết quả
			result.push_back(i);
		}
	}

	// Trả về mảng kết quả
	return result;
}
std::vector<std::tuple<double, double>> insert(std::vector<std::tuple<double, double>> P, int pdoubt_index, std::tuple<double, double> p_hat_plus) {
	// Chèn phần tử mới vào vị trí thứ pdoubt_index
	P.insert(P.begin() + pdoubt_index, p_hat_plus);
	return P;
}
vector<tuple<double, double>> delete_tuple_by_index(vector<tuple<double, double>>& v, int index) {
	// Tạo một vector mới để lưu kết quả
	vector<tuple<double, double>> v_result;

	// Duyệt qua vector gốc
	for (int i = 0; i < v.size(); i++) {
		// Nếu chỉ số không phải là index cần xóa
		if (i != index) {
			// Thêm phần tử vào vector mới
			v_result.push_back(v[i]);
		}
	}

	// Trả về vector mới
	return v_result;
}
vector<tuple<double, double>> dao_vi_tri_ptest(vector<tuple<double, double>>& Ptest, int ptest_index) {
	// Lưu trữ phần tử có chỉ số là ptest_index
	tuple<double, double> ptest_to_move = Ptest[ptest_index];

	// Di chuyển tất cả các phần tử sau ptest_index xuống một vị trí
	for (int i = ptest_index; i < Ptest.size() - 1; i++) {
		Ptest[i] = Ptest[i + 1];
	}

	// Chèn ptest_to_move vào cuối mảng
	Ptest.back() = ptest_to_move;

	return Ptest;
}
template <typename T>
bool allclose(const std::vector<std::vector<T>>& v1, const std::vector<std::vector<T>>& v2, double rtol = 1e-5, double atol = 1e-8) {
	// Kiểm tra xem hai mảng có cùng kích thước hay không
	if (v1.size() != v2.size()) {
		return false;
	}

	// Duyệt qua hai mảng
	for (int i = 0; i < v1.size(); i++) {
		// Kiểm tra xem hai mảng con có cùng kích thước hay không
		if (v1[i].size() != v2[i].size()) {
			return false;
		}

		// Duyệt qua hai mảng con
		for (int j = 0; j < v1[i].size(); j++) {
			// Nếu hai giá trị không gần bằng nhau trong phạm vi sai số cho trước thì trả về false
			if (std::abs(v1[i][j] - v2[i][j]) > rtol * std::abs(v1[i][j]) + atol) {
				return false;
			}
		}
	}

	// Nếu hai mảng gần bằng nhau trong phạm vi sai số cho trước thì trả về true
	return true;
}
template <typename T>
bool allclose(const std::vector<T>& v1, const std::vector<T>& v2, double rtol = 1e-5, double atol = 1e-8) {
	// Kiểm tra xem hai mảng có cùng kích thước hay không
	if (v1.size() != v2.size()) {
		return false;
	}

	// Duyệt qua hai mảng
	for (int i = 0; i < v1.size(); i++) {
		// Nếu hai giá trị không gần bằng nhau trong phạm vi sai số cho trước thì trả về false
		if (std::abs(v1[i] - v2[i]) > rtol * std::abs(v1[i]) + atol) {
			return false;
		}
	}

	// Nếu hai mảng gần bằng nhau trong phạm vi sai số cho trước thì trả về true
	return true;
}
int tim_chi_so_phan_tu_tuple(const vector<tuple<double, double>>& P, const tuple<double, double>& pdoubt) {
	// Khởi tạo biến chỉ số
	int index = -1;

	// Duyệt qua mảng
	for (int i = 0; i < P.size(); i++) {
		// So sánh phần tử tuple tại vị trí i với pdoubt
		if (P[i] == pdoubt) {
			// Nếu tìm thấy, cập nhật chỉ số
			index = i;
			break;
		}
	}

	// Trả về chỉ số
	return index;
}
vector<vector<double>> operator+(const vector<vector<double>>& m1, const vector<vector<double>>& m2) {
	// Kiểm tra kích thước của 2 ma trận
	if (m1.size() != m2.size() || m1[0].size() != m2[0].size()) {
		// Trả về ma trận rỗng nếu 2 ma trận không có cùng kích thước
		return vector<vector<double>>();
	}

	// Tạo ma trận kết quả
	vector<vector<double>> result(m1.size());
	for (size_t i = 0; i < m1.size(); i++) {
		result[i].resize(m1[i].size());
	}

	// Cộng 2 ma trận
	for (size_t i = 0; i < m1.size(); i++) {
		for (size_t j = 0; j < m1[i].size(); j++) {
			result[i][j] = m1[i][j] + m2[i][j];
		}
	}

	return result;
}
std::vector<std::vector<double>> multiplyMatrixWithANumber(const std::vector<std::vector<double>>& matrix, double scalar) {
	// Tạo ma trận kết quả
	std::vector<std::vector<double>> result(matrix.size());
	for (size_t i = 0; i < matrix.size(); i++) {
		result[i].resize(matrix[i].size());
	}

	// Nhân ma trận với số
	for (size_t i = 0; i < matrix.size(); i++) {
		for (size_t j = 0; j < matrix[i].size(); j++) {
			result[i][j] = matrix[i][j] * scalar;
		}
	}

	return result;
}
std::vector<std::tuple<double, double>> deleteTuple(std::vector<std::tuple<double, double>>& v, std::tuple<double, double>& t) {
	// Tìm vị trí của tuple cần xoá
	auto it = std::find(v.begin(), v.end(), t);

	// Xoá tuple nếu tìm thấy
	if (it != v.end()) {
		v.erase(it);
	}

	return v;
}
std::vector<std::vector<double>> deleteVectorRow(std::vector<std::vector<double>>& matrix, std::vector<double>& pdoubt) {
	// Tìm các hàng cần xoá
	std::vector<int> rows_to_delete;
	for (int i = 0; i < matrix.size(); i++) {
		if (std::equal(matrix[i].begin(), matrix[i].end(), pdoubt.begin())) {
			rows_to_delete.push_back(i);
		}
	}

	// Xoá các hàng đã tìm thấy
	for (int i = rows_to_delete.size() - 1; i >= 0; i--) {
		matrix.erase(matrix.begin() + rows_to_delete[i]);
	}

	return matrix;
}
std::tuple<double, double> convertMatrixToTuple(const std::vector<std::vector<double>> matrix) {
	// Kiểm tra xem ma trận có đúng kích thước không
	if (matrix.size() != 1 || matrix[0].size() != 2) {
		std::cerr << "Ma tran phai co kich thuoc 1x2" << std::endl;
		return std::tuple<double, double>(0.0, 0.0);
	}

	// Lấy phần tử đầu tiên của ma trận
	double a = matrix[0][0];

	// Lấy phần tử thứ hai của ma trận
	double b = matrix[0][1];

	// Trả về tuple
	return std::tuple<double, double>(a, b);
}
double dotProduct(const std::vector<std::vector<double>>& v1, const std::vector<std::vector<double>>& v2) {
	// Kiểm tra xem 2 vector có cùng kích thước không
	if (v1.size() != v2.size() || v1[0].size() != v2[0].size()) {
		std::cerr << "2 vector phải có cùng kích thước" << std::endl;
		return -1;
	}

	// Tính tích vô hướng của 2 vector
	double dot_product = 0.0;
	for (size_t i = 0; i < v1.size(); i++) {
		for (size_t j = 0; j < v1[0].size(); j++) {
			dot_product += v1[i][j] * v2[i][j];
		}
	}

	return dot_product;
}
std::vector<std::vector<double>> convertTupleToMatrix(std::tuple<double, double>& tuple) {
	// Xác định số cột (số phần tử trong tuple)
	size_t cols = std::tuple_size< std::tuple<double, double>>::value;

	// Khởi tạo ma trận với kích thước 1 hàng và số cột đã xác định
	std::vector<std::vector<double>> matrix(1, std::vector<double>(cols, 0.0));
	matrix[0][0] = std::get<0>(tuple);
	matrix[0][1] = std::get<1>(tuple);


	return matrix;
}
std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& A) {
	// Khởi tạo ma trận chuyển vị với kích thước đã xác định
	size_t rows = A.size();
	size_t cols = A[0].size();
	std::vector<std::vector<double>> transposed_matrix(cols, std::vector<double>(rows, 0.0));

	// Sao chép dữ liệu từ ma trận ban đầu sang ma trận chuyển vị
	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < cols; j++) {
			transposed_matrix[j][i] = A[i][j];
		}
	}

	return transposed_matrix;
}
double getMax(const std::vector<std::vector<double>>& A) {
	// Khởi tạo biến max để lưu trữ giá trị max
	double max = -INFINITY;

	// Duyệt qua từng phần tử của ma trận
	for (size_t i = 0; i < A.size(); i++) {
		for (size_t j = 0; j < A[i].size(); j++) {
			// Nếu phần tử hiện tại lớn hơn max thì cập nhật max
			if (A[i][j] > max) {
				max = A[i][j];
			}
		}
	}

	return max;
}
std::vector<std::tuple<double, double>> to_matrix(const Point* points, int size) {
	std::vector<std::tuple<double, double>> X(size);
	for (int i = 0; i < size; i++) {
		X[i] = std::make_tuple(points[i].x, points[i].y);
	}
	return X;
}
template <typename T>
double norm(const std::tuple<T, T>& t) {
	return std::sqrt(std::get<0>(t) * std::get<0>(t) + std::get<1>(t) * std::get<1>(t));
}

template <typename T>
std::vector<std::vector<T>> divide_matrix(const std::vector<std::vector<T>>& matrix, T divisor) {
	std::vector<std::vector<T>> result(matrix.size());
	for (int i = 0; i < matrix.size(); i++) {
		result[i] = std::vector<T>(matrix[i].size());
		for (int j = 0; j < matrix[i].size(); j++) {
			result[i][j] = matrix[i][j] / divisor;
		}
	}
	return result;
}
std::vector<std::vector<double>> multiplyMatrices(const std::vector<std::vector<double>> matrix1, const std::vector<std::vector<double>> matrix2) {
	int rows1 = matrix1.size();
	int cols1 = matrix1[0].size();
	int rows2 = matrix2.size();
	int cols2 = matrix2[0].size();

	// Kiểm tra xem có thể nhân hai ma trận này với nhau không
	if (cols1 != rows2) {
		std::cerr << "Khong the nhan hai ma tran nay voi nhau." << std::endl;
		return std::vector<std::vector<double>>();
	}

	// Khởi tạo ma trận kết quả với kích thước đúng
	std::vector<std::vector<double>> result(rows1, std::vector<double>(cols2, 0.0));

	for (int i = 0; i < rows1; ++i) {
		for (int j = 0; j < cols2; ++j) {
			for (int k = 0; k < cols1; ++k) {
				result[i][j] += matrix1[i][k] * matrix2[k][j];
			}
		}
	}

	return result;
}
//chuyển vị tuple X sang dạng mảng 2 chiều
vector<vector<double>> transposeTupleVector(const  vector< tuple<double, double>>& tuple_vector) {
	// Xác định số hàng (số lượng phần tử trong vector)
	size_t rows = tuple_vector.size();

	// Xác định số cột (số phần tử trong mỗi tuple)
	size_t cols = rows > 0 ? tuple_size< tuple<double, double>>::value : 0;

	// Khởi tạo ma trận chuyển vị với kích thước đã xác định
	vector< vector<double>> transposed_matrix(cols, vector<double>(rows, 0.0));

	// Sao chép dữ liệu từ vector vào ma trận chuyển vị
	for (size_t i = 0; i < rows; i++) {
		for (size_t j = 0; j < cols; j++) {
			transposed_matrix[j][i] = j == 0 ? get<0>(tuple_vector[i]) : get<1>(tuple_vector[i]);
		}
	}

	return transposed_matrix;
}

// Hàm trừ hai tuple kiểu  tuple<double, double>
tuple<double, double> subtract_tuples(const  tuple<double, double>& t1, const  tuple<double, double>& t2) {
	double result1 = get<0>(t1) - get<0>(t2);
	double result2 = get<1>(t1) - get<1>(t2);
	return  make_tuple(result1, result2);
}

Point* OuterConvexApproximation(Point* in_poly, int& n_poly, double δ = 0.0) {
	//khởi tạo mảng xoay R, các hàm sin cos ở dưới cũng dùng thư viện cmath
	double alpha = -M_PI / 2;
	//khởi tạo R
	double** R = new double* [2];
	for (int i = 0; i < 2; i++) {
		R[i] = new double[2];
	}
	R[0][0] = cos(alpha);
	R[0][1] = sin(alpha);
	R[1][0] = -sin(alpha);
	R[1][1] = cos(alpha);

	// Tạo một mảng các đối tượng Point
	//nhớ dùng xong phải gán lại về bằng -999 hết
	Point* D = new Point[NMAX];
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
	for (size_t i = 0; i < n_poly; i++) {
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

	Point r1 = Point(max_x, max_y);
	Point r2 = Point(min_x, max_y);
	Point r3 = Point(min_x, min_y);
	Point r4 = Point(max_x, min_y);

	Point* P = new Point[NMAX];
	for (int i = 0; i < NMAX; i++) {
		P[i].x = NMIN;
		P[i].y = NMIN;
	}
		P[0] = r1;
		P[1] = r2;
		P[2] = r3;
		P[3] = r4;

	Point* Pdoubt = P;
	Point* Ptest = P;
	int count = 0, count4 = 0, count1 = 0, count2 = 0, count3 = 0;
	while (size_array_point(Pdoubt) > 0) {
		count += 1;
		Point pdoubt = Pdoubt[0];
		//lấy được chỉ số của pdoubt trong Ptest
		int pdoubt_index_idx = find_point_index(Ptest, pdoubt);
		
		int pdoubt_minus_index = (pdoubt_index_idx + size_array_point(Ptest) - 1) % size_array_point(Ptest);

		// Lấy chỉ số của mảng liền sau
		int pdoubt_plus_index = (pdoubt_index_idx + 1) % size_array_point(Ptest);

		// Lấy phần tử tại vị trí liền trước
		Point pdoubt_minus = Ptest[pdoubt_minus_index];

		// Lấy phần tử tại vị trí liền sau
		Point pdoubt_plus = Ptest[pdoubt_plus_index];

		//bắt đầu tính dp
		Point result_sub_pminus_plus = subtract_points(pdoubt_minus, pdoubt_plus);
		//Point vector_result_sub_pminus_plus = new Point[1];
		//vector_result_sub_pminus_plus[0] = result_sub_pminus_plus;

		//vector<vector<double>> transposed_matrix = transposeTupleVector(vector_result_sub_pminus_plus);
		//chuyển đổi thành ma trận cột 2x1, viết hàm chuyển đổi Point thành ma trận cột 
		// Chuyển đổi Point p thành ma trận 2 hàng 1 cột
		double** transposed_matrix = convert_point_to_matrix(result_sub_pminus_plus);
		//thực hiện nhân ma trận chuyển vị trên với R
		double** result_mul_matrix = multiply_matrix(R, transposed_matrix, 2, 2, 1);

		//std::vector<std::vector<double>> result_mul_matrix = multiplyMatrices(R, transposed_matrix);
		//tính chuẩn norm 2 của kết quả này
		double norm_sub_pminus_pplus = norm_2(result_sub_pminus_plus);
		//dp ở đây có kích thước 2x1, là ma trận cột 
		double** dp = divide_matrix_by_double_and_return_new_matrix(result_mul_matrix, 2, 1, norm_sub_pminus_pplus);


		//chuyển in_poly về vector X cho dễ tính
		//vector<tuple<double, double>> n_poly_vector = GetEnoughPoints(in_poly, n_poly);
		//=================làm tiếp ở đây=====================
		double** X_transpose = convert_point_to_matrix(in_poly, n_poly);

		//vector<vector<double>> X_transpose = transposeTupleVector(n_poly_vector);
		//dp đang là ma trận 2x1, vậy nếu muốn nhân với ma trận có kích thước 2x50 thì phải chuyển vị ma trận dp đi
		
		double** dp_transpose = transpose_dp(dp);

		
		//vector<vector<double>> mul_dp_xtranspose = multiplyMatrices(dp_transpose, X_transpose);
		double** mul_dp_xtranspose = multiply_matrix(dp_transpose, X_transpose, 1, 2, n_poly);
		double βdp = get_max_value(mul_dp_xtranspose,1,n_poly);

		////convert pdoubt_plus sang vector<vector<double>>  để có thể tái sử dụng hàm nhân ma trận
		//vector<vector<double>>  pdoubt_plus_vector = convertTupleToMatrix(pdoubt_plus);

		//cần tính tích vô hướng trong điều kiện dưới
		if (βdp == dot_product(dp_transpose, pdoubt_plus)) {
			count1 += 1;
			//convert dp to tuple để thêm vào tập D
			/*vector<vector<double>> dp_transpose = transpose(dp);
			tuple<double, double> dp_tuple = convertMatrixToTuple(dp_transpose);*/

			//chuyển đổi dp_transpose về kiểu Point để nạp vào tập D
			Point dp_transpose_point = convert_double_to_point(dp_transpose);
			D[D_size] = dp_transpose_point;
			D_size++;

			P = deleteTuple(P, pdoubt);
			Pdoubt = deleteTuple(Pdoubt, pdoubt);
			Ptest = deleteTuple(Ptest, pdoubt);
		}
		//chuyển đổi pdoubt sang dạng vector<double<double>>
		else if (dotProduct(dp_transpose, convertTupleToMatrix(pdoubt)) - βdp > δ) {
			count2++;
			//tính λp
			double λp = (βdp - dotProduct(dp_transpose, convertTupleToMatrix(pdoubt_minus)))
				/ (dotProduct(dp_transpose, convertTupleToMatrix(pdoubt)) - dotProduct(dp_transpose, convertTupleToMatrix(pdoubt_minus)));

			vector<vector<double>> p_hat_minus = multiplyMatrixWithANumber(convertTupleToMatrix(pdoubt_minus), (1 - λp))
				+ multiplyMatrixWithANumber(convertTupleToMatrix(pdoubt), λp);

			vector<vector<double>> p_hat_plus = multiplyMatrixWithANumber(convertTupleToMatrix(pdoubt_plus), (1 - λp))
				+ multiplyMatrixWithANumber(convertTupleToMatrix(pdoubt), λp);
			//cần phải chuyển vị dp về dạng ma trận hàng
			vector<vector<double>> dp_transpose = transpose(dp);
			tuple<double, double> dp_tuple = convertMatrixToTuple(dp_transpose);
			D.push_back(dp_tuple);
			int pdoubt_indexp = tim_chi_so_phan_tu_tuple(P, pdoubt);
			int pdoubt_index = tim_chi_so_phan_tu_tuple(Pdoubt, pdoubt);
			Pdoubt = deleteTuple(Pdoubt, pdoubt);
			if (allclose(convertTupleToMatrix(pdoubt), p_hat_minus) == true && allclose(convertTupleToMatrix(pdoubt), p_hat_plus) == true) {
				count4++;
				int ptest_index = tim_chi_so_phan_tu_tuple(Ptest, pdoubt);
				Ptest = dao_vi_tri_ptest(Ptest, ptest_index);
			}
			else {
				P = delete_tuple_by_index(P, pdoubt_indexp);
				int ptest_index = tim_chi_so_phan_tu_tuple(Ptest, pdoubt);
				Ptest = delete_tuple_by_index(Ptest, ptest_index);
				if (allclose(p_hat_plus, convertTupleToMatrix(pdoubt_plus))) {
					cout << "1" << endl;
				}
				else {
					P = insert(P, pdoubt_indexp, convertMatrixToTuple(p_hat_plus));
					Pdoubt = insert(Pdoubt, pdoubt_index, convertMatrixToTuple(p_hat_plus));
					Ptest = insert(Ptest, ptest_index, convertMatrixToTuple(p_hat_plus));
				}
				if (allclose(p_hat_minus, convertTupleToMatrix(pdoubt_minus))) {
					cout << "2" << endl;
				}
				else {
					P = insert(P, pdoubt_indexp, convertMatrixToTuple(p_hat_minus));
					Pdoubt = insert(Pdoubt, pdoubt_index, convertMatrixToTuple(p_hat_minus));
					Ptest = insert(Ptest, ptest_index, convertMatrixToTuple(p_hat_minus));
				}

			}
		}
		else {
			count3++;
			Pdoubt = deleteTuple(Pdoubt, pdoubt);
			vector<int> ptest_index = find_list_index(Ptest, pdoubt);
			int first_index = ptest_index.at(0);
			Ptest = dao_vi_tri_ptest(Ptest, first_index);
		}
	}
	std::vector<Point> output = to_point_vector(P);
	// Gán các giá trị của mảng points2 cho mảng points1 bắt đầu từ vị trí số 0
	// Gán các giá trị của mảng points2 cho mảng points1 bắt đầu từ vị trí số 0
	copy_vector_to_array(in_poly, output, output.size());
	n_poly = output.size();

	return P;
}


int main() {

	Point* points = new Point[50];

	// Sinh ngẫu nhiên giá trị x, y cho các phần tử trong mảng
	default_random_engine generator(time(NULL));
	uniform_real_distribution<double> distribution(-50.0, 50.0);
	for (int i = 0; i < 50; i++) {
		points[i].x = distribution(generator);
		points[i].y = distribution(generator);
	}
	int n_poly = 50;
	vector<tuple<double, double>> result = OuterConvexApproximation(points, n_poly, 0);
	// In ra giá trị x, y của các phần tử trong mảng
	cout << "======================================================" << endl;
	cout << "tap hop cac diem point:" << endl;
	for (int i = 0; i < 50; i++) {
		cout << "[" << points[i].x << ", " << points[i].y << "]" << endl;
	}
	cout << "====================================================" << endl;
	cout << "bao loi tim duoc la:";
	PrintVector(result);
	// Giải phóng bộ nhớ
	delete[] points;

	return 0;
}


