﻿#define _USE_MATH_DEFINES
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
#define MAXN 100

const double EPS = 1E-8;
using namespace std;
struct Point {
	double x, y;
	Point() {}
	Point(double x, double y) : x(x), y(y) {}
};
int sig(double d) { return (d > EPS) - (d < -EPS); }
bool point_same(Point& a, Point& b) {
	return sig(a.x - b.x) == 0 && sig(a.y - b.y) == 0;
}


//================phần bên dưới dùng để đưa vào thư viện==================
void copy_vector_to_array(Point* dest, const vector<Point>& src, int size) {
	for (int i = 0; i < size; i++) {
		dest[i].x = src[i].x;
		dest[i].y = src[i].y;
	}
}

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

vector<tuple<double, double>> OuterConvexApproximation_and_index(Point* in_poly, int& n_poly, int* points_to_convex_ind) {
	int n_input = n_poly;
	cout << "n_input  = " << n_input;
	Point* input_poly = new Point[51];
	for (int i = 0; i < n_input; i++) {
		input_poly[i].x = i;
		input_poly[i].y = i;
	}
	for (int i = 0; i < n_input; i++) {
		input_poly[i].x = in_poly[i].x;
		input_poly[i].y = in_poly[i].y;
	}
	double δ = 0.0;
	
	for (int i = 0; i < n_input; i++) {
		cout << "(" << input_poly[i].x << ", " << input_poly[i].y << ") " << endl;
	}
	//khởi tạo mảng xoay R, các hàm sin cos ở dưới cũng dùng thư viện cmath
	vector<vector<double>> R(2, vector<double>(2));
	double alpha = -M_PI / 2;
	// Gán giá trị cho vector
	R[0][0] = cos(alpha);
	R[0][1] = sin(alpha);
	R[1][0] = -sin(alpha);
	R[1][1] = cos(alpha);


	vector<tuple<double, double>> D;

	D.emplace_back(1.0, 0.0);
	D.emplace_back(0.0, 1.0);
	D.emplace_back(-1.0, 0.0);
	D.emplace_back(0.0, -1.0);


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
	tuple<double, double> r1 = { max_x, max_y };
	tuple<double, double> r2 = { min_x, max_y };
	tuple<double, double> r3 = { min_x, min_y };
	tuple<double, double> r4 = { max_x, min_y };

	vector<tuple<double, double>> P = { r1, r2, r3, r4 };
	vector<tuple<double, double>> Pdoubt = P;
	vector<tuple<double, double>> Ptest = P;
	int count = 0, count4 = 0, count1 = 0, count2 = 0, count3 = 0;
	while (Pdoubt.size() > 0) {
		count += 1;
		tuple<double, double> pdoubt = Pdoubt[0];
		auto it = find(Ptest.begin(), Ptest.end(), pdoubt);
		//lấy được chỉ số của pdoubt trong Ptest
		int pdoubt_index_idx = it - Ptest.begin();
		// Lấy chỉ số của mảng liền trước
		int pdoubt_minus_index = (pdoubt_index_idx + Ptest.size() - 1) % Ptest.size();

		// Lấy chỉ số của mảng liền sau
		int pdoubt_plus_index = (pdoubt_index_idx + 1) % Ptest.size();

		// Lấy phần tử tại vị trí liền trước
		tuple<double, double> pdoubt_minus = Ptest[pdoubt_minus_index];

		// Lấy phần tử tại vị trí liền sau
		tuple<double, double> pdoubt_plus = Ptest[pdoubt_plus_index];
		//bắt đầu tính dp
		tuple<double, double> result_sub_pminus_plus = subtract_tuples(pdoubt_minus, pdoubt_plus);
		vector<tuple<double, double>> vector_result_sub_pminus_plus;
		vector_result_sub_pminus_plus.push_back(result_sub_pminus_plus);
		vector<vector<double>> transposed_matrix = transposeTupleVector(vector_result_sub_pminus_plus);
		//vector transpose này là mảng 1 chiều có kiểu là (a, b), đem nhân với ma trận R có chiều (c,d)
		 //                                                                                         (e,f)
		 //cần đổi R sang kiểu vector 2 chiều vector<vector<double>>
		std::vector<std::vector<double>> result_mul_matrix = multiplyMatrices(R, transposed_matrix);
		double norm_sub_pminus_pplus = norm(result_sub_pminus_plus);
		std::vector<std::vector<double>> dp = divide_matrix(result_mul_matrix, norm_sub_pminus_pplus);
		//chuyển in_poly về vector X cho dễ tính
		vector<tuple<double, double>> n_poly_vector = GetEnoughPoints(in_poly, n_poly);

		vector<vector<double>> X_transpose = transposeTupleVector(n_poly_vector);
		//dp đang là ma trận 2x1, vậy nếu muốn nhân với ma trận có kích thước 2x50 thì phải chuyển vị ma trận dp đi
		vector<vector<double>> dp_transpose = transpose(dp);
		vector<vector<double>> mul_dp_xtranspose = multiplyMatrices(dp_transpose, X_transpose);
		double βdp = getMax(mul_dp_xtranspose);

		//convert pdoubt_plus sang vector<vector<double>>  để có thể tái sử dụng hàm nhân ma trận
		vector<vector<double>>  pdoubt_plus_vector = convertTupleToMatrix(pdoubt_plus);

		//cần tính tích vô hướng trong điều kiện dưới
		if (βdp == dotProduct(dp_transpose, pdoubt_plus_vector)) {
			count1 += 1;
			//convert dp to tuple để thêm vào tập D
			vector<vector<double>> dp_transpose = transpose(dp);
			tuple<double, double> dp_tuple = convertMatrixToTuple(dp_transpose);
			D.push_back(dp_tuple);
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

	//========================
	for (int i = 0; i < n_poly; i++) {
		for (int j = 0; j < n_input; j++) {
			if (point_same(in_poly[i], input_poly[j])) {
				points_to_convex_ind[i] = j;
				break;
			}
		}
	}
	return P;
}


int main() {

	Point* points = new Point[9];

	// Sinh ngẫu nhiên giá trị x, y cho các phần tử trong mảng
	default_random_engine generator(time(NULL));
	uniform_real_distribution<double> distribution(-50.0, 50.0);
	for (int i = 0; i < 9; i++) {
		points[i].x = distribution(generator);
		points[i].y = distribution(generator);
	}
	int points_to_convex_ind[9] = { -1, -1, -1, -1, -1, -1, -1, -1, -1 };
	int n_poly = 9;
	vector<tuple<double, double>> result = OuterConvexApproximation_and_index(points, n_poly, points_to_convex_ind);
	// In ra giá trị x, y của các phần tử trong mảng
	cout << "======================================================" << endl;
	//cout << "tap hop cac diem point:" << endl;
	//for (int i = 0; i < 100; i++) {
	//	cout << "[" << points[i].x << ", " << points[i].y << "]" << endl;
	//}
	cout << "====================================================" << endl;
	cout << "bao loi tim duoc la:";
	PrintVector(result);
	// Giải phóng bộ nhớ
	cout << "====================================================" << endl;

	cout << "chi so cua cac diem convex so voi tap goc:" << endl;
	for (int i = 0; i < 9; i++) cout << points_to_convex_ind[i] << " ";
	delete[] points;

	return 0;
}


