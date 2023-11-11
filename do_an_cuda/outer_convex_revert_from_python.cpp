﻿#define _USE_MATH_DEFINES

#include <iostream>


#define NMAX 512
#define NMIN -999999
using namespace std;

struct Point {
    double x, y;

    Point() {}

    Point(double x, double y) : x(x), y(y) {}
};

//copy phần tử của mảng này chuyển sang mảng khác
void copy_points(Point *src, Point *dst, int &n_src, int &n_dst) {
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
int *find_all_point(Point *P, int n, Point p) {
    // Khởi tạo mảng kết quả
    int *result = new int[n];
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
void insert_point_to_index(Point *P, int &n, int index, Point p) {
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
void delete_point_by_index(Point *P, int &n, int index) {
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
void move_point_to_end(Point *P, int n, int index) {
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
bool allclose(const Point &p1, const Point &p2, double rtol = 1e-5, double atol = 1e-8) {
    // Tính độ lệch tương đối giữa hai phần tử x
    double rel_diff_x = fabs(p1.x - p2.x) / (atol + rtol * fmax(fabs(p1.x), fabs(p2.x)));

    // Tính độ lệch tương đối giữa hai phần tử y
    double rel_diff_y = fabs(p1.y - p2.y) / (atol + rtol * fmax(fabs(p1.y), fabs(p2.y)));

    // Kiểm tra xem cả hai độ lệch tương đối đều nhỏ hơn hoặc bằng 1.0
    return rel_diff_x <= 1.0 && rel_diff_y <= 1.0;
}

//xoá 1 Point trong Point* không dùng index
void delete_point(Point *P, int &n, Point pdoubt) {
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
int find_index(Point *P, int n, Point pdoubt) {
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
double **add_two_matrix(double **A, double **B, int m, int n) {
    // Duyệt qua tất cả các phần tử của hai ma trận
    double **C = new double *[m];
    for (int i = 0; i < n; i++) {
        C[i] = new double[n];
    }
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

    // Trả về ma trận tổng
    return C;
}

//nhân ma trận double** với 1 số double
double **multiply_matrix_with_double(double **matrix, int n, int m, double scalar) {
    // Khai báo ma trận kết quả
    double **result = new double *[n];
    for (int i = 0; i < n; i++) {
        result[i] = new double[m];
    }

    // Lặp qua các phần tử của ma trận cần nhân và nhân với số double cần nhân
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            result[i][j] = matrix[i][j] * scalar;
        }
    }

    return result;
}

//xoá Point trong Point*
void deletePoint(Point *arr, int &n, Point p) {
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
Point convert_double_to_point(double **matrix) {
    // Khởi tạo điểm
    Point point = {0, 0};

    // Gán giá trị cho điểm
    point.x = matrix[0][0];
    point.y = matrix[0][1];

    // Trả về điểm
    return point;
}

//tích vo hướng của ma trận và point
double dot_product(double **matrix, Point point) {
    // Khởi tạo tích vô hướng
    double dot_product = 0;
    dot_product += matrix[0][0] * point.x + matrix[0][1] * point.y;
    return dot_product;
}

//lấy ra phần tử lớn nhất trong ma trận
double get_max_value(double **mul_dp_xtranspose, int rows, int n_poly) {
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

double **transpose_matrix(double **A, int cols, int rows) {
    // Khai báo ma trận chuyển vị
    double **B = new double *[cols];
    for (int i = 0; i < cols; i++) {
        B[i] = new double[rows];
    }

    // Duyệt qua tất cả các phần tử của ma trận A
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            // Gán giá trị của phần tử A[i][j] cho vị trí B[j][i]
            B[j][i] = A[i][j];
        }
    }

    // Trả về ma trận chuyển vị
    return B;
}

//ma trận chuyển vị của con trỏ double**, chuyển dọc sang ngang
double **transpose_dp(double **matrix) {
    // Khởi tạo ma trận chuyển vị
    double **transposed_matrix = new double *[1];
    transposed_matrix[0] = new double[2];

    // Gán giá trị cho ma trận chuyển vị
    transposed_matrix[0][0] = matrix[0][0];
    transposed_matrix[0][1] = matrix[1][0];

    // Trả về ma trận chuyển vị
    return transposed_matrix;
}


//chuyển Point sang double**, áp dụng cho chuyển đổi mảng in_poly thành X
double **convert_point_to_matrix(Point *points, int n) {
    // Khởi tạo ma trận 2 chiều
    double **matrix = new double *[n];
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

// chia ma trận cho một số trả về một ma trận mới
double **divide_matrix_by_double_and_return_new_matrix(double **matrix, int rows, int cols, double d) {
    // Tạo mảng double** mới
    double **new_matrix = new double *[rows];
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

//tính chuẩn euclid của Point
double norm_2(Point p) {
    // Tính bình phương của từng phần tử của Point
    double x2 = p.x * p.x;
    double y2 = p.y * p.y;

    // Trả về căn bậc hai của tổng bình phương
    return sqrt(x2 + y2);
}

//nhân hai ma trận với nhau
double **multiply_matrix(double **A, double **B, int m, int n, int p) {
    // Khởi tạo ma trận kết quả
    double **C = new double *[m];
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

//trừ hai điểm cho nhau
Point subtract_points(Point p1, Point p2) {
    // Hàm trừ 2 Point cho nhau

    Point p3;
    p3.x = p1.x - p2.x;
    p3.y = p1.y - p2.y;

    return p3;
}

//chuyển Point sang ma trận cột chỉ dành cho ma trận trả về kích thước 2x1
double **convert_point_to_matrix(Point p) {
    // Khởi tạo ma trận 2 hàng 1 cột

    double **matrix = new double *[2];
    matrix[0] = new double[1];
    matrix[1] = new double[1];

    // Gán giá trị cho ma trận
    matrix[0][0] = p.x;
    matrix[1][0] = p.y;

    return matrix;
}

//chuyển Point sang ma trận cột chỉ dành cho ma trận trả về kích thước 1x2
double **convert_point_to_row_matrix(Point p) {
    // Khởi tạo ma trận 2 hàng 1 cột

    double **matrix = new double *[1];
    matrix[0] = new double[2];

    // Gán giá trị cho ma trận
    matrix[0][0] = p.x;
    matrix[0][1] = p.y;

    return matrix;
}

//tìm chỉ số của point trong Point*
int find_point_index(Point *Ptest, Point pdoubt) {
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


Point *OuterConvexApproximation(Point *in_poly, int &n_poly) {
    //global
    double δ = 0.0;
    //khởi tạo mảng xoay R, các hàm sin cos ở dưới cũng dùng thư viện cmath
    double alpha = -M_PI / 2;
    //khởi tạo R
    //đã free
    double **R = new double *[2];
    for (int i = 0; i < 2; i++) {
        R[i] = new double[2];
    }
    R[0][0] = cos(alpha);
    R[0][1] = sin(alpha);
    R[1][0] = -sin(alpha);
    R[1][1] = cos(alpha);

    // Tạo một mảng các đối tượng Point
    //nhớ dùng xong phải gán lại về bằng -999 hết
    Point *D = new Point[NMAX];
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

    //khởi tạo P
    Point *P = new Point[NMAX];
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
    Point *Pdoubt = new Point[NMAX];
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
    Point *Ptest = new Point[NMAX];
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
        Point result_sub_pminus_plus = subtract_points(pdoubt_minus, pdoubt_plus);

        // Chuyển đổi Point p thành ma trận 2 hàng 1 cột
        double **transposed_matrix = convert_point_to_matrix(result_sub_pminus_plus);
        //thực hiện nhân ma trận chuyển vị trên với R
        double **result_mul_matrix = multiply_matrix(R, transposed_matrix, 2, 2, 1);

        //tính chuẩn norm 2 của kết quả này
        double norm_sub_pminus_pplus = norm_2(result_sub_pminus_plus);
        //dp ở đây có kích thước 2x1, là ma trận cột
        double **dp = divide_matrix_by_double_and_return_new_matrix(result_mul_matrix, 2, 1, norm_sub_pminus_pplus);
        //=========lỗi từ đây==============
        //chuyển X về ma trận cho dễ làm việc
        double **X = convert_point_to_matrix(in_poly, n_poly);
        double **X_transpose = transpose_matrix(X, 2, n_poly);
        //dp_tranpose 1x2
        double **dp_transpose = transpose_dp(dp);
        double **mul_dp_xtranspose = multiply_matrix(dp_transpose, X_transpose, 1, 2, n_poly);
        //tính Bdp
        double βdp = get_max_value(mul_dp_xtranspose, 1, n_poly);
        //cần tính tích vô hướng trong điều kiện dưới
        if (βdp == dot_product(dp_transpose, pdoubt_plus)) {
            count1 += 1;
            //chuyển đổi dp_transpose về kiểu Point để nạp vào tập D
            Point dp_transpose_point = convert_double_to_point(dp_transpose);
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

            double **pdoubt_minus_convert_to_matrix = convert_point_to_row_matrix(pdoubt_minus);
            double **pdoubt_convert_to_matrix = convert_point_to_row_matrix(pdoubt);
            double **A = multiply_matrix_with_double(pdoubt_minus_convert_to_matrix, 1, 2, (1 - λp));
            double **B = multiply_matrix_with_double(pdoubt_convert_to_matrix, 1, 2, λp);
            double **p_hat_minus = add_two_matrix(A, B, 1, 2);
            //            double **p_hat_minus = add_two_matrix(
            //                    multiply_matrix_with_double(convert_point_to_row_matrix(pdoubt_minus), 1, 2, (1 - λp)),
            //                    multiply_matrix_with_double(convert_point_to_row_matrix(pdoubt), 1, 2, λp), 1, 2);

            //1x2


            double **pdoubt_plus_convert_to_matrix = convert_point_to_row_matrix(pdoubt_plus);
            double **C = multiply_matrix_with_double(pdoubt_plus_convert_to_matrix, 1, 2, (1 - λp));
            double **p_hat_plus = add_two_matrix(C, B, 1, 2);
            //            double **p_hat_plus = add_two_matrix(
            //                    multiply_matrix_with_double(convert_point_to_row_matrix(pdoubt_plus), 1, 2, (1 - λp)),
            //                    multiply_matrix_with_double(convert_point_to_row_matrix(pdoubt), 1, 2, λp), 1, 2);
            Point dp_transpose_point = convert_double_to_point(dp_transpose);
            D[D_size] = dp_transpose_point;
            D_size++;

            int pdoubt_indexp = find_index(P, size_P, pdoubt);

            int pdoubt_index = find_index(Pdoubt, size_Pdoubt, pdoubt);

            delete_point(Pdoubt, size_Pdoubt, pdoubt);
            if (allclose(pdoubt, convert_double_to_point(p_hat_minus))
                && allclose(pdoubt, convert_double_to_point(p_hat_plus))) {

                count4++;
                int ptest_index = find_point_index(Ptest, pdoubt);
                move_point_to_end(Ptest, size_Ptest, ptest_index);
            } else {
                delete_point_by_index(P, size_P, pdoubt_indexp);
                int ptest_index = find_point_index(Ptest, pdoubt);
                delete_point_by_index(Ptest, size_Ptest, ptest_index);

                if (allclose(convert_double_to_point(p_hat_plus), pdoubt_plus)) {
                    //cout << "1" << endl;
                } else {
                    Point tmp = convert_double_to_point(p_hat_plus);
                    insert_point_to_index(P, size_P, pdoubt_indexp, convert_double_to_point(p_hat_plus));
                    insert_point_to_index(Pdoubt, size_Pdoubt, pdoubt_index, convert_double_to_point(p_hat_plus));
                    insert_point_to_index(Ptest, size_Ptest, ptest_index, convert_double_to_point(p_hat_plus));
                }
                if (allclose(convert_double_to_point(p_hat_minus), pdoubt_minus)) {
                    //cout << "2" << endl;
                } else {
                    insert_point_to_index(P, size_P, pdoubt_indexp, convert_double_to_point(p_hat_minus));
                    insert_point_to_index(Pdoubt, size_Pdoubt, pdoubt_index, convert_double_to_point(p_hat_minus));
                    insert_point_to_index(Ptest, size_Ptest, ptest_index, convert_double_to_point(p_hat_minus));

                }

            }
            for (int i = 0; i < 2; i++) delete p_hat_plus[i];
            for (int i = 0; i < 2; i++) delete p_hat_minus[i];
            delete[] pdoubt_minus_convert_to_matrix;
            delete[] pdoubt_plus_convert_to_matrix;
            delete[] pdoubt_convert_to_matrix;
            delete[] A;
            delete[] B;
            delete[] C;
        } else {
            count3++;
            delete_point(Pdoubt, size_Pdoubt, pdoubt);
            int *ptest_index = find_all_point(Ptest, size_Ptest, pdoubt);
            int first_index = ptest_index[0];
            move_point_to_end(Ptest, size_Ptest, first_index);
        }
        delete transposed_matrix[0];
        delete transposed_matrix[1];
        delete[] transposed_matrix;

        delete result_mul_matrix[0];
        delete result_mul_matrix[1];
        delete[] result_mul_matrix;

        delete dp[0];
        delete dp[1];
        delete[] dp;
        for (int i = 0; i < n_poly; i++) {
            delete X[i];
        }
        delete[] X;

        delete X_transpose[0];
        delete X_transpose[1];
        delete[] X_transpose;

        delete dp_transpose[0];
        delete[] dp_transpose;

        delete mul_dp_xtranspose[0];
        delete[] mul_dp_xtranspose;
    }

    copy_points(P, in_poly, size_P, n_poly);
    //cần phải giải phóng hết bộ nhớ đi
    for (int i = 0; i < 2; i++) {
        delete[] R[i];
    }
    delete[] R;
    delete[] D;
    delete[] P;
    delete[] Pdoubt;
    delete[] Ptest;
    return in_poly;
}


int main() {

    Point *points = new Point[50];

    points[0] = {-36.8462, -4.13499};
    points[1] = {-28.1041, 17.8865};
    points[2] = {43.4693, 1.94164};
    points[3] = {-46.5428, 2.97002};
    points[4] = {-49.2302, -43.3158};
    points[5] = {18.6773, 43.0436};
    points[6] = {2.69288, 15.3919};
    points[7] = {20.1191, 26.2198};
    points[8] = {-45.2535, -17.1766};
    points[9] = {25.641, -13.4661};
    points[10] = {48.255, 25.3356};
    points[11] = {-42.7314, 38.4707};
    points[12] = {-6.35886, -2.22682};
    points[13] = {-22.5093, -33.3493};
    points[14] = {39.7656, -43.9436};
    points[15] = {0.452289, -18.0967};
    points[16] = {-0.602331, -40.9267};
    points[17] = {-42.6251, -11.5858};
    points[18] = {41.3817, -3.55542};
    points[19] = {-44.9916, 27.0205};
    points[20] = {-37.4635, 18.8455};
    points[21] = {12.9543, 22.5412};
    points[22] = {38.8572, -19.3678};
    points[23] = {1.32737, 34.5982};
    points[24] = {34.1511, -8.46054};
    points[25] = {-3.20826, -32.1672};
    points[26] = {7.16548, -46.6946};
    points[27] = {-0.151988, 24.8293};
    points[28] = {39.0737, 34.204};
    points[29] = {-28.7248, -36.9573};
    points[30] = {-22.5412, -8.57067};
    points[31] = {20.982, -26.0089};
    points[32] = {-18.246, 15.2059};
    points[33] = {18.1346, -11.2275};
    points[34] = {-35.2467, 34.5576};
    points[35] = {45.5409, -35.1848};
    points[36] = {-9.12333, 6.48987};
    points[37] = {-1.14855, 46.1095};
    points[38] = {-30.0243, 12.9269};
    points[39] = {15.1254, 30.3073};
    points[40] = {-2.35682, -29.675};
    points[41] = {40.1673, -35.7979};
    points[42] = {-8.9687, 38.5648};
    points[43] = {-33.7801, -13.4661};
    points[44] = {-36.4891, -4.46927};
    points[45] = {-4.76998, 43.1674};
    points[46] = {-28.4752, 40.8922};
    points[47] = {36.086, 0.595588};
    points[48] = {31.7561, -3.7755};
    points[49] = {13.2739, 32.4697};


    int n_poly = 50;
    cout << "======================================================" << endl;
    cout << "tap hop cac diem point:" << endl;
    for (int i = 0; i < n_poly; i++) {
        cout << "[" << points[i].x << ", " << points[i].y << "]" << endl;
    }
    cout << "====================================================" << endl;

    Point *result = OuterConvexApproximation(points, n_poly);
    // In ra giá trị x, y của các phần tử trong mảng
    cout << "======================================================" << endl;
    cout << "tap hop bao loi:" << n_poly << " diem" << endl;
    for (int i = 0; i < n_poly; i++) {
        cout << "[" << result[i].x << ", " << result[i].y << "]" << endl;
    }
    cout << "====================================================" << endl;

    // Giải phóng bộ nhớ
    delete[] points;

    return 0;
}


