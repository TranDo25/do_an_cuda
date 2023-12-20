//
// Created by ADMIN on 18/12/10023.
//
#include<iostream>
#include<cmath>

using namespace std;
const double EPS = 1E-8;
const double PI = 3.14159265358979323846;

struct Point {
    double x, y;

    Point() {}

    Point(double x, double y) : x(x), y(y) {}
};
int sig(double d) { return (d > EPS) - (d < -EPS); }
bool point_same(Point& a, Point& b) {
    return sig(a.x - b.x) == 0 && sig(a.y - b.y) == 0;
}
//=============================
struct Edge {
    Point first_point;
    Point second_point;
};
//thêm hai hàm mới bên dưới này vào nữa
int findMinPointIndex(const Point points[], int n) {
    if (n <= 0) {
        return -1;
    }

    int minIndex = 0;

    for (int i = 1; i < n; ++i) {
        if (points[i].y < points[minIndex].y) {
            minIndex = i;
        }
        else if (points[i].y == points[minIndex].y && points[i].x < points[minIndex].x) {
            minIndex = i;
        }
    }

    return minIndex;
}

void moveMinPointsToFront(Point points[], int n) {
    int minIndex = findMinPointIndex(points, n);
    Point points_premitive[20];
    for (int i = 0; i < n; i++) {
        points_premitive[i] = points[i];
    }
    if (minIndex != -1) {
        Point temp[20];
        int tempIndex = 0;

        // Lưu trữ các phần tử từ minIndex đến cuối mảng vào mảng tạm thời.
        for (int i = minIndex; i < n; ++i) {
            temp[tempIndex++] = points[i];
        }

        // Di chuyển các phần tử từ 0 đến minIndex - 1 lên một vị trí.
        for (int i = 0; i < minIndex; ++i) {
            points[i + n - minIndex] = points_premitive[i];
        }

        // Đặt các phần tử có hoành độ và tung độ nhỏ nhất lên đầu mảng.
        for (int i = 0; i < tempIndex; ++i) {
            points[i] = temp[i];
        }
    }
}
// Hàm đảo ngược các phần tử trong mảng Point
void reversePoints(Point* points, int size) {
    int start = 0;
    int end = size - 1;

    while (start < end) {
        // Swap giữa points[start] và points[end]
        Point temp = points[start];
        points[start] = points[end];
        points[end] = temp;

        // Di chuyển lần lượt đến phần tử thứ hai và thứ hai từ cuối
        ++start;
        --end;
    }
}
//copy phần tử của mảng này chuyển sang mảng khác
void copy_points(Point* src, Point* dst, int& n_src, int& n_dst) {
    // Sao chép từng phần tử của mảng
    for (int i = 0; i < n_src; i++) {
        dst[i].x = src[i].x;
        dst[i].y = src[i].y;

    }
    //thay đổi lại số phần tử của mảng n_poly
    n_dst = n_src;
}
// Hàm kiểm tra xem một Point đã tồn tại trong mảng hay chưa
bool pointExists(const Point* uniquePoints, int size, const Point& p) {
    for (int i = 0; i < size; ++i) {
        if (uniquePoints[i].x == p.x && uniquePoints[i].y == p.y) {
            return true;
        }
    }
    return false;
}

// Hàm lọc ra các Point phân biệt từ mảng các Edge
void uniquePoints(const Edge* edges, int size, Point* uniquePoints, int& uniqueSize) {
    for (int i = 0; i < size; ++i) {
        // Kiểm tra first_point
        if (!pointExists(uniquePoints, uniqueSize, edges[i].first_point)) {
            uniquePoints[uniqueSize] = edges[i].first_point;
            ++uniqueSize;
        }

        // Kiểm tra second_point
        if (!pointExists(uniquePoints, uniqueSize, edges[i].second_point)) {
            uniquePoints[uniqueSize] = edges[i].second_point;
            ++uniqueSize;
        }
    }
}

//chen 1 canh vao phan tu
void insert_edge(Edge* edges, int& n, const Edge& edge, int index) {

    if (index == -1) {
        index = n;
    }
    for (int i = n - 1; i >= index; i--) {
        edges[i + 1] = edges[i];
    }

    edges[index] = edge;
    n++;
}

//xoa canh theo chi so
void erase_edge(Edge* edges, int& n, int index) {
    if (index < 0 || index >= n) {
        return;
    }

    for (int i = index + 1; i < n; i++) {
        edges[i - 1] = edges[i];
    }

    n--;
}

//tim chi so cua 1 canh trong mang
int find_edge_index(Edge* edges, int n, Edge& edge) {
    int index = -1;
    for (int i = 0; i < n; i++) {
        if (
            edges[i].first_point.x == edge.first_point.x &&
            edges[i].first_point.y == edge.first_point.y &&
            edges[i].second_point.x == edge.second_point.x &&
            edges[i].second_point.y == edge.second_point.y
            ) {
            index = i;
            break;
        }
    }
    return index;
}

double f_max(double x, double y) {
    if (x > y) {
        return x;
    }
    else {
        return y;
    }
}

// chia ma trận cho một số trả về một ma trận mới
void divide_colmatrix_by_double(double matrix[2][1], double d) {

    // Chia từng phần tử của mảng
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 1; j++) {
            matrix[i][j] = matrix[i][j] / d;
        }
    }
}

//nhân hai ma trận với nhau
void multiply_R_matrix(double R[2][2], double B[2][1], double C[2][1], int m = 2, int n = 2, int p = 1) {
    // Khởi tạo ma trận kết quả

    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            C[i][j] = 0;
        }
    }

    // Nhân hai ma trận
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < p; j++) {
            for (int k = 0; k < n; k++) {
                C[i][j] += R[i][k] * B[k][j];
            }
        }
    }
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

// Hàm tìm giá trị lớn nhất trong mảng
double findMaxValue(double arr[], int size) {
    if (size <= 0) {
        return 0.0; // Trả về giá trị mặc định hoặc có thể điều chỉnh tùy theo yêu cầu
    }

    // Tìm giá trị lớn nhất
    double maxVal = arr[0];
    for (int i = 1; i < size; ++i) {
        if (arr[i] > maxVal) {
            maxVal = arr[i];
        }
    }

    return maxVal;
}
// Hàm xoá tất cả các cạnh có giá trị bằng edoubt khỏi mảng các cạnh
void delete_edges(Edge* edges, int& n, Edge edoubt) {
    // Khởi tạo biến đếm để lưu vị trí của các cạnh cần xoá
    int count = 0;

    // Duyệt qua tất cả các cạnh trong mảng edges
    for (int i = 0; i < n; i++) {
        // So sánh cạnh hiện tại với cạnh cần xoá
        if (
            (edges[i].first_point.x == edoubt.first_point.x &&
                edges[i].first_point.y == edoubt.first_point.y &&
                edges[i].second_point.x == edoubt.second_point.x &&
                edges[i].second_point.y == edoubt.second_point.y) ||
            (edges[i].first_point.x == edoubt.second_point.x &&
                edges[i].first_point.y == edoubt.second_point.y &&
                edges[i].second_point.x == edoubt.first_point.x &&
                edges[i].second_point.y == edoubt.first_point.y)
            ) {
            // Không sao chép cạnh này vào mảng mới
            // Điều này giống như là loại bỏ cạnh này khỏi mảng
            count++;
        }
        else {
            // Nếu không phải là cạnh cần xoá, sao chép nó vào vị trí mới
            edges[i - count] = edges[i];
        }
    }

    // Giảm kích thước của mảng edges
    n -= count;
}

void convert_point_to_colmatrix(Point p, double matrix[2][1]) {
    // Gán giá trị cho ma trận
    matrix[0][0] = p.x;
    matrix[1][0] = p.y;
}

//ma trận chuyển vị của con trỏ double**, chuyển ngang sang dọc
void transpose_to_colmatrix(double matrix[1][2], double transposed_matrix[2][1]) {
    // Gán giá trị cho ma trận chuyển vị
    transposed_matrix[0][0] = matrix[0][0];
    transposed_matrix[1][0] = matrix[0][1];
}

//ma trận chuyển vị của con trỏ double**, chuyển dọc sang ngang
void transpose_to_rowmatrix(double matrix[2][1], double transposed_matrix[1][2]) {
    // Gán giá trị cho ma trận chuyển vị
    transposed_matrix[0][0] = matrix[0][0];
    transposed_matrix[0][1] = matrix[1][0];
}

//tích vo hướng của ma trận và point
double dot_product(double matrix1[1][2], double matrix2[2][1]) {
    // Khởi tạo tích vô hướng
    double dot_product = 0;
    dot_product += matrix1[0][0] * matrix2[0][0] + matrix1[0][1] * matrix2[1][0];
    return dot_product;
}

//chuyển Point sang double**, áp dụng cho chuyển đổi mảng in_poly thành X
void convert_many_points_to_matrix(Point* points, int n, double X[][2]) {
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

//tính chuẩn euclid của Point
double norm_2(Point p1, Point p2) {
    // Tính bình phương của từng phần tử của Point
    Point tmp = { 0, 0 };
    tmp.x = p1.x - p2.x;
    tmp.y = p1.y - p2.y;
    double x2 = tmp.x * tmp.x;
    double y2 = tmp.y * tmp.y;

    // Trả về căn bậc hai của tổng bình phương
    return sqrt(x2 + y2);
}


void convert_point_to_rowmatrix(Point p, double matrix[1][2]) {
    // Gán giá trị cho ma trận
    matrix[0][0] = p.x;
    matrix[0][1] = p.y;
}

//trừ hai điểm cho nhau
void subtract_points(Point p1, Point p2, Point& p3) {
    // Hàm trừ 2 Point cho nhau
    p3.x = p1.x - p2.x;
    p3.y = p1.y - p2.y;
}

double find_max_x(Point X[], int n) {
    double max_x = X[0].x;

    for (int i = 1; i < n; i++) {
        if (X[i].x > max_x) {
            max_x = X[i].x;
        }
    }

    // Trả về giá trị lớn nhất
    return max_x;
}

double find_min_x(Point X[], int n) {
    double min_x = X[0].x;
    for (int i = 1; i < n; i++) {
        if (X[i].x < min_x) {
            min_x = X[i].x;
        }
    }

    return min_x;
}

double find_max_y(Point X[], int n) {
    double max_y = X[0].y;
    for (int i = 1; i < n; i++) {
        if (X[i].y > max_y) {
            max_y = X[i].y;
        }
    }

    // Trả về giá trị lớn nhất
    return max_y;
}

double find_min_y(Point X[], int n) {
    double min_y = X[0].y;
    for (int i = 1; i < n; i++) {
        if (X[i].y < min_y) {
            min_y = X[i].y;
        }
    }
    return min_y;
}

// Hàm tìm các phần tử Point có giá trị hoành độ bằng x
void find_points_by_x(Point* X, int& n, double min_x, Point points[], int& size_point) {

    // Duyệt qua tất cả các phần tử trong mảng X
    int count = 0;
    for (int i = 0; i < n; i++) {
        // Kiểm tra xem giá trị hoành độ của phần tử hiện tại có bằng min_x hay không
        if (X[i].x == min_x) {
            // Thêm phần tử hiện tại vào mảng
            points[count].x = X[i].x;
            points[count].y = X[i].y;
            count++;
        }
    }
    size_point = count;
}

// Hàm tìm các phần tử Point có giá trị hoành độ bằng y
void find_points_by_y(Point* X, int n, double min_y, Point points[], int& size_point) {

    // Duyệt qua tất cả các phần tử trong mảng X
    int count = 0;
    for (int i = 0; i < n; i++) {
        // Kiểm tra xem giá trị hoành độ của phần tử hiện tại có bằng min_x hay không
        if (X[i].y == min_y) {
            // Thêm phần tử hiện tại vào mảng
            points[count].x = X[i].x;
            points[count].y = X[i].y;
            count++;
        }
    }
    size_point = count;
}

void inner_convex_approximation_and_index(Point* in_poly, int& n_poly, int* points_to_convex_ind) {
    int n_input = n_poly;
    Point input_poly[20];
    for (int i = 0; i < n_input; i++) {
        input_poly[i].x = in_poly[i].x;
        input_poly[i].y = in_poly[i].y;
    }
    //=========================hay thay the loi Jarvis vao duoi===========
    //global
    double delta = 0.0;
    //khởi tạo mảng xoay R, các hàm sin cos ở dưới cũng dùng thư viện cmath
    double alpha = -PI / 2;
    //khởi tạo R
    double R[2][2];

    R[0][0] = cos(alpha);
    R[0][1] = sin(alpha);
    R[1][0] = -sin(alpha);
    R[1][1] = cos(alpha);

    // Khởi tạo giá trị ban đầu
    double min_x = in_poly[0].x;
    double max_x = in_poly[0].x;
    double min_y = in_poly[0].y;
    double max_y = in_poly[0].y;

    // Duyệt qua mảng
    for (int i = 1; i < n_poly; i++) {
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
    Point min_X[100];
    Point max_X[100];
    Point min_Y[100];
    Point max_Y[100];
    int size_min_X = 0;
    int size_max_X = 0;
    int size_min_Y = 0;
    int size_max_Y = 0;

    //tim cac diem co hoanh do nho nhat
    find_points_by_x(in_poly, n_poly, min_x, min_X, size_min_X);
    //tim cac diem co hoanh do lon nhat
    find_points_by_x(in_poly, n_poly, max_x, max_X, size_max_X);
    //tim cac diem co tung do nho nhat
    find_points_by_y(in_poly, n_poly, min_y, min_Y, size_min_Y);
    //tim cac diem co tung do lon nhat
    find_points_by_y(in_poly, n_poly, max_y, max_Y, size_max_Y);

    double q12 = find_max_y(max_X, size_max_X);
    double q21 = find_min_x(max_Y, size_max_Y);
    double q32 = find_min_y(min_X, size_min_X);
    double q41 = find_max_x(min_Y, size_min_Y);

    Point q1 = { max_x, q12 };
    Point q2 = { q21, max_y };
    Point q3 = { min_x, q32 };
    Point q4 = { q41, min_y };

    Point Xcomma[100];
    Xcomma[0] = q1;
    Xcomma[1] = q2;
    Xcomma[2] = q3;
    Xcomma[3] = q4;
    int size_X_comma = 4;

    Edge edge0 = { q1, q2 };
    Edge edge1 = { q2, q3 };
    Edge edge2 = { q3, q4 };
    Edge edge3 = { q4, q1 };
    //khoi tao E_test
    Edge E_test[100];
    E_test[0] = edge0;
    E_test[1] = edge1;
    E_test[2] = edge2;
    E_test[3] = edge3;
    int size_E_test = 4;

    //khoi tao Edoubt
    Edge Edoubt[100];
    for (int i = 0; i < size_E_test; i++) {
        Edoubt[i] = E_test[i];
    }
    int size_Edoubt = size_E_test;
    while (size_Edoubt > 0) {
        Edge edoubt = Edoubt[0];
        Point p = edoubt.first_point;
        Point p_plus = edoubt.second_point;

        //tinh dpp_plus theo cong thuc (44)
        Point result_sub_pplus_p = { 0, 0 };
        subtract_points(p_plus, p, result_sub_pplus_p);
        double transposed_matrix[2][1];
        convert_point_to_colmatrix(result_sub_pplus_p, transposed_matrix);
        double result_mul_matrix[2][1];
        multiply_R_matrix(R, transposed_matrix, result_mul_matrix);
        double norm_sub_pplus_p = norm_2(result_sub_pplus_p);
        double dpp_plus[2][1];
        divide_matrix_by_double_and_return_new_matrix(result_mul_matrix, dpp_plus, 2, 1, norm_sub_pplus_p);
        //chuyển in_poly về ma trận cho dễ làm việc
        double X[100][2];
        convert_many_points_to_matrix(in_poly, n_poly, X);

        Point Xp_pplus[100];
        int size_Xpp_plus = 0;
        double dpp_plus_T[1][2];
        transpose_to_rowmatrix(dpp_plus, dpp_plus_T);
        double p_T[2][1];
        convert_point_to_colmatrix(p, p_T);
        double dpp_plus_mul_p_T = dot_product(dpp_plus_T, p_T);
        for (int i = 0; i < n_poly; i++) {
            double item_X_T[2][1];
            item_X_T[0][0] = X[i][0];
            item_X_T[1][0] = X[i][1];

            double tmp = dot_product(dpp_plus_T, item_X_T);
            if (tmp > dpp_plus_mul_p_T) {
                Xp_pplus[size_Xpp_plus].x = X[i][0];
                Xp_pplus[size_Xpp_plus].y = X[i][1];
                size_Xpp_plus++;
            }
        }
        if (size_Xpp_plus != 0) {
            double X_p_pplus0[1][2];
            X_p_pplus0[0][0] = Xp_pplus[0].x;
            X_p_pplus0[0][1] = Xp_pplus[0].y;
            double X_p_pplus0_T[2][1];
            transpose_to_colmatrix(X_p_pplus0, X_p_pplus0_T);
            double beta_p_pplus = dot_product(dpp_plus_T, X_p_pplus0_T);
            for (int i = 0; i < size_Xpp_plus; i++) {
                double item_X_p_pplus[1][2];
                item_X_p_pplus[0][0] = Xp_pplus[i].x;
                item_X_p_pplus[0][1] = Xp_pplus[i].y;
                double item_X_p_pplus_T[2][1];
                transpose_to_colmatrix(item_X_p_pplus, item_X_p_pplus_T);
                double tmp = dot_product(dpp_plus_T, item_X_p_pplus_T);
                if (tmp > beta_p_pplus) {
                    beta_p_pplus = tmp;
                }
            }

            Point Bpp_plus[100];
            int size_Bpp_plus = 0;
            for (int i = 0; i < size_Xpp_plus; i++) {
                double item_Xpp_plus[1][2];
                item_Xpp_plus[0][0] = Xp_pplus[i].x;
                item_Xpp_plus[0][1] = Xp_pplus[i].y;
                double item_Xpp_plus_T[2][1];
                transpose_to_colmatrix(item_Xpp_plus, item_Xpp_plus_T);
                if (dot_product(dpp_plus_T, item_Xpp_plus_T) == beta_p_pplus) {
                    Bpp_plus[size_Bpp_plus].x = Xp_pplus[i].x;
                    Bpp_plus[size_Bpp_plus].y = Xp_pplus[i].y;
                    size_Bpp_plus++;
                }

            }
            if (beta_p_pplus - dot_product(dpp_plus_T, p_T) <= delta) {
                delete_edges(Edoubt, size_Edoubt, edoubt);
            }
            else {
                Point p_hat[100];
                int size_p_hat = 0;
                Point p_hat_calc[100];
                int size_p_hat_calc = 0;
                for (int i = 0; i < size_Bpp_plus; i++) {
                    p_hat_calc[size_p_hat_calc].x = Bpp_plus[i].x;
                    p_hat_calc[size_p_hat_calc].y = Bpp_plus[i].y;
                    size_p_hat_calc++;
                }
                Point X_calc[100];
                int size_X_calc = 0;
                for (int i = 0; i < size_Bpp_plus; i++) {
                    X_calc[size_X_calc].x = Bpp_plus[i].x;
                    X_calc[size_X_calc].y = Bpp_plus[i].y;
                    size_X_calc++;
                }
                Point p_hat_calc_sub_p = { 0, 0 };
                subtract_points(p_hat_calc[0], p, p_hat_calc_sub_p);
                //cần tính max của chuẩn || x-p ||
                //mang nay dung de tim max cua norm(x_calc - p)
                double arr_norm_xcalc_p[100];

                int size_arr_norm_xcalc_p = 0;
                for (int i = 0; i < size_Bpp_plus; i++) {
                    Point x_calc_sub_p = { 0, 0 };
                    subtract_points(X_calc[i], p, x_calc_sub_p);
                    double res_norm = norm_2(x_calc_sub_p);
                    arr_norm_xcalc_p[size_arr_norm_xcalc_p++] = res_norm;
                }
                double max_arr_norm_xcalc_p = findMaxValue(arr_norm_xcalc_p, size_arr_norm_xcalc_p);
                for (int i = 0; i < size_Bpp_plus; i++) {
                    Point tmp = { 0, 0 };
                    subtract_points(Bpp_plus[i], p, tmp);
                    if (norm_2(tmp) == max_arr_norm_xcalc_p) {
                        p_hat[size_p_hat++] = Bpp_plus[i];
                    }
                }
                Point p_hat_point = p_hat[0];
                if (allclose(p_hat_point, p) && !allclose(p_hat_point, p_plus)) {
                    int edoubt_index = find_edge_index(Edoubt, size_Edoubt, edoubt);
                    erase_edge(Edoubt, size_Edoubt, edoubt_index);
                }
                else if (allclose(p_hat_point, p_plus) && !allclose(p_hat_point, p)) {
                    int edoubt_index = find_edge_index(Edoubt, size_Edoubt, edoubt);
                    erase_edge(Edoubt, size_Edoubt, edoubt_index);
                }
                else if (allclose(p_hat_point, p_plus) && allclose(p_hat_point, p)) {
                    //nothing to do
                }
                else {
                    // tinh d[p, p^]
                    double dpp_hat[2][1];
                    Point tmp = { 0, 0 };
                    subtract_points(p_hat_point, p, tmp);
                    double tmp2[1][2];
                    convert_point_to_rowmatrix(tmp, tmp2);
                    double tmp3[2][1];
                    transpose_to_colmatrix(tmp2, tmp3);
                    multiply_R_matrix(R, tmp3, dpp_hat);
                    divide_colmatrix_by_double(dpp_hat, norm_2(p_hat_point, p));

                    //tinh d[p^, p+]
                    double dp_hat_p_plus[2][1];
                    Point tmp7 = { 0, 0 };
                    subtract_points(p_plus, p_hat_point, tmp7);
                    double tmp8[1][2];
                    convert_point_to_rowmatrix(tmp7, tmp8);
                    double tmp9[2][1];
                    transpose_to_colmatrix(tmp8, tmp9);
                    multiply_R_matrix(R, tmp9, dp_hat_p_plus);
                    divide_colmatrix_by_double(dp_hat_p_plus, norm_2(p_plus, p_hat_point));

                    //tinh d[p^, p.T]
                    double dp_hat_p_T[2][1];
                    Point tmp4 = { 0, 0 };
                    subtract_points(p, p_hat_point, tmp4);
                    double tmp5[1][2];
                    convert_point_to_rowmatrix(tmp4, tmp5);
                    double tmp6[2][1];
                    transpose_to_colmatrix(tmp5, tmp6);
                    multiply_R_matrix(R, tmp6, dp_hat_p_T);
                    divide_colmatrix_by_double(dp_hat_p_T, norm_2(p, p_hat_point));


                    //tim tap X[p, p^] va X[p^, p+]. Tuy nhien lai thay trong code python khong dung
                    //den hai tap nay, nen khong can thuc thi cho mat thoi gian

                    //them p_hat_point vao XComma
                    Xcomma[size_X_comma++] = p_hat_point;

                    Edge A = { p, p_hat_point };
                    Edge B = { p_hat_point, p_plus };
                    //xoa canh cu, them 2 canh moi vao mang E_test
                    int e_index = find_edge_index(E_test, size_E_test, edoubt);
                    erase_edge(E_test, size_E_test, e_index);
                    insert_edge(E_test, size_E_test, A, e_index);
                    insert_edge(E_test, size_E_test, B, e_index + 1);
                    //xoa canh cu, them 2 canh moi vao mang E_test
                    int edoubt_index = find_edge_index(Edoubt, size_Edoubt, edoubt);
                    erase_edge(Edoubt, size_Edoubt, edoubt_index);
                    insert_edge(Edoubt, size_Edoubt, A, edoubt_index);
                    insert_edge(Edoubt, size_Edoubt, B, edoubt_index + 1);

                }
            }


        }
        else {
            delete_edges(Edoubt, size_Edoubt, edoubt);
        }
    }

    int uniqueSize = 0;
    Point in_uniquePoints[20];
    uniquePoints(E_test, size_E_test, in_uniquePoints, uniqueSize);
    reversePoints(in_uniquePoints, uniqueSize);
    moveMinPointsToFront(in_uniquePoints, uniqueSize);
    copy_points(in_uniquePoints, in_poly, uniqueSize, n_poly);

    //================================================================

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
        Point points[100] = {
                //cac diem co hoanh do nho nhat
                Point(-3, 1),
                Point(-3, -1),
                // cac diem co hoanh do lon nhat
                Point(4, -1),
                Point(4, 1),
    
                //cac diem co tung do nho nhat
                Point(-2, -1),
    
                //cac diem co tung do lon nhat
                Point(1, 4),
                Point(1, 2),
                Point(2, 1),
                Point(3, 2),
                Point(-1, 1),
                Point(-2, 3),
    
        };
    //    Point points[100] = {
    //            Point( 12.71359434 , -21.73167608),
    //            Point(-20.38584267 , -57.00750014),
    //            Point(-75.60994969 ,  40.07504469),
    //            Point(-51.28840354 , -78.34603557),
    //            Point( 31.01026119 ,   9.79174368),
    //            Point(-93.54753399 ,   4.56224875),
    //            Point(-80.78339439 ,  95.87260967),
    //            Point( 90.61767492 ,  -1.93642264),
    //            Point( 88.71974025 ,-141.23911972),
    //            Point(  1.34299907 ,-146.06373453),
    //            Point(-18.40068414 ,-135.56561234),
    //            Point(-95.44704724 ,  32.50642382),
    //            Point(-71.60161099 , -50.45836824),
    //            Point( -8.24640395 ,-100.21427597),
    //            Point( 35.94412166 ,  -5.68853525),
    //            Point(  1.34060032 , -38.00406922),
    //            Point( 75.17118834 ,  26.8814041 ),
    //            Point(-86.45155878 ,  67.45278846),
    //    };
    //Point points[100] = {
    //        Point(-1.08203696, -123.03327827),
    //        Point(-20.86190579,  -84.63800741),
    //        Point(-85.30525529,  127.19585274),
    //        Point(50.9714922 ,   24.44729803),
    //        Point(101.88982227,   56.0780549),
    //        Point(12.24647246,   68.01561449),
    //        Point(-29.06184458,   65.33139573),
    //        Point(83.16683891,  -11.53539228),
    //        Point(123.73996755, -111.00221745),
    //        Point(5.39768613,  149.6938521),
    //        Point(38.28235857,  121.2395262),
    //        Point(120.5365791 , -105.65799319),
    //        Point(-65.76238712,  125.46268328),
    //        Point(-39.19309983,  -12.2523361),
    //        Point(28.29766166,   -1.24934286),
    //        Point(-5.63865253,  -21.65336384),
    //        Point(94.82452389, -105.04749459),
    //        Point(48.04456786,   68.3858024),
    //};
    int n_poly = 11;
    for (int i = 0; i < n_poly; i++) {
        std::cout << "[" << points[i].x << ", " << points[i].y << "]\n";
    }
    int point_to_convex_indx[20] = { -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1 };
    inner_convex_approximation_and_index(points, n_poly, point_to_convex_indx);

    std::cout << "==============Inner Convex Approximation============" << std::endl;
    std::cout << "tap hop bao loi:" << n_poly << " diem" << std::endl;
    for (int i = 0; i < n_poly; i++) {
        std::cout << "[" << points[i].x << ", " << points[i].y << "]" << std::endl;
    }
    std::cout << "====================================================" << std::endl;

    std::cout << "chi so cac phan tu trong mang la:" << std::endl;


    for (int i = 0; i < 20; i++) {
        std::cout << point_to_convex_indx[i] << " ";
    }


}