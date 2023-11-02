#include <vector>
#include<tuple>
#include <iostream>
using namespace std;

// Hàm xóa 1 phần tử tuple theo index trong 1 vector<tuple<double, double>>
vector<tuple<double, double>> xóa_tuple(vector<tuple<double, double>>& v, int index) {
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

// Ví dụ sử dụng
int main() {
    // Khởi tạo vector
    vector<tuple<double, double>> v = { {1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0} };

    // Xác định chỉ số của phần tử cần xóa
    int index = 1;

    // Xóa phần tử có chỉ số là index
    v = xóa_tuple(v, index);

    // In ra vector sau khi xóa
    for (auto tuple : v) {
        cout << get<0>(tuple) << " " << get<1>(tuple) << endl;
    }

    return 0;
}
