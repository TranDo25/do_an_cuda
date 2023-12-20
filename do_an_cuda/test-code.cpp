#include <iostream>

struct Point {
    int x;
    int y;
};

struct Edge {
    Point first_point;
    Point second_point;
};

// Hàm kiểm tra xem một Point đã tồn tại trong mảng hay chưa
bool pointExists(const Point *uniquePoints, int size, const Point &p) {
    for (int i = 0; i < size; ++i) {
        if (uniquePoints[i].x == p.x && uniquePoints[i].y == p.y) {
            return true;
        }
    }
    return false;
}

// Hàm lọc ra các Point phân biệt từ mảng các Edge
void uniquePoints(const Edge *edges, int size, Point *uniquePoints, int &uniqueSize) {
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

int main() {
    int size = 5;
    Edge edges[5] = {
            {{1, 2},  {3,  4}},
            {{5, 6},  {7,  8}},
            {{1, 2},  {3,  4}},
            {{9, 10}, {11, 12}},
            {{1, 2},  {3,  4}}
    };

    Point in_uniquePoints[size * 2]; // Tối đa có thể là size * 2 vì mỗi Edge có 2 điểm
    int uniqueSize = 0;

    uniquePoints(edges, size, in_uniquePoints, uniqueSize);

    std::cout << "Unique Points:" << std::endl;
    for (int i = 0; i < uniqueSize; ++i) {
        std::cout << "(" << in_uniquePoints[i].x << "," << in_uniquePoints[i].y << ")" << std::endl;
    }

    return 0;
}
