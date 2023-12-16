//#include <iostream>
//
//struct Point {
//    double x, y;
//};
//
//
//int findMinPointIndex(const Point points[], int n) {
//    if (n <= 0) {
//        return -1;
//    }
//
//    int minIndex = 0;
//
//    for (int i = 1; i < n; ++i) {
//        if (points[i].y < points[minIndex].y) {
//            minIndex = i;
//        }
//        else if (points[i].y == points[minIndex].y && points[i].x < points[minIndex].x) {
//            minIndex = i;
//        }
//    }
//
//    return minIndex;
//}
//
//void moveMinPointsToFront(Point points[], int n) {
//    int minIndex = findMinPointIndex(points, n);
//    Point points_premitive[9];
//    for (int i = 0; i < n; i++) {
//        points_premitive[i] = points[i];
//    }
//    if (minIndex != -1) {
//        Point temp[20];
//        int tempIndex = 0;
//
//        // Lưu trữ các phần tử từ minIndex đến cuối mảng vào mảng tạm thời.
//        for (int i = minIndex; i < n; ++i) {
//            temp[tempIndex++] = points[i];
//        }
//
//        // Di chuyển các phần tử từ 0 đến minIndex - 1 lên một vị trí.
//        for (int i = 0; i < minIndex; ++i) {
//            points[i + n - minIndex] = points_premitive[i];
//        }
//
//        // Đặt các phần tử có hoành độ và tung độ nhỏ nhất lên đầu mảng.
//        for (int i = 0; i < tempIndex; ++i) {
//            points[i] = temp[i];
//        }
//    }
//}
//
//int main() {
//    int n = 5;
//    Point points[20] = { {-11, -2}, {3, 4}, {0, 1}, {5, 6}, {-7, 8} };
//
//    std::cout << "Mảng điểm trước khi di chuyển: " << std::endl;
//    for (const auto& point : points) {
//        std::cout << "(" << point.x << ", " << point.y << ") ";
//    }
//    std::cout << std::endl;
//
//    // Di chuyển toàn bộ các phần tử có hoành độ và tung độ nhỏ nhất về sau lên đầu mảng.
//    moveMinPointsToFront(points, n);
//
//    std::cout << "Mảng điểm sau khi di chuyển: " << std::endl;
//    for (const auto& point : points) {
//        std::cout << "(" << point.x << ", " << point.y << ") ";
//    }
//    std::cout << std::endl;
//
//    return 0;
//}
