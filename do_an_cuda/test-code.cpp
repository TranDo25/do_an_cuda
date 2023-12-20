//#include <iostream>
//
//struct Point {
//    int x;
//    int y;
//};
//
//// Hàm đảo ngược các phần tử trong mảng Point
//void reversePoints(Point *points, int size) {
//    int start = 0;
//    int end = size - 1;
//
//    while (start < end) {
//        // Swap giữa points[start] và points[end]
//        Point temp = points[start];
//        points[start] = points[end];
//        points[end] = temp;
//
//        // Di chuyển lần lượt đến phần tử thứ hai và thứ hai từ cuối
//        ++start;
//        --end;
//    }
//}
//
//int main() {
//    const int size = 5;
//    Point points[size] = {
//            {1, 2},
//            {3, 4},
//            {5, 6},
//            {7, 8},
//            {9, 10}
//    };
//
//    std::cout << "Original Points:" << std::endl;
//    for (int i = 0; i < size; ++i) {
//        std::cout << "(" << points[i].x << "," << points[i].y << ")" << std::endl;
//    }
//
//    reversePoints(points, size);
//
//    std::cout << "\nReversed Points:" << std::endl;
//    for (int i = 0; i < size; ++i) {
//        std::cout << "(" << points[i].x << "," << points[i].y << ")" << std::endl;
//    }
//
//    return 0;
//}
