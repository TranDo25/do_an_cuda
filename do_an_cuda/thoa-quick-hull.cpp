#include <iostream>
using namespace std;

struct Point
{
    double x, y;
};

//================================================================================
double findMaxY(Point points[], int n) {
    double maxY = points[0].y;
    for (int i = 1; i < n; i++) {
        if (points[i].y > maxY) {
            maxY = points[i].y;
        }
    }
    return maxY;
}
double findMinY(Point points[], int n) {
    double minY = points[0].y;
    for (int i = 1; i < n; i++) {
        if (points[i].y < minY) {
            minY = points[i].y;
        }
    }
    return minY;
}
double findMaxX(Point points[], int n) {
    double maxX = points[0].y;
    for (int i = 1; i < n; i++) {
        if (points[i].x > maxX) {
            maxX = points[i].x;
        }
    }
    return maxX;
}
double findMinX(Point points[], int n) {
    double minX = points[0].x;
    for (int i = 1; i < n; i++) {
        if (points[i].x < minX) {
            minX = points[i].x;
        }
    }
    return minX;
}
bool comparePoints(const Point& p1, const Point& p2)
{
    return (p1.x < p2.x) || (p1.x == p2.x && p1.y < p2.y);
}
void findPointsByY(Point point[], int size, double maxY, Point foundPoints[], int& foundPointsCount) {

    // Duyệt qua mảng
    for (int i = 0; i < size; i++) {
        // Nếu tọa độ y của phần tử tại vị trí i bằng với maxY
        if (point[i].y == maxY) {
            // Thêm phần tử vào mảng foundPoints
            foundPoints[foundPointsCount++] = point[i];
        }
    }

}
void findPointsByX(Point* points, int n, double x, Point result[], int& n_result) {


    // Duyệt qua mảng points
    for (int i = 0; i < n; i++) {
        // Nếu hoành độ của phần tử thứ i bằng x
        if (points[i].x == x) {
            // Lưu phần tử đó vào mảng result
            result[n_result++] = points[i];
        }
    }
}
void sortPointsByXDescending(Point* points, int n) {
    // Sắp xếp mảng theo hoành độ giảm dần
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            // Nếu hoành độ của phần tử thứ i nhỏ hơn hoành độ của phần tử thứ j
            if (points[i].x < points[j].x) {
                // Đổi chỗ hai phần tử
                Point temp = points[i];
                points[i] = points[j];
                points[j] = temp;
            }
        }
    }
}
void sortPointsByXAscending(Point* points, int n) {
    // Sắp xếp mảng theo hoành độ tăng dần
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            // Nếu hoành độ của phần tử thứ i nhỏ hơn hoành độ của phần tử thứ j
            if (points[i].x > points[j].x) {
                // Đổi chỗ hai phần tử
                Point temp = points[i];
                points[i] = points[j];
                points[j] = temp;
            }
        }
    }
}
void sortPointsByYDescending(Point* points, int n) {
    // Sắp xếp mảng theo tung độ giảm dần
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            // Nếu tung độ của phần tử thứ i nhỏ hơn tung độ của phần tử thứ j
            if (points[i].y < points[j].y) {
                // Đổi chỗ hai phần tử
                Point temp = points[i];
                points[i] = points[j];
                points[j] = temp;
            }
        }
    }
}
void sortPointsByYAscending(Point* points, int n) {
    // Sắp xếp mảng theo tung độ tăng dần
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            // Nếu tung độ của phần tử thứ i lớn hơn tung độ của phần tử thứ j
            if (points[i].y > points[j].y) {
                // Đổi chỗ hai phần tử
                Point temp = points[i];
                points[i] = points[j];
                points[j] = temp;
            }
        }
    }
}
void getPoints1(Point* points, int n, Point q1, Point qq1, Point* result, int& resultSize) {


    // Duyệt qua mảng points
    for (int i = 0; i < n; i++) {
        // Nếu hoành độ của phần tử thứ i nhỏ hơn hoành độ của q1 và tung độ của phần tử thứ i lớn hơn tung độ của qq1
        if (points[i].x <= q1.x && points[i].y >= qq1.y) {
            // Lưu phần tử đó vào mảng result
            result[resultSize++] = points[i];
        }
    }
}
void getPoints2(Point* points, int n, Point qq2, Point q2, Point* result, int& resultSize) {

    // Duyệt qua mảng points
    for (int i = 0; i < n; i++) {
        // Nếu hoành độ của phần tử thứ i nhỏ hơn hoành độ của q1 và tung độ của phần tử thứ i lớn hơn tung độ của qq1
        if (points[i].x <= qq2.x && points[i].y <= q2.y) {
            // Lưu phần tử đó vào mảng result
            result[resultSize++] = points[i];
        }
    }
}
void getPoints3(Point* points, int n, Point q3, Point qq3, Point* result, int& resultSize) {


    // Duyệt qua mảng points
    for (int i = 0; i < n; i++) {
        // Nếu hoành độ của phần tử thứ i nhỏ hơn hoành độ của q1 và tung độ của phần tử thứ i lớn hơn tung độ của qq1
        if (points[i].x >= q3.x && points[i].y <= qq3.y) {
            // Lưu phần tử đó vào mảng result
            result[resultSize++] = points[i];
        }
    }

}
void getPoints4(Point* points, int n, Point qq4, Point q4, Point* result, int& resultSize) {


    // Duyệt qua mảng points
    for (int i = 0; i < n; i++) {
        // Nếu hoành độ của phần tử thứ i nhỏ hơn hoành độ của q1 và tung độ của phần tử thứ i lớn hơn tung độ của qq1
        if (points[i].x >= qq4.x && points[i].y >= q4.y) {
            // Lưu phần tử đó vào mảng result
            result[resultSize++] = points[i];
        }
    }
}
void findPointsWithYGreaterThan(Point points[], int n, double y, Point results[], int& n_result) {
    for (int i = 0; i < n; i++) {
        if (points[i].y > y) {
            results[n_result++] = points[i];
        }
    }
}
void findPointsWithYSmallerThan(Point points[], int n, double y, Point results[], int& n_result) {
    for (int i = 0; i < n; i++) {
        if (points[i].y < y) {
            results[n_result++] = points[i];
        }
    }
}
void findPointsWithXSmallerThan(Point points[], int n, double x, Point results[], int& n_result) {
    for (int i = 0; i < n; i++) {
        if (points[i].x < x) {
            results[n_result++] = points[i];
        }
    }
}
void findPointsWithXGreaterThan(Point points[], int n, double x, Point results[], int& n_result) {
    for (int i = 0; i < n; i++) {
        if (points[i].x > x) {
            results[n_result++] = points[i];
        }
    }
}
//Point* readPoints(const string& filename)
//{
//	ifstream file(filename);
//	Point* points;
//	double x, y;
//	int count = 0;
//	while (file >> x >> y)
//	{
//		points[count++] = { x, y };
//	}
//
//	return points;
//}
void calculateSquareDifference(Point points[], int n, double q1, double qq1, double results[]) {
    double tmp1[1000];
    double tmp2[1000];
    for (int i = 0; i < n; i++) {
        tmp1[i] = (points[i].x - q1) * (points[i].x - q1);
    }
    for (int i = 0; i < n; i++) {
        tmp2[i] = (points[i].y - qq1) * (points[i].y - qq1);
    }
    for (int i = 0; i < n; i++) {
        results[i] = tmp1[i] + tmp2[i];
    }
}
double findMax(double array[], int size) {
    double max = array[0];
    for (int i = 1; i < size; i++) {
        if (array[i] > max) {
            max = array[i];
        }
    }
    return max;
}
void findEqualElements(double key1[], int size, double maxset1, int results[]) {
    int count = 0;
    for (int i = 0; i < size; i++) {
        if (key1[i] == maxset1) {
            results[count] = i;
            count++;
        }
    }
}
void findOHull1(Point set1[], int n_set1, Point q1, Point qq1, Point arrangedPoints[], int& n_arrangedPoints)
{
    if (n_set1 == 0)
    {
        //n_arrangedPoints = 0;
        return;
    }
    double key1[1000];
    int n_key1 = n_set1;
    calculateSquareDifference(set1, n_set1, q1.x, qq1.y, key1);
    double maxset1 = findMax(key1, n_key1);
    int arr_index_newpoint1[100];
    findEqualElements(key1, n_key1, maxset1, arr_index_newpoint1);
    Point new_point1 = set1[arr_index_newpoint1[0]];
    Point new_set11[1000], new_set12[1000];
    int n_new_set11 = 0;
    int n_new_set12 = 0;
    findPointsWithYGreaterThan(set1, n_set1, new_point1.y, new_set11, n_new_set11);
    findPointsWithXSmallerThan(set1, n_set1, new_point1.x, new_set12, n_new_set12);

    findOHull1(new_set11, n_new_set11, q1, new_point1, arrangedPoints, n_arrangedPoints);
    arrangedPoints[n_arrangedPoints++] = new_point1;
    findOHull1(new_set12, n_new_set12, new_point1, qq1, arrangedPoints, n_arrangedPoints);

}
void findOHull2(Point set2[], int n_set2, Point q2, Point qq2, Point arrangedPoints[], int& n_arrangedPoints)
{
    if (n_set2 == 0)
    {
        //n_arrangedPoints = 0;
        return;
    }
    double key2[1000];
    int n_key2 = n_set2;
    calculateSquareDifference(set2, n_set2, qq2.x, q2.y, key2);
    double maxset2 = findMax(key2, n_key2);
    int arr_index_newpoint2[200];

    findEqualElements(key2, n_key2, maxset2, arr_index_newpoint2);
    Point new_point2 = set2[arr_index_newpoint2[0]];
    Point new_set21[1000], new_set22[1000];
    int n_new_set21 = 0;
    int n_new_set22 = 0;
    findPointsWithXSmallerThan(set2, n_set2, new_point2.x, new_set21, n_new_set21);
    findPointsWithYSmallerThan(set2, n_set2, new_point2.y, new_set22, n_new_set22);

    findOHull2(new_set21, n_new_set21, q2, new_point2, arrangedPoints, n_arrangedPoints);
    arrangedPoints[n_arrangedPoints++] = new_point2;
    findOHull2(new_set22, n_new_set22, new_point2, qq2, arrangedPoints, n_arrangedPoints);

}
void findOHull3(Point set3[], int n_set3, Point q3, Point qq3, Point arrangedPoints[], int& n_arrangedPoints)
{
    if (n_set3 == 0)
    {
        //n_arrangedPoints = 0;
        return;
    }
    double key3[1000];
    int n_key3 = n_set3;
    calculateSquareDifference(set3, n_set3, q3.x, qq3.y, key3);
    double maxset3 = findMax(key3, n_key3);
    int arr_index_newpoint3[100];
    findEqualElements(key3, n_key3, maxset3, arr_index_newpoint3);
    Point new_point3 = set3[arr_index_newpoint3[0]];
    Point new_set31[1000], new_set32[1000];
    int n_new_set31 = 0;
    int n_new_set32 = 0;
    findPointsWithYSmallerThan(set3, n_set3, new_point3.y, new_set31, n_new_set31);
    findPointsWithXGreaterThan(set3, n_set3, new_point3.x, new_set32, n_new_set32);

    findOHull3(new_set31, n_new_set31, q3, new_point3, arrangedPoints, n_arrangedPoints);
    arrangedPoints[n_arrangedPoints++] = new_point3;
    findOHull3(new_set32, n_new_set32, new_point3, qq3, arrangedPoints, n_arrangedPoints);
}
void findOHull4(Point set4[], int n_set4, Point q4, Point qq4, Point arrangedPoints[], int& n_arrangedPoints)
{
    if (n_set4 == 0)
    {
        return;
    }
    double key4[1000];
    int n_key4 = n_set4;
    calculateSquareDifference(set4, n_set4, qq4.x, q4.y, key4);
    double maxset4 = findMax(key4, n_key4);
    int arr_index_newpoint4[100];
    findEqualElements(key4, n_key4, maxset4, arr_index_newpoint4);
    Point new_point4 = set4[arr_index_newpoint4[0]];
    Point new_set41[1000], new_set42[1000];
    int n_new_set41 = 0;
    int n_new_set42 = 0;
    findPointsWithXGreaterThan(set4, n_set4, new_point4.x, new_set41, n_new_set41);
    findPointsWithYGreaterThan(set4, n_set4, new_point4.y, new_set42, n_new_set42);

    findOHull4(new_set41, n_new_set41, q4, new_point4, arrangedPoints, n_arrangedPoints);
    arrangedPoints[n_arrangedPoints++] = new_point4;
    findOHull4(new_set42, n_new_set42, new_point4, qq4, arrangedPoints, n_arrangedPoints);

}
void findConvexHull(Point points[], int point_size, Point arranged_points[], int& n_arranged_point)
{

    double maxY = findMaxY(points, point_size);
    double minY = findMinY(points, point_size);
    double maxX = findMaxX(points, point_size);
    double minX = findMinX(points, point_size);

    Point rightPoints[1000];
    int n_rightPoints = 0;
    findPointsByX(points, point_size, maxX, rightPoints, n_rightPoints);

    Point leftPoints[1000];
    int n_leftPoints = 0;
    findPointsByX(points, point_size, minX, leftPoints, n_leftPoints);

    Point topPoints[1000];
    int n_topPoints = 0;
    findPointsByY(points, point_size, maxY, topPoints, n_topPoints);

    Point bottomPoints[1000];
    int n_bottomPoints = 0;
    findPointsByY(points, point_size, minY, bottomPoints, n_bottomPoints);



    //top
    Point top[1000];
    int n_top = 0;
    if (n_topPoints == 1) {
        top[0] = topPoints[0];
        n_top = 1;
    }
    else {
        sortPointsByXAscending(topPoints, n_topPoints);
        n_top = 2;
        top[0] = topPoints[0];
        top[1] = topPoints[n_topPoints - 1];
    }


    //bottom

    Point bottom[1000];
    int n_bottom = 0;
    if (n_bottomPoints == 1) {
        bottom[0] = bottomPoints[0];
        n_bottom = 1;
    }
    else {
        sortPointsByXDescending(bottomPoints, n_bottomPoints);
        n_bottom = 2;
        bottom[0] = bottomPoints[0];
        bottom[1] = bottomPoints[n_bottomPoints - 1];
    }

    //right
    Point right[1000];
    int n_right = 0;
    if (n_rightPoints == 1) {
        right[0] = rightPoints[0];
        n_right = 1;
    }
    else {
        sortPointsByYDescending(rightPoints, n_rightPoints);
        n_right = 2;
        right[0] = rightPoints[0];
        right[1] = rightPoints[n_rightPoints - 1];
    }

    //left
    Point left[1000];
    int n_left = 0;
    if (n_leftPoints == 1) {
        left[0] = leftPoints[0];
        n_left = 1;
    }
    else {
        sortPointsByYAscending(leftPoints, n_leftPoints);
        n_left = 2;
        left[0] = leftPoints[0];
        left[1] = leftPoints[n_leftPoints - 1];
    }
    Point q1, qq1, q2, qq2, q3, qq3, q4, qq4;
    if (n_top == 1) {
        q1 = top[0];
        qq4 = top[0];
    }
    else {
        q1 = top[0];
        qq4 = top[1];
    }
    q4 = right[0];

    if (n_right == 1) {
        qq3 = right[0];
    }
    else {
        qq3 = right[1];
    }
    q3 = bottom[0];

    if (n_bottom == 1) {
        qq2 = bottom[0];
    }
    else {
        qq2 = bottom[1];
    }
    q2 = left[0];

    if (n_left == 1) {
        qq1 = left[0];
    }
    else {
        qq1 = left[1];
    }

    Point set1[1000], set2[1000], set3[1000], set4[1000];
    int n_set1 = 0;
    int n_set2 = 0;
    int n_set3 = 0;
    int n_set4 = 0;
    getPoints1(points, point_size, q1, qq1, set1, n_set1);
    getPoints2(points, point_size, qq2, q2, set2, n_set2);
    getPoints3(points, point_size, q3, qq3, set3, n_set3);
    getPoints4(points, point_size, qq4, q4, set4, n_set4);



    Point new_arranged_points[1000];
    new_arranged_points[0] = q1;
    new_arranged_points[1] = q1;
    int n_new_arranged_points = 2;
    findOHull1(set1, n_set1, q1, qq1, new_arranged_points, n_new_arranged_points);
    new_arranged_points[n_new_arranged_points++] = qq1;
    new_arranged_points[n_new_arranged_points++] = q2;
    findOHull2(set2, n_set2, q2, qq2, new_arranged_points, n_new_arranged_points);
    new_arranged_points[n_new_arranged_points++] = qq2;

    new_arranged_points[n_new_arranged_points++] = q3;
    findOHull3(set3, n_set3, q3, qq3, new_arranged_points, n_new_arranged_points);
    new_arranged_points[n_new_arranged_points++] = qq3;

    new_arranged_points[n_new_arranged_points++] = q4;
    findOHull4(set4, n_set4, q4, qq4, new_arranged_points, n_new_arranged_points);
    new_arranged_points[n_new_arranged_points++] = qq4;
    for (int i = 0; i < n_new_arranged_points; i++) {
        arranged_points[i] = new_arranged_points[i];
    }
    n_arranged_point = n_new_arranged_points;

}


int main()
{
    // tao diem
    Point points[10000];
    points[0] = { 1.0, 1.0 };
    points[1] = { 2.0, 2.0 };
    points[2] = { 3.0, 1.5 };
    points[3] = { 4.0, 5.0 };
    points[4] = { 5.0, 4.0 };
    points[5] = { 6.0, 3.0 };
    points[6] = { 7.0, 2.0 };
    points[7] = { 8.0, 3.5 };
    points[8] = { 9.0, 1.0 };
    int n_points = 9;
    Point convexHull[1000];
    int n_convexHull = 0;
    findConvexHull(points, n_points, convexHull, n_convexHull);

    cout << "Convex Hull Points: "<<n_convexHull<<" diem" << endl;
    for (int i = 0; i<n_convexHull; i++)
    {
        cout << "(" << convexHull[i].x << ", " << convexHull[i].y << ")" << endl;
    }
    return 0;
}
