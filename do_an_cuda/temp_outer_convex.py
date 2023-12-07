#import các thư viện cần thiết
import numpy as np
import math
import matplotlib.pyplot as plt
np.random.seed(42)

#Khởi tạo 
N_square = 10000

k=16
count = 0
count1 = 0
count2 = 0
count3 = 0
count4 = 0

def outer_convex_approximation(X, δ):
    #Rotation angle
    alpha = -math.pi/2

    #Rotation Matrix
    R = np.array([[math.cos(alpha), math.sin(alpha)],[-math.sin(alpha), math.cos(alpha)]])
#     print("Tập X: ", X)
    # Step I: Determine D, βd, and P
    D = [(1, 0), (0, 1), (-1, 0), (0, -1)] #khởi tạo D theo công thức (7)
#     βd = [np.max(np.dot(X, d)) for d in D] #khởi tạo βd theo công thức (6)
    #công thức (8)
    min_x = np.min((X[:, 0]))
    max_x = np.max(X[:, 0])
    min_y = np.min((X[:, 1]))
    max_y = np.max(X[:, 1])
    
    
    #các giá trị r1, r2, r3, r4 đặt theo công thức (10)
    r1 = (max_x, max_y)
    r2 = (min_x, max_y)
    r3 = (min_x, min_y)
    r4 = (max_x, min_y)
    
#     print(r1, r2, r3, r4)
    
    P = np.array([r1, r2, r3, r4]) #Khởi tạo P theo công thức (9)

    #STEP II: gán Pdoubt = P
    Pdoubt = P[:]
    Ptest = P[:]
    global count, count4, count1, count2, count3
    # Step III: Iterate until Pdoubt is empty
    #Đặt biến count để xem số lần chạy

    while Pdoubt.shape[0] > 0: #Nếu số đỉnh Pdoubt > 0, tiếp tục lặp
        count += 1
        pdoubt = Pdoubt[0] #Lấy một giá trị pdoubt ∈ Pdoubt
        
        #Xác định số liền trước và liền sau ngược chiều kim đồng hồ của pdoubt
        pdoubt_index_idx = np.where(np.all(Ptest == pdoubt, axis=1))[0][0]
        pdoubt_minus = Ptest[(pdoubt_index_idx - 1 ) % len(Ptest)] #liền trước
        pdoubt_plus = Ptest[(pdoubt_index_idx + 1) % len(Ptest)] #liền sau
        
        
        #Xác định dp
        dp = np.dot(R, ((pdoubt_minus - pdoubt_plus).T)) / np.linalg.norm(pdoubt_minus - pdoubt_plus)
        
        
        #xác định βdp
        βdp = np.max(np.dot(dp, X.T))
    
        
        #kiểm tra công thức (16)
        
        if βdp == np.dot(dp, pdoubt_plus):
            count1 += 1 
            #Loại bỏ các phần tử khỏi mảng theo (17)
            #print("trường hợp 1")
            D.append(dp)  # Thêm dp vào danh sách D

            
            P = np.delete(P, np.where(np.all(P == pdoubt, axis=1))[0], axis=0)
            Pdoubt = np.delete(Pdoubt, np.where(np.all(Pdoubt == pdoubt, axis=1))[0], axis=0)
            Ptest = np.delete(Ptest, np.where(np.all(Ptest == pdoubt, axis=1))[0], axis=0)

            
        #Kiểm tra (18) và (19)
        elif np.dot(dp, pdoubt.T) - βdp > δ:
            count2 += 1
            #print("Trường hợp 2")
            #Khởi tạo (20)
            λp = (βdp - np.dot(dp, pdoubt_minus.T)) / (np.dot(dp, pdoubt.T) - np.dot(dp, pdoubt_minus.T)) #λp
            p_hat_minus = (1 - λp) * (pdoubt_minus.T) + λp * (pdoubt.T) #p^
            p_hat_plus = (1 - λp) * (pdoubt_plus.T) + λp * (pdoubt.T) #p^+
            
            #print("p_hat_minus", p_hat_minus)
            #print("p_hat_plus", p_hat_plus)
            #Thêm xóa theo (21)
            D.append((dp[0], dp[1]))  # Thêm dp vào danh sách D
            # Tìm chỉ số của pdoubt trong mảng Pdoubt
            pdoubt_indexp = np.where(np.all(P == pdoubt, axis=1))[0][0]
            
            #thêm (23)
            # Tìm chỉ số của pdoubt trong mảng Pdoubt
            pdoubt_index = np.where(np.all(Pdoubt == pdoubt, axis=1))[0][0]

            # Xóa pdoubt khỏi mảng Pdoubt và thêm p_hat_minus và p_hat_plus
            Pdoubt = np.delete(Pdoubt, pdoubt_index, axis=0)
                    
            if np.allclose(p_hat_plus, pdoubt) and np.allclose(p_hat_minus, pdoubt):
                count4 += 1
                ptest_index = np.where(np.all(Ptest == pdoubt, axis=1))[0][0]
                Ptest = np.concatenate((Ptest[:ptest_index], Ptest[ptest_index+1:], [Ptest[ptest_index]]))

            else:
                # Xóa pdoubt khỏi mảng P
                P = np.delete(P, pdoubt_indexp, axis=0)
                
                #test với Ptest
                ptest_index = np.where(np.all(Ptest == pdoubt, axis=1))[0][0]
                Ptest = np.delete(Ptest, ptest_index, axis=0)
                
                if np.allclose(p_hat_plus, pdoubt_plus):
                    print("1")
                else:
                    P = np.insert(P, pdoubt_indexp, p_hat_plus, axis=0)
                    Pdoubt = np.insert(Pdoubt, pdoubt_index, p_hat_plus, axis=0)
                    Ptest = np.insert(Ptest, ptest_index, p_hat_plus, axis=0)
                if np.allclose(p_hat_minus, pdoubt_minus):
                    print("2")
                else:
                    P = np.insert(P, pdoubt_indexp, p_hat_minus, axis=0)
                    Pdoubt = np.insert(Pdoubt, pdoubt_index, p_hat_minus, axis=0)
                    Ptest = np.insert(Ptest, ptest_index, p_hat_minus, axis=0)


        #các trường hợp còn lại
        else:
            count3 += 1
            #print("trường hợp 3")
            Pdoubt = np.delete(Pdoubt, np.where(np.all(Pdoubt == pdoubt, axis=1))[0][0], axis=0) #xóa pboubt khỏi Pboubt
            
            ptest_index = np.where(np.all(Ptest == pdoubt, axis=1))[0]
            Ptest = np.concatenate((Ptest[:ptest_index[0]], Ptest[ptest_index[0]+1:], [Ptest[ptest_index[0]]]))
        
    # Step V: Return D, βd, and P
    return D, P


# Test the algorithm with data
p1 = np.array([-130.658, -128])
p2 = np.array([-87.522, -128])
p3 = np.array([48.95, -104.547])
p4 = np.array([97.871, -89.647])
p5 = np.array([124.452, -69.508])
p6 = np.array([143.815, -35.205])
p7 = np.array([168, 62.629])
p8 = np.array([168, 96])
p9 = np.array([127.717, 128])
p10 = np.array([82.477, 128])
p11 = np.array([-61.85, 109.243])
p12 = np.array([-102.638, 98.195])
p13 = np.array([-124.269, 81.396])
p14 = np.array([-140.644, 46.529])
p15 = np.array([-168, -67.84])
p16 = np.array([-168, -99.849])

#List of set of cut points
P = np.array([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10,p11, p12, p13,p14, p15, p16, p1])


#Khởi tạo X

#Rotation angle
alpha = -math.pi/2

#Rotation Matrix
R = np.array([[math.cos(alpha), math.sin(alpha)],[-math.sin(alpha), math.cos(alpha)]])

#Rotate the set of cut points
po = np.dot(P, R)


# Determine min, max of p
#công thức (8)
x_min = np.nanmin(po[:, 0])
x_max = np.nanmax(po[:, 0])
y_min = np.nanmin(po[:, 1])
y_max = np.nanmax(po[:, 1])

xy_min = [x_min, y_min]
xy_max = [x_max , y_max]


#Khởi tạo X
X = np.random.uniform(low=xy_min, high=xy_max, size=(N_square, 2))
def outside( points, a, b ):
    res = np.array([])
    nx = b[1] - a[1]
    ny = a[0] - b[0]
    d =  nx * a[0] + ny * a[1]
    l =  nx * points[:,0] + ny* points[:,1]
    res = np.append(res, l-d >= 0)
    return res


#Delete the points outside the polygon to creat a test set
for j in range( k ):
    X = np.delete(X, np.where(outside(X[:], po[j], po[j+1])), 0)
X = X[0:50]


δ = 0.0
D, P = outer_convex_approximation(X, δ)

print("D:", D)
# print("βd:", βd)
print("P: ",P)
print("Số đỉnh: ", P.shape[0])
print("Số lượt chạy qua step III: ", count)

#print("Số trường hợp 1: ",count1)
#print("Số trường hợp 2: ",count2)
#print("Số trường hợp 3: ",count3)
#print("count4: ",count4)
# Plotting
from scipy.spatial import ConvexHull

# Tạo các mảng x và y từ tập P
x = P[:, 0]
y = P[:, 1]

# Tìm bao lồi của tập điểm
points = np.column_stack((x, y))
hull = ConvexHull(points)

# Lấy các điểm trên bao lồi
convex_points = points[hull.vertices]

# Thêm điểm đầu tiên vào cuối cùng để tạo thành đa giác hoàn chỉnh
convex_points = np.vstack((convex_points, convex_points[0]))

# Plotting
fig, ax = plt.subplots(figsize=(8, 8))
ax.plot(X[:, 0], X[:, 1], 'r.')
ax.plot(P[:, 0], P[:, 1], 'r.')
ax.plot(convex_points[:, 0], convex_points[:, 1], 'k-', linewidth=2)
ax.set_xlabel('X', fontsize=12)
ax.set_ylabel('Y', fontsize=12)
ax.set_title('Outer Convex Approximation', fontsize=14)
plt.grid(True)
plt.show()
