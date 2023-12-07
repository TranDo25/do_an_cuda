#import các thư viện cần thiết
import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.spatial import distance
from scipy.spatial import ConvexHull
np.random.seed(42)

#Khởi tạo 
N_square = 10000

k=16

count = 0
count1 = 0
count2 = 0
count3 = 0
count4 = 0
count5 = 0
count6 = 0


def inner_convex_approximation(X, δ):
#     print("X: ",X)
    #Rotation angle
    alpha = -math.pi/2

    #Rotation Matrix
    R = np.array([[math.cos(alpha), math.sin(alpha)],[-math.sin(alpha), math.cos(alpha)]])
    #công thức (44)
    min_x = np.min(X[:, 0])
    max_x = np.max(X[:, 0])
    min_y = np.min(X[:, 1])
    max_y = np.max(X[:, 1])

    
    min_X = X[X[:, 0] == min_x]
    max_X = X[X[:, 0] == max_x]
    min_Y = X[X[:, 1] == min_y]
    max_Y = X[X[:, 1] == max_y]
    
    q12 = np.max(max_X[:, 1])
    q21 = np.min(max_Y[:, 0])
    q32 = np.min(min_X[:, 1])
    q41 = np.max(min_Y[:, 0])
    
    q1 = (max_x, q12)
    q2 = (q21, max_y)
    q3 = (min_x, q32)
    q4 = (q41, min_y)
    
    Xcomma = np.array([q1, q2, q3, q4])
    E_test = np.array([[q1, q2], [q2, q3], [q3, q4], [q4, q1]])
        
    Edoubt = E_test[:]
    
    global count, count1, count2, count3, count4, count5, count6
    while Edoubt.shape[0] > 0:
        count += 1
        #print("Lượt chạy thứ: ", count)
        edoubt = Edoubt[0]
        p, p_plus = edoubt[0], edoubt[1]
        
        dpp_plus = np.dot(R, ((p_plus - p).T)) / np.linalg.norm(p_plus - p)
        
        Xpp_plus = []
        for x in X:
            if np.dot(dpp_plus, x.T) > np.dot(dpp_plus, p.T):
                Xpp_plus.append(x)
        
        Xpp_plus_fix = np.array(Xpp_plus)
        

        if Xpp_plus_fix.size != 0:
            
            βpp_plus = np.max(np.dot(dpp_plus, Xpp_plus_fix.T))
            
            
            Bpp_plus = []
            for x in Xpp_plus_fix:
                if np.dot(dpp_plus, x.T) == βpp_plus:
                    Bpp_plus = np.append(Bpp_plus, x)
            
            
            
            if βpp_plus - np.dot(dpp_plus, p.T) <= δ:
                count1 += 1
                Edoubt = np.delete(Edoubt, np.where(np.all(Edoubt == edoubt, axis=1))[0], axis=0)
            
            else:
                p_hat = []
                p_hat_calc = Bpp_plus[:]
                x_calc = Bpp_plus[:]
               
                if np.linalg.norm(p_hat_calc - p) == np.max(np.linalg.norm(x_calc - p)):
                        p_hat = np.append(p_hat, p_hat_calc)
                
                
                if np.allclose(p_hat, p) and not np.allclose(p_hat, p_plus):
                    count2 += 1
                    print("1")
                    edoubt_index = np.where(np.all(Edoubt == edoubt, axis=(1,2)))[0][0]
                    Edoubt = np.delete(Edoubt, edoubt_index, axis=0)
                elif np.allclose(p_hat, p_plus) and not np.allclose(p_hat, p):
                    count3 += 1
                    print("2")
                    edoubt_index = np.where(np.all(Edoubt == edoubt, axis=(1,2)))[0][0]
                    Edoubt = np.delete(Edoubt, edoubt_index, axis=0)
                elif np.allclose(p_hat, p_plus) and np.allclose(p_hat, p):
                    count6 +=1
                else:
                    count4 += 1
                    dpp_hat = np.dot(R, ((p_hat - p).T)) / np.linalg.norm(p_hat - p)

                    dp_hat_p_plus = np.dot(R, ((p_plus - p_hat).T)) / np.linalg.norm(p_plus - p_hat)

                    dp_hat_p_T = np.dot(R, (((p.T) - p_hat).T)) / np.linalg.norm((p.T) - p_hat)

                    Xpp_hat_find = []
                    for x in Xpp_plus:
                        if np.dot(dpp_hat, x.T) > np.dot(dpp_hat, p.T):
                            Xpp_hat_find.append(x)
                    Xpp_hat = np.array(Xpp_hat_find)

                    Xp_hat_p_plus_find = []
                    for x in Xpp_plus:
                        if np.dot(dp_hat_p_plus, x.T) > np.dot(dp_hat_p_T, p_hat.T):
                            Xp_hat_p_plus_find.append(x)
                    Xp_hat_p_plus = np.array(Xp_hat_p_plus_find)


                    Xcomma = np.vstack((Xcomma, p_hat))
                    
                    
                    A = np.array([p, p_hat])
                    B = np.array([p_hat, p_plus])

                    #Xác định vị trí và xóa cạnh [p, p+] và thêm 2 cạnh [p+, p^] và [p^, p] vào E
                    e_index = np.where(np.all(E_test == edoubt, axis=(1,2)))[0][0]
                    E_test = np.delete(E_test, e_index, axis=0)
                    E_test = np.insert(E_test, e_index, [A, B], axis=0)

                    #Xác định vị trí và xóa cạnh [p, p+] và thêm 2 cạnh [p+, p^] và [p^, p]
                    edoubt_index = np.where(np.all(Edoubt == edoubt, axis=(1,2)))[0][0]
                    Edoubt = np.delete(Edoubt, edoubt_index, axis=0)
                    Edoubt = np.insert(Edoubt, edoubt_index, [A, B], axis=0)

        
        else:
            count5 += 1
            Edoubt = np.delete(Edoubt, np.where(np.all(Edoubt == edoubt, axis=1))[0], axis=0)
        
    # Step IV: Return X', E
    Xcomma = np.unique(Xcomma, axis=0)
    hull = ConvexHull(Xcomma)

    # Lấy các cạnh của đa giác lồi
    E = []
    for simplex in hull.simplices:
        E.append(Xcomma[simplex])

    E = np.array(E)
    E_test = np.unique(E_test, axis=0)
    return Xcomma, E
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
X = X[0:9000]


δ = 0.0

Xcomma, E = inner_convex_approximation(X, δ)

print("Xcomma", Xcomma)
print("E", E)
print("Số đỉnh tập Xcomma: ", Xcomma.shape[0])
print("Số cạnh tập E: ", E.shape[0])


print("Số lần lặp qua step III: ", count)

# X_uni = np.unique(Xcomma, axis=0)
# print("Số đỉnh: ", X_uni.shape[0])
# Plotting
# from scipy.spatial import ConvexHull

# Tạo các mảng x và y từ tập P
x = Xcomma[:, 0]
y = Xcomma[:, 1]

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
ax.plot(Xcomma[:, 0], Xcomma[:, 1], 'r.')
ax.plot(convex_points[:, 0], convex_points[:, 1], 'k-', linewidth=2)
ax.set_xlabel('X', fontsize=12)
ax.set_ylabel('Y', fontsize=12)
ax.set_title('Inner Convex Approximation', fontsize=14)
plt.grid(True)
plt.show()
