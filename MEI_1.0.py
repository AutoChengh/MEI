# 这版融合了（1）碰撞检测；（2）只有两车最近距离点仍在靠近时，才计算SSM数值
# 包括的SSMs有：(1)ACT; (2)RTTC(2DTTC); (3)MEI
import numpy as np
import math


# 检查碰撞的辅助函数1
def get_projection_offsets(length_1, width_1, heading_1, length_2, width_2, heading_2):
    dx_1 = length_1 / 2
    dy_1 = width_1 / 2
    dx_2 = length_2 / 2
    dy_2 = width_2 / 2
    vertices_1 = np.array([
        [dx_1, dy_1],
        [-dx_1, dy_1],
        [-dx_1, -dy_1],
        [dx_1, -dy_1]
    ])
    vertices_2 = np.array([
        [dx_2, dy_2],
        [-dx_2, dy_2],
        [-dx_2, -dy_2],
        [dx_2, -dy_2]
    ])
    rotation_matrix_1 = np.array([
        [np.cos(heading_1), -np.sin(heading_1)],
        [np.sin(heading_1), np.cos(heading_1)]
    ])
    rotation_matrix_2 = np.array([
        [np.cos(heading_2), -np.sin(heading_2)],
        [np.sin(heading_2), np.cos(heading_2)]
    ])
    rotated_vertices_1 = np.dot(vertices_1, rotation_matrix_1.T)
    rotated_vertices_2 = np.dot(vertices_2, rotation_matrix_2.T)
    axes = []
    for i in range(2):
        edge = rotated_vertices_1[i] - rotated_vertices_1[i - 1]
        axis = np.array([-edge[1], edge[0]])
        axes.append(axis / np.linalg.norm(axis))
    for i in range(2):
        edge = rotated_vertices_2[i] - rotated_vertices_2[i - 1]
        axis = np.array([-edge[1], edge[0]])
        axes.append(axis / np.linalg.norm(axis))
    projections_1 = [np.dot(rotated_vertices_1, axis) for axis in axes]
    max_and_min_projections_1 = []
    for i in projections_1:
        max_and_min_projections_1.append([np.min(i), np.max(i)])
    projections_2 = [np.dot(rotated_vertices_2, axis) for axis in axes]
    max_and_min_projections_2 = []
    for i in projections_2:
        max_and_min_projections_2.append([np.min(i), np.max(i)])
    return axes, max_and_min_projections_1, max_and_min_projections_2

# 检查碰撞的辅助函数2
def check_collisions_between_series(A_series, B_series, axes, max_and_min_projections_1, max_and_min_projections_2):
    # A_series 包含车辆的中心点坐标
    # 用numpy并行计算A_series的投影
    proj_A_min_max = []
    proj_B_min_max = []
    for i in range(4):
        proj_A = np.dot(A_series[:, :2], axes[i])
        proj_A_min_max.append([proj_A + max_and_min_projections_1[i][0], proj_A + max_and_min_projections_1[i][1]])
        proj_B = np.dot(B_series[:, :2], axes[i])
        proj_B_min_max.append([proj_B + max_and_min_projections_2[i][0], proj_B + max_and_min_projections_2[i][1]])
    proj_A_min_max = np.array(proj_A_min_max)
    proj_B_min_max = np.array(proj_B_min_max)
    # 检查投影是否重叠
    # 对A B中的每两辆车，在4条轴上检查是否有重叠
    # 用numpy并行计算
    if_collision = []
    for i in range(4):
        if_collision.append(np.logical_not(
            (proj_A_min_max[i][1][:, np.newaxis] < proj_B_min_max[i][0][np.newaxis, :]) | (
                        proj_B_min_max[i][1][np.newaxis, :] < proj_A_min_max[i][0][:, np.newaxis])))
    if_collision = np.array(if_collision)
    if_collision = np.all(if_collision, axis=0)
    return if_collision

# 检查碰撞的辅助函数3
def get_collision_A_and_B(x_A, y_A, x_B, y_B, h_A, h_B, l_A, w_A, l_B, w_B):
    x_A_array, y_A_array = np.array([x_A]), np.array([y_A])
    x_B_array, y_B_array = np.array([x_B]), np.array([y_B])
    axes, max_and_min_projections_1, max_and_min_projections_2 = get_projection_offsets(l_A, w_A, h_A,
                                                                                        l_B, w_B, h_B)
    A_series = np.array([x_A_array, y_A_array]).T
    B_series = np.array([x_B_array, y_B_array]).T
    if_collision = check_collisions_between_series(A_series, B_series, axes, max_and_min_projections_1,
                                                   max_and_min_projections_2)
    return if_collision



# 计算RTTC的辅助函数1
def is_ray_intersect_segment(ray_origin_x, ray_origin_y, ray_direction_x, ray_direction_y,
                             segment_start_x, segment_start_y, segment_end_x, segment_end_y):
    ray_origin = np.array([ray_origin_x, ray_origin_y])
    ray_direction = np.array([ray_direction_x, ray_direction_y])
    segment_start = np.array([segment_start_x, segment_start_y])
    segment_end   = np.array([segment_end_x,   segment_end_y])

    v1 = ray_origin - segment_start
    v2 = segment_end - segment_start
    v3 = np.array([-ray_direction[1], ray_direction[0]])
    v3 = v3 / np.linalg.norm(v3)

    dot = np.dot(v2, v3)
    if abs(dot) < 1e-10:
        if abs(np.cross(v1, v2)) < 1e-10:
            t0 = np.dot(segment_start - ray_origin, ray_direction)
            t1 = np.dot(segment_end - ray_origin, ray_direction)
            if t0 >= 0 and t1 >= 0:
                return min(t0, t1)
            if t0 < 0 and t1 < 0:
                return None
            return 0
        return None

    t1 = np.cross(v2, v1) / dot
    t2 = np.dot(v1, v3) / dot

    if 0 <= t2 <= 1:
        return t1
    return None

# 计算RTTC的辅助函数2
def compute_RTTC(x_A, y_A, v_A, h_A, l_A, w_A, x_B, y_B, v_B, h_B, l_B, w_B):
    rotate_matrix_A = np.array([[np.cos(h_A), np.sin(h_A)], [-np.sin(h_A), np.cos(h_A)]])
    rotate_matrix_B = np.array([[np.cos(h_B), np.sin(h_B)], [-np.sin(h_B), np.cos(h_B)]])
    #20250221修改
    bbox_A = np.array([x_A, y_A]) + np.dot(
        np.array([[l_A / 2, w_A / 2], [l_A / 2, -w_A / 2], [-l_A / 2, -w_A / 2], [-l_A / 2, w_A / 2]]), rotate_matrix_A)
    bbox_B = np.array([x_B, y_B]) + np.dot(
        np.array([[l_B / 2, w_B / 2], [l_B / 2, -w_B / 2], [-l_B / 2, -w_B / 2], [-l_B / 2, w_B / 2]]), rotate_matrix_B)
    v_A_array = np.array([v_A * np.cos(h_A), v_A * np.sin(h_A)])
    v_B_array = np.array([v_B * np.cos(h_B), v_B * np.sin(h_B)])
    v_rel = v_A_array - v_B_array
    DTC = np.nan
    for i in range(4):
        dist_has_nagative = False
        for j in range(4):
            dist = is_ray_intersect_segment(
                bbox_A[i][0], bbox_A[i][1],
                v_rel[0], v_rel[1],
                bbox_B[j][0], bbox_B[j][1],
                bbox_B[(j + 1) % 4][0], bbox_B[(j + 1) % 4][1]
            )
            if dist is not None:
                if np.isnan(DTC):
                    DTC = dist
                if dist > 0:
                    DTC = min(DTC, dist)
                if dist < 0:
                    dist_has_nagative = True
                if dist_has_nagative and dist > 0:
                    return 0
    for i in range(4):
        dist_has_nagative = False
        for j in range(4):
            dist = is_ray_intersect_segment(
                bbox_B[i][0], bbox_B[i][1],
                -v_rel[0], -v_rel[1],
                bbox_A[j][0], bbox_A[j][1],
                bbox_A[(j + 1) % 4][0], bbox_A[(j + 1) % 4][1]
            )
            if dist is not None:
                if np.isnan(DTC):
                    DTC = dist
                if dist >= 0:
                    DTC = min(DTC, dist)
                if dist < 0:
                    if DTC < 0:
                        DTC = max(DTC, dist)
                    dist_has_nagative = True
                if dist_has_nagative and dist > 0:
                    return 0
    if not np.isnan(DTC) and np.linalg.norm(v_rel) > 1e-12:
        RTTC = DTC / np.linalg.norm(v_rel)
        return RTTC, DTC, np.linalg.norm(v_rel)
    return np.nan, np.nan, np.nan



# 计算ACT的辅助函数1
def get_rect_corners(x, y, h, l, w):
    # 矩形四个角的坐标相对中心点
    corners = [
        [-l / 2 * math.cos(h) - w / 2 * math.sin(h), -l / 2 * math.sin(h) + w / 2 * math.cos(h)],  # A1
        [l / 2 * math.cos(h) - w / 2 * math.sin(h), l / 2 * math.sin(h) + w / 2 * math.cos(h)],  # A2
        [l / 2 * math.cos(h) + w / 2 * math.sin(h), l / 2 * math.sin(h) - w / 2 * math.cos(h)],  # A3
        [-l / 2 * math.cos(h) + w / 2 * math.sin(h), -l / 2 * math.sin(h) - w / 2 * math.cos(h)],  # A4
    ]
    # 计算绝对坐标
    corners = [[x + corner[0], y + corner[1]] for corner in corners]
    return corners

# 计算ACT的辅助函数2 计算两点之间的距离
def distance(p1, p2):
    return math.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)

# 计算ACT的辅助函数3 计算点到线段的最短距离
def point_to_segment_distance(p, v1, v2):
    # 向量计算
    line_len = distance(v1, v2)
    if line_len == 0:  # 避免除以零
        return distance(p, v1)
    # 计算点到线段的投影
    t = max(0, min(1, ((p[0] - v1[0]) * (v2[0] - v1[0]) + (p[1] - v1[1]) * (v2[1] - v1[1])) / line_len ** 2))
    # 计算最近点
    closest = [v1[0] + t * (v2[0] - v1[0]), v1[1] + t * (v2[1] - v1[1])]
    return distance(p, closest), closest

# 计算ACT的辅助函数4  计算矩形之间的最短距离
def get_shortest_distance(corners_A, corners_B):
    min_distance = float('inf')
    closest_points = None
    closest_A = None
    closest_B = None
    # 计算每个角点到另一个矩形四条边的最短距离
    for i in range(4):
        for j in range(4):
            # A的角点到B的边
            p1 = corners_A[i]
            for k in range(4):
                p2 = corners_B[k]
                p3 = corners_B[(k + 1) % 4]
                dist, closest = point_to_segment_distance(p1, p2, p3)
                if dist < min_distance:
                    min_distance = dist
                    closest_A = p1
                    closest_B = closest
            # B的角点到A的边
            p1 = corners_B[i]
            for k in range(4):
                p2 = corners_A[k]
                p3 = corners_A[(k + 1) % 4]
                dist, closest = point_to_segment_distance(p1, p2, p3)
                if dist < min_distance:
                    min_distance = dist
                    closest_A = closest
                    closest_B = p1
    return min_distance, closest_A, closest_B

# 计算ACT的辅助函数5
def compute_shortest_distance(x_A, y_A, v_A, h_A, l_A, w_A, x_B, y_B, v_B, h_B, l_B, w_B):
    # 计算两个矩形的四个角点
    corners_A = get_rect_corners(x_A, y_A, h_A, l_A, w_A)
    corners_B = get_rect_corners(x_B, y_B, h_B, l_B, w_B)
    # 获取碰撞检测结果
    # 计算最短距离
    min_distance, closest_A, closest_B = get_shortest_distance(corners_A, corners_B)
    delta_x = closest_B[0] - closest_A[0]
    delta_y = closest_B[1] - closest_A[1]
    norm_delta = np.sqrt(delta_x ** 2 + delta_y ** 2)
    if norm_delta != 0:
        unit_vector = np.array([delta_x / norm_delta, delta_y / norm_delta])
        velocity_diff = np.array([v_B * np.cos(h_B) - v_A * np.cos(h_A), v_B * np.sin(h_B) - v_A * np.sin(h_A)])
        v_closest = -np.dot(unit_vector, velocity_diff)
    else:
        v_closest = 0
    return min_distance, (closest_A, closest_B), v_closest



# 计算EI的辅助函数1
def compute_v_Br(x_A, y_A, v_A, h_A, x_B, y_B, v_B, h_B):
    delta_x = x_B - x_A
    delta_y = y_B - y_A
    norm_delta = np.sqrt(delta_x ** 2 + delta_y ** 2)
    if norm_delta != 0:
        unit_vector = np.array([delta_x / norm_delta, delta_y / norm_delta])
        velocity_diff = np.array([v_B * np.cos(h_B) - v_A * np.cos(h_A), v_B * np.sin(h_B) - v_A * np.sin(h_A)])
        v_Br = -np.dot(unit_vector, velocity_diff)
    else:
        v_Br = 0
    return v_Br

# 计算EI的辅助函数2
def compute_InDepth(x_A, y_A, v_A, h_A, l_A, w_A, x_B, y_B, v_B, h_B, l_B, w_B):
    v_diff = np.array([v_B * np.cos(h_B) - v_A * np.cos(h_A), v_B * np.sin(h_B) - v_A * np.sin(h_A)])
    theta_B_prime = v_diff / np.linalg.norm(v_diff)
    delta = np.array([x_B - x_A, y_B - y_A])
    D_t1 = np.linalg.norm(delta - np.dot(delta, theta_B_prime) * theta_B_prime)

    AA1 = np.array([l_A / 2 * np.cos(h_A) - w_A / 2 * -np.sin(h_A), l_A / 2 * np.sin(h_A) - w_A / 2 * np.cos(h_A)])
    AA2 = np.array([l_A / 2 * np.cos(h_A) + w_A / 2 * -np.sin(h_A), l_A / 2 * np.sin(h_A) + w_A / 2 * np.cos(h_A)])
    AA3 = np.array([-l_A / 2 * np.cos(h_A) - w_A / 2 * -np.sin(h_A), -l_A / 2 * np.sin(h_A) - w_A / 2 * np.cos(h_A)])
    AA4 = np.array([-l_A / 2 * np.cos(h_A) + w_A / 2 * -np.sin(h_A), -l_A / 2 * np.sin(h_A) + w_A / 2 * np.cos(h_A)])
    d_A1 = np.linalg.norm(AA1 - np.dot(AA1, theta_B_prime) * theta_B_prime)
    d_A2 = np.linalg.norm(AA2 - np.dot(AA2, theta_B_prime) * theta_B_prime)
    d_A3 = np.linalg.norm(AA3 - np.dot(AA3, theta_B_prime) * theta_B_prime)
    d_A4 = np.linalg.norm(AA4 - np.dot(AA4, theta_B_prime) * theta_B_prime)
    d_As = np.array([d_A1, d_A2, d_A3, d_A4])
    d_A_max = np.max(d_As)

    BB1 = np.array([l_B / 2 * np.cos(h_B) - w_B / 2 * -np.sin(h_B), l_B / 2 * np.sin(h_B) - w_B / 2 * np.cos(h_B)])
    BB2 = np.array([l_B / 2 * np.cos(h_B) + w_B / 2 * -np.sin(h_B), l_B / 2 * np.sin(h_B) + w_B / 2 * np.cos(h_B)])
    BB3 = np.array([-l_B / 2 * np.cos(h_B) - w_B / 2 * -np.sin(h_B), -l_B / 2 * np.sin(h_B) - w_B / 2 * np.cos(h_B)])
    BB4 = np.array([-l_B / 2 * np.cos(h_B) + w_B / 2 * -np.sin(h_B), -l_B / 2 * np.sin(h_B) + w_B / 2 * np.cos(h_B)])
    d_B1 = np.linalg.norm(BB1 - np.dot(BB1, theta_B_prime) * theta_B_prime)
    d_B2 = np.linalg.norm(BB2 - np.dot(BB2, theta_B_prime) * theta_B_prime)
    d_B3 = np.linalg.norm(BB3 - np.dot(BB3, theta_B_prime) * theta_B_prime)
    d_B4 = np.linalg.norm(BB4 - np.dot(BB4, theta_B_prime) * theta_B_prime)
    d_Bs = np.array([d_B1, d_B2, d_B3, d_B4])
    d_B_max = np.max(d_Bs)

    MFD = D_t1 - (d_A_max + d_B_max)
    InDepth = D_SAFE - MFD

    return InDepth


# 计算ACT、EI数值
def compute_real_time_metrics(x_A, y_A, v_A, h_A, l_A, w_A, x_B, y_B, v_B, h_B, l_B, w_B):

    v_Br = compute_v_Br(x_A, y_A, v_A, h_A, x_B, y_B, v_B, h_B)
    collision_result = get_collision_A_and_B(x_A, y_A, x_B, y_B, h_A, h_B, l_A, w_A, l_B, w_B)

    if not collision_result:

        shortest_distance, closest_points, v_closest = (
            compute_shortest_distance(x_A, y_A, v_A, h_A, l_A, w_A, x_B, y_B, v_B, h_B, l_B, w_B))

        RTTC, DTC, v_norm = compute_RTTC(x_A, y_A, v_A, h_A, l_A, w_A, x_B, y_B, v_B, h_B, l_B, w_B)

        if v_closest > 0:
            InDepth= compute_InDepth(x_A, y_A, v_A, h_A, l_A, w_A, x_B, y_B, v_B, h_B, l_B, w_B)
            InDepth = np.nan if InDepth is None else InDepth
            if InDepth >= 0:
                ACT = shortest_distance / v_closest
                MEI = InDepth / RTTC if not np.isnan(RTTC) and RTTC != 0 else np.nan

            else:
                MEI = ACT = np.nan
        else:
            InDepth = MEI = ACT = np.nan

    else:
         InDepth = MEI = ACT = v_closest = shortest_distance = RTTC = DTC = v_norm = np.nan


    return ACT, v_closest, shortest_distance, InDepth, MEI, RTTC, DTC, v_norm


# 定义常量
D_SAFE = 0  # EI的安全区域参数，暂时默认为0（不考虑安全冗余）


# 主函数
def main():
    # 自车A参数实时参数传入（示例）
    x_A = 0  # 绝对坐标，单位是m
    y_A = 0
    v_A = 0.1  # 单位是m/s
    h_A = 0  # 航向角，弧度制，范围是[-pi, pi]，例如1.57是90°
    l_A = 10  # 车长
    w_A = 2.5  # 车宽

    # 周车B参数实时参数传入（示例）
    x_B = -2
    y_B = 8
    v_B = 5
    h_B = -1
    l_B = 4.8
    w_B = 1.8

    # 包括的SSMs有：(1)ACT; (2)RTTC(2DTTC); (3)MEI
    # 计算SSMs数值
    ACT, v_closest, Shortest_D, InDepth, MEI, RTTC, DTC, v_norm = compute_real_time_metrics(x_A, y_A, v_A, h_A, l_A, w_A, x_B, y_B, v_B, h_B, l_B, w_B)

    # 输出结果
    # (1) ACT,【越小】越危险，越接近0越危险，取值范围是[0,+∞]，取危险场景可以筛选ACT小于3s的所有事件。
    print(f"ACT: {ACT:.4f} s")
    print(f"v_closest: {v_closest:.4f} m/s")  # ACT计算过程量
    print(f"Shortest_D: {Shortest_D:.4f} m")  # ACT计算过程量

    # (2) RTTC(2DTTC),【越小】越危险，越接近0越危险，取值范围是[0,+∞]，取危险场景可以筛选RTTC小于3s的所有事件。
    print(f"RTTC: {RTTC:.4f} s")
    print(f"DTC: {DTC:.4f} m")  # RTTC计算过程量
    print(f"v_norm: {v_norm:.4f} m/s")  # RTTC计算过程量

    # (3) MEI，【越大】越危险，取值范围是[0,+∞]，取危险场景可以筛选EI大于3 [m/s]的所有事件。
    print(f"InDepth: {InDepth:.4f} m")  # EI、MEI计算过程量
    print(f"MEI: {MEI:.4f} m/s")


if __name__ == "__main__":
    main()
