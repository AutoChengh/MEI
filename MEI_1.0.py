# This version integrates: (1) collision detection; (2) SSM values are calculated only when the closest points between two vehicles are approaching
# Included SSMs: (1) ACT; (2) RTTC (2DTTC); (3) MEI
import numpy as np
import math


# Collision detection helper function 1
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

# Collision detection helper function 2
def check_collisions_between_series(A_series, B_series, axes, max_and_min_projections_1, max_and_min_projections_2):    # A_series contains the central coordinates of the vehicle
    # Calculate projections of A_series in parallel using numpy
    proj_A_min_max = []
    proj_B_min_max = []
    for i in range(4):
        proj_A = np.dot(A_series[:, :2], axes[i])
        proj_A_min_max.append([proj_A + max_and_min_projections_1[i][0], proj_A + max_and_min_projections_1[i][1]])
        proj_B = np.dot(B_series[:, :2], axes[i])
        proj_B_min_max.append([proj_B + max_and_min_projections_2[i][0], proj_B + max_and_min_projections_2[i][1]])
    proj_A_min_max = np.array(proj_A_min_max)
    proj_B_min_max = np.array(proj_B_min_max)    # Check if projections overlap
    # For each pair of vehicles in A and B, check for overlaps along 4 axes
    # Calculate in parallel using numpy
    if_collision = []
    for i in range(4):
        if_collision.append(np.logical_not(
            (proj_A_min_max[i][1][:, np.newaxis] < proj_B_min_max[i][0][np.newaxis, :]) | (
                        proj_B_min_max[i][1][np.newaxis, :] < proj_A_min_max[i][0][:, np.newaxis])))
    if_collision = np.array(if_collision)
    if_collision = np.all(if_collision, axis=0)
    return if_collision

# Collision detection helper function 3
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



# RTTC calculation helper function 1
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

# RTTC calculation helper function 2
def compute_RTTC(x_A, y_A, v_A, h_A, l_A, w_A, x_B, y_B, v_B, h_B, l_B, w_B):
    rotate_matrix_A = np.array([[np.cos(h_A), np.sin(h_A)], [-np.sin(h_A), np.cos(h_A)]])
    rotate_matrix_B = np.array([[np.cos(h_B), np.sin(h_B)], [-np.sin(h_B), np.cos(h_B)]])
    #Modified on 2025/02/21
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



# ACT calculation helper function 1
def get_rect_corners(x, y, h, l, w):
    # Coordinates of four corners of the rectangle relative to the center point
    corners = [
        [-l / 2 * math.cos(h) - w / 2 * math.sin(h), -l / 2 * math.sin(h) + w / 2 * math.cos(h)],  # A1
        [l / 2 * math.cos(h) - w / 2 * math.sin(h), l / 2 * math.sin(h) + w / 2 * math.cos(h)],  # A2
        [l / 2 * math.cos(h) + w / 2 * math.sin(h), l / 2 * math.sin(h) - w / 2 * math.cos(h)],  # A3
        [-l / 2 * math.cos(h) + w / 2 * math.sin(h), -l / 2 * math.sin(h) - w / 2 * math.cos(h)],  # A4
    ]
    # Calculate absolute coordinates
    corners = [[x + corner[0], y + corner[1]] for corner in corners]
    return corners

# ACT calculation helper function 2: Calculate distance between two points
def distance(p1, p2):
    return math.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)

# ACT calculation helper function 3: Calculate shortest distance from point to line segment
def point_to_segment_distance(p, v1, v2):    # Vector calculation
    line_len = distance(v1, v2)
    if line_len == 0:  # Avoid division by zero
        return distance(p, v1)    # Calculate projection of point onto line segment
    t = max(0, min(1, ((p[0] - v1[0]) * (v2[0] - v1[0]) + (p[1] - v1[1]) * (v2[1] - v1[1])) / line_len ** 2))
    # Calculate closest point
    closest = [v1[0] + t * (v2[0] - v1[0]), v1[1] + t * (v2[1] - v1[1])]
    return distance(p, closest), closest

# ACT calculation helper function 4: Calculate shortest distance between rectangles
def get_shortest_distance(corners_A, corners_B):
    min_distance = float('inf')
    closest_points = None
    closest_A = None
    closest_B = None
    # Calculate shortest distance from each corner to the four edges of the other rectangle
    for i in range(4):
        for j in range(4):
            # From corner of A to edges of B
            p1 = corners_A[i]
            for k in range(4):
                p2 = corners_B[k]
                p3 = corners_B[(k + 1) % 4]
                dist, closest = point_to_segment_distance(p1, p2, p3)
                if dist < min_distance:
                    min_distance = dist
                    closest_A = p1
                    closest_B = closest
            # From corner of B to edges of A
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

# ACT calculation helper function 5
def compute_shortest_distance(x_A, y_A, v_A, h_A, l_A, w_A, x_B, y_B, v_B, h_B, l_B, w_B):
    # Calculate four corner points of two rectangles
    corners_A = get_rect_corners(x_A, y_A, h_A, l_A, w_A)
    corners_B = get_rect_corners(x_B, y_B, h_B, l_B, w_B)    # Get collision detection results
    # Calculate shortest distance
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



# EI calculation helper function 1
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

# EI calculation helper function 2
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


# Calculate ACT and EI values
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


# Define constants
D_SAFE = 0  # Safety zone parameter for EI, temporarily set to 0 (no safety redundancy considered)


# Main function
def main():    
    # Ego vehicle (Vehicle A) parameters input in real-time (example)
    x_A = 0  # Absolute coordinate, unit: m
    y_A = 0
    v_A = 0.1  # Unit: m/s
    h_A = 0  # Heading angle, in radians, range [-pi, pi], e.g., 1.57 is 90°
    l_A = 10  # Vehicle length
    w_A = 2.5  # Vehicle width    
    # Surrounding vehicle (Vehicle B) parameters input in real-time (example)
    x_B = -2
    y_B = 8
    v_B = 5
    h_B = -1
    l_B = 4.8
    w_B = 1.8    
    
    # Included SSMs: (1)ACT; (2)RTTC(2DTTC); (3)MEI
    # Calculate SSM values
    ACT, v_closest, Shortest_D, InDepth, MEI, RTTC, DTC, v_norm = compute_real_time_metrics(x_A, y_A, v_A, h_A, l_A, w_A, x_B, y_B, v_B, h_B, l_B, w_B)    # Output results
    
    # (1) ACT: Lower values indicate higher danger, especially as it approaches 0. Range is [0,+∞]. 
    #     For dangerous scenarios, consider filtering events where ACT is less than 3s.
    print(f"ACT: {ACT:.4f} s")    
    print(f"v_closest: {v_closest:.4f} m/s")  # Intermediate value for ACT calculation
    print(f"Shortest_D: {Shortest_D:.4f} m")  # Intermediate value for ACT calculation

    # (2) RTTC(2DTTC): Lower values indicate higher danger, especially as it approaches 0. Range is [0,+∞].
    #     For dangerous scenarios, consider filtering events where RTTC is less than 3s.
    print(f"RTTC: {RTTC:.4f} s")
    print(f"DTC: {DTC:.4f} m")  # Intermediate value for RTTC calculation
    print(f"v_norm: {v_norm:.4f} m/s")  # Intermediate value for RTTC calculation

    # (3) MEI: Higher values indicate higher danger. Range is [0,+∞].
    #     For dangerous scenarios, consider filtering events where MEI is greater than 3 [m/s].
    print(f"InDepth: {InDepth:.4f} m")  # Intermediate value for EI and MEI calculation
    print(f"MEI: {MEI:.4f} m/s")


if __name__ == "__main__":
    main()
