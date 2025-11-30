import cv2
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors


def analyze_red_square_matrix(image_path):
    """
    完整分析红色正方形方阵并可视化结果
    """
    image = cv2.imread(image_path)
    if image is None:
        return None

    image_rgb = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
    original_image = image.copy()

    # HSV颜色空间检测红色
    hsv = cv2.cvtColor(image, cv2.COLOR_BGR2HSV)

    lower_red1 = np.array([0, 120, 70])
    upper_red1 = np.array([10, 255, 255])
    lower_red2 = np.array([170, 120, 70])
    upper_red2 = np.array([180, 255, 255])

    mask1 = cv2.inRange(hsv, lower_red1, upper_red1)
    mask2 = cv2.inRange(hsv, lower_red2, upper_red2)
    red_mask = mask1 + mask2

    # 形态学操作清理噪声;
    kernel = np.ones((3, 3), np.uint8)
    red_mask = cv2.morphologyEx(red_mask, cv2.MORPH_CLOSE, kernel)
    red_mask = cv2.morphologyEx(red_mask, cv2.MORPH_OPEN, kernel)

    # 查找掩码中的轮廓;
    contours, _ = cv2.findContours(red_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)

    # 遍历轮廓, 筛选出正方形区域;
    squares_info = []
    for i, contour in enumerate(contours):
        area = cv2.contourArea(contour)
        if area < 30:
            continue

        x, y, w, h = cv2.boundingRect(contour)
        aspect_ratio = w / float(h)

        if 0.2 <= aspect_ratio <= 2:
            center_x = x + w // 2
            center_y = y + h // 2
            squares_info.append({
                'id': i,
                'center': (center_x, center_y),
                'bbox': (x, y, w, h),
                'contour': contour,
                'area': area
            })
            # print(i, ": center: ", (center_x, center_y))

    if not squares_info:
        print("未检测到红色正方形")
        return None

    # ========= 用 DBSCAN 过滤掉边缘零星红点，只保留主阵列 =========
    coords = np.array([s['center'] for s in squares_info], dtype=np.float32)

    if len(coords) >= 10:
        # 估计 spot 间距：使用最近邻距离的中位数
        nbrs = NearestNeighbors(n_neighbors=2).fit(coords)
        distances, _ = nbrs.kneighbors(coords)
        # 每个点到最近邻的距离（忽略自身距离 0）
        nn_dists = distances[:, 1]
        typical_spacing = np.median(nn_dists)

        # eps 设置为典型间距的 1.5 倍（可以视情况微调）
        eps = typical_spacing * 1.5
        # 至少有多少个邻居算一个“密集点”，阵列比较密，可以用 4~6
        min_samples = 2

        db = DBSCAN(eps=eps, min_samples=min_samples).fit(coords)
        labels = db.labels_

        # 去掉噪声点（label = -1），只在剩余 cluster 里选“最大簇”
        valid = labels != -1
        if np.any(valid):
            valid_labels = labels[valid]
            # 统计每个簇的大小，选出点数最多的那个簇
            counts = np.bincount(valid_labels)
            main_label = np.argmax(counts)

            filtered_squares = [
                s for s, label in zip(squares_info, labels)
                if label == main_label
            ]

            # 确保过滤后仍然有足够多的点
            if len(filtered_squares) >= 10:
                squares_info = filtered_squares
                print(f"使用 DBSCAN 过滤后保留 {len(squares_info)} 个正方形（主阵列）")
            else:
                print("DBSCAN 过滤后点数太少，退回使用全部点")
        else:
            print("DBSCAN 没有找到有效簇，退回使用全部点")
    else:
        print("点太少，跳过 DBSCAN 过滤")

    # ---------- Corner detection: region + distance to image corner ----------
    # 图会有一点点倾斜, 所以这样选择会好一些, 而不是直接找x和y轴最小值;
    h, w = image.shape[:2]

    # helper: squared distance
    def dist2(p, q):
        return (p[0] - q[0]) ** 2 + (p[1] - q[1]) ** 2

    # region masks
    def in_top_left(c):    return c[0] < w * 0.3 and c[1] < h * 0.3
    def in_top_right(c):   return c[0] > w * 0.7 and c[1] < h * 0.3
    def in_bottom_left(c): return c[0] < w * 0.3 and c[1] > h * 0.7

    tl_candidates = [s for s in squares_info if in_top_left(s['center'])]
    tr_candidates = [s for s in squares_info if in_top_right(s['center'])]
    bl_candidates = [s for s in squares_info if in_bottom_left(s['center'])]

    # 若某个角区域没有候选，就退回到全局最接近该角点的方块
    if not tl_candidates:
        tl_candidates = squares_info
    if not tr_candidates:
        tr_candidates = squares_info
    if not bl_candidates:
        bl_candidates = squares_info

    # image corners: (x, y)
    tl_corner = (0, 0)
    tr_corner = (w - 1, 0)
    bl_corner = (0, h - 1)

    upper_left_square = min(tl_candidates,
                            key=lambda s: dist2(s['center'], tl_corner))
    upper_right_square = min(tr_candidates,
                             key=lambda s: dist2(s['center'], tr_corner))
    bottom_left_square = min(bl_candidates,
                             key=lambda s: dist2(s['center'], bl_corner))

    upper_left_center = upper_left_square['center']
    upper_right_center = upper_right_square['center']
    bottom_left_center = bottom_left_square['center'] 

    print(f"检测到 {len(squares_info)} 个红色正方形")
    print(f"左上角正方形中心坐标: {upper_left_center}")
    print(f"左上角正方形边界框: {upper_left_square['bbox']}")
    print(f"左上角正方形面积: {upper_left_square['area']:.2f}")

    print(f"右上角正方形中心坐标: {upper_right_center}")
    print(f"右上角正方形边界框: {upper_right_square['bbox']}")
    print(f"右上角正方形面积: {upper_right_square['area']:.2f}")

    print(f"左下角正方形中心坐标: {bottom_left_center}")
    print(f"左下角正方形边界框: {bottom_left_square['bbox']}")
    print(f"左下角正方形面积: {bottom_left_square['area']:.2f}")

    # --------------------- Visualization ---------------------
    # 可视化结果
    plt.figure(figsize=(15, 10))

    # 原始图像
    plt.subplot(2, 3, 1)
    plt.imshow(image_rgb)
    plt.title('original image')
    plt.axis('off')

    # 红色掩码
    plt.subplot(2, 3, 2)
    plt.imshow(red_mask, cmap='gray')
    plt.title('Red zone detection')
    plt.axis('off')

    # 所有检测到的正方形
    plt.subplot(2, 3, 3)
    result_image = image_rgb.copy()
    for square in squares_info:
        x, y, w, h = square['bbox']
        cv2.rectangle(result_image, (x, y), (x + w, y + h), (0, 255, 0), 2)
        cv2.circle(result_image, square['center'], 3, (255, 0, 0), -1)
    plt.imshow(result_image)
    plt.title(f'Detected {len(squares_info)} red squares')
    plt.axis('off')

    # 左上角正方形特写
    plt.subplot(2, 3, 4)
    x, y, w, h = upper_left_square['bbox']
    # 扩大区域以便更好观察
    margin = 10
    y_start = max(0, y - margin)
    y_end = min(image.shape[0], y + h + margin)
    x_start = max(0, x - margin)
    x_end = min(image.shape[1], x + w + margin)
    roi = image_rgb[y_start:y_end, x_start:x_end]
    plt.imshow(roi)
    plt.title('The square area in the upper left corner')
    plt.axis('off')

    # 标记左上角中心点
    plt.subplot(2, 3, 5)
    marked_image = image_rgb.copy()
    # 绘制所有正方形中心
    for square in squares_info:
        cv2.circle(marked_image, square['center'], 3, (0, 255, 0), -1)
    # 突出显示左上角中心
    cv2.circle(marked_image, upper_left_center, 8, (255, 0, 0), -1)
    cv2.putText(marked_image, f"({upper_left_center[0]}, {upper_left_center[1]})",
                (upper_left_center[0] + 10, upper_left_center[1] - 10),
                cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255, 0, 0), 2)
    plt.imshow(marked_image)
    plt.title('The central point mark in the upper left corner')
    plt.axis('off')

    # 轮廓检测结果
    plt.subplot(2, 3, 6)
    contour_image = np.zeros_like(image_rgb)
    cv2.drawContours(contour_image, [upper_left_square['contour']], -1, (255, 0, 0), 2)
    # 绘制最小外接矩形
    rect = cv2.minAreaRect(upper_left_square['contour'])
    box = cv2.boxPoints(rect)
    box = box.astype(np.int32)
    cv2.drawContours(contour_image, [box], 0, (0, 255, 0), 2)
    # 绘制中心点
    cv2.circle(contour_image, upper_left_center, 5, (0, 0, 255), -1)
    plt.imshow(contour_image)
    plt.title('Outline and the smallest enclosing rectangle')
    plt.axis('off')

    plt.tight_layout()
    plt.show()

    
    return upper_left_center


# 使用示例
# result = analyze_red_square_matrix("/data/workdir/panw/py_singlecell/Spatial/MAGIC-seq/data/T9-sample1-7/T9-sample1-7-spot.png")

# result = analyze_red_square_matrix("/data/workdir/panw/py_singlecell/Spatial/MAGIC-seq/data/T9-sample1-1-spot.png")
result = analyze_red_square_matrix("/data/workdir/panw/py_singlecell/Spatial/MAGIC-seq/data/T9-sample1-2-spot.png")
# result = analyze_red_square_matrix("/data/workdir/panw/py_singlecell/Spatial/MAGIC-seq/data/T9-sample1-3-spot.png")
