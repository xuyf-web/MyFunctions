import numpy as np
from math import radians, sin, cos, sqrt, asin

def nearest_position(stn_lon, stn_lat, lon2d, lat2d):
    """
    获取最临近格点坐标索引
    """
    
    # 计算差异矩阵
    difflon = stn_lon - lon2d
    difflat = stn_lat - lat2d
    
    # 计算距离平方
    distances_squared = np.square(difflon) + np.square(difflat)
    
    # 获取最小距离的索引
    min_index = np.unravel_index(np.argmin(distances_squared), distances_squared.shape)
    
    return min_index

def nearest_positions(stn_lon, stn_lat, lon2d, lat2d, num=4):
    """
    获取最临近的指定数量num个格点坐标索引
    """
    
    # 计算差异矩阵
    difflon = stn_lon - lon2d
    difflat = stn_lat - lat2d
    
    # 计算距离平方
    distances_squared = np.square(difflon) + np.square(difflat)
    
    # 获取最小的指定数量num个距离的索引
    flat_indices = np.argsort(distances_squared.ravel())[:num]
    indices = np.unravel_index(flat_indices, distances_squared.shape)
    
    return indices

def weighted_average(stn_lon, stn_lat, lon2d, lat2d, data2d, num=4):
    """
    计算近邻格点的加权平均值
    """
    
    # 获取最临近的num个格点坐标索引
    indices = nearest_positions(stn_lon, stn_lat, lon2d, lat2d, num=num)
    
    # 提取num个格点的数据
    flat_indices = np.ravel_multi_index(indices, data2d.shape)
    points_data = data2d.flat[flat_indices]

    # 计算num个格点到站点的距离
    difflon = stn_lon - lon2d.flat[flat_indices]
    difflat = stn_lat - lat2d.flat[flat_indices]
    distances_squared = np.square(difflat) + np.square(difflon)

    # 计算权重（距离的倒数可以作为权重）
    weights = 1 / (distances_squared + 1e-10)  # 避免除零错误
    weights /= np.sum(weights)  # 归一化权重

    # 加权平均
    weighted_avg = np.sum(points_data * weights)

    return weighted_avg

def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified with longitude and latitude in decimal degrees)
    """
    # 将十进制度数转换为弧度
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

    # 地球半径（千米）
    R = 6371.0

    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a))
    
    return R * c


if __name__ == '__main__':
    # 生成测试数据
    lon2d, lat2d = np.meshgrid(np.linspace(-180, 180, 361), np.linspace(-90, 90, 181))
    data2d = np.random.rand(181, 361)
    stn_lon = 116.3971
    stn_lat = 39.9074
    
    # 提取最临近格点坐标索引
    min_index = nearest_position(stn_lon, stn_lat, lon2d, lat2d)
    print('nearest index = ', min_index)
    
    # 提取最临近的指定数量num个格点坐标索引
    indices = nearest_positions(stn_lon, stn_lat, lon2d, lat2d, num=4)
    print('nearest indices = ', indices)
    
    # 计算加权平均
    weighted_avg = weighted_average(stn_lon, stn_lat, lon2d, lat2d, data2d, num=4)
    print(weighted_avg)
    
    # 计算最近点的距离
    nearest_lon = lon2d[min_index]
    nearest_lat = lat2d[min_index]
    distance = haversine(stn_lon, stn_lat, nearest_lon, nearest_lat)
    print(f"Distance to nearest point: {distance:.2f} km")