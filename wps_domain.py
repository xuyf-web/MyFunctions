import re
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
from shapely.geometry import Polygon
from pyproj import Proj, CRS, Transformer

def get_param(pattern, content):
    match = re.search(pattern, content)
    if match:
        return float(match.group(1))
    else:
        raise ValueError(f'Parameter {pattern} not found in namelist.wps')

def get_params(pattern, content):
    match = re.search(pattern, content)
    if match:
        return [int(x.strip()) for x in match.group(1).split(',') if x.strip()]
    else:
        raise ValueError(f'Parameter {pattern} not found in namelist.wps')
    
def parse_namelist(namelist_path):
    with open(namelist_path, 'r') as file:
        namelist_content = file.read()
        params = {
            'max_dom'          : int(get_param(r'max_dom\s*=\s*(\d+)', namelist_content)),
            'parent_id'        : get_params(r'parent_id\s*=\s*([\d,\s]+)', namelist_content),
            'parent_grid_ratio': get_params(r'parent_grid_ratio\s*=\s*([\d,\s]+)', namelist_content),
            'i_parent_start'   : get_params(r'i_parent_start\s*=\s*([\d,\s]+)', namelist_content),
            'j_parent_start'   : get_params(r'j_parent_start\s*=\s*([\d,\s]+)', namelist_content),
            'e_we'             : get_params(r'e_we\s*=\s*([\d,\s]+)', namelist_content),
            'e_sn'             : get_params(r'e_sn\s*=\s*([\d,\s]+)', namelist_content),
            'dx'               : get_param(r'dx\s*=\s*(\d+)', namelist_content),
            'dy'               : get_param(r'dy\s*=\s*(\d+)', namelist_content),
            'ref_lat'          : get_param(r'ref_lat\s*=\s*([-+]?\d*\.\d+|\d+)', namelist_content),
            'ref_lon'          : get_param(r'ref_lon\s*=\s*([-+]?\d*\.\d+|\d+)', namelist_content),
            'truelat1'         : get_param(r'truelat1\s*=\s*([-+]?\d*\.\d+|\d+)', namelist_content),
            'truelat2'         : get_param(r'truelat2\s*=\s*([-+]?\d*\.\d+|\d+)', namelist_content),
            'stand_lon'        : get_param(r'stand_lon\s*=\s*([-+]?\d*\.\d+|\d+)', namelist_content)
        }
    return params

def init_projection(params):
    # 检查参数是否完整
    if not all(key in params for key in ['ref_lat', 'ref_lon', 'truelat1', 'truelat2', 'stand_lon']):
        raise ValueError('Missing required parameters in namelist.wps')
    
    # 定义平面投影参数
    proj_params = {
        'proj': 'lcc',
        'lat_1': params['truelat1'],
        'lat_2': params['truelat2'],
        'lat_0': params['ref_lat'],
        'lon_0': params['ref_lon'],
        'lon_1': params['stand_lon'],
        'x_0': 0,
        'y_0': 0,
        'datum': 'WGS84'
    }

    # 创建平面投影坐标系
    projected_coordinate = Proj(proj_params)
    
    return projected_coordinate

def compute_corner(x_left, y_bottom, x_right, y_top, projection, dx, dy):
    # 定义坐标点
    corners = {
        'west_south': (x_left, y_bottom),
        'east_north': (x_right, y_top),
        'west_north': (x_left, y_top),
        'east_south': (x_right, y_bottom),
    }
    
    # 使用投影计算每个点的经纬度
    results = {}
    for key, (x, y) in corners.items():
        lon, lat = projection(x, y, inverse=True)
        results[f'lon_{key}'] = lon
        results[f'lat_{key}'] = lat
    
    # 添加额外参数
    results.update({
        'dx': dx,
        'dy': dy,
        'x_left': x_left,
        'y_bottom': y_bottom,
    })

    return results

def init_domains(params):
    # 初始化投影坐标系
    projection = init_projection(params)

    # 获取主网格中心点投影坐标
    ref_lat, ref_lon = float(params['ref_lat']), float(params['ref_lon'])
    x_center, y_center = projection(ref_lon, ref_lat)
    
    # 初始化网格列表
    domains = []

    # 遍历每个网格层级
    for i in range(params['max_dom']):
        parent_id = int(params['parent_id'][i])
        grid_ratio = int(params['parent_grid_ratio'][i])
        e_we = int(params['e_we'][i])
        e_sn = int(params['e_sn'][i])
        dx = float(params['dx']) / (grid_ratio if i > 0 else 1)
        dy = float(params['dy']) / (grid_ratio if i > 0 else 1)

        # 如果是主网格（d01）
        if i == 0:
            x_left = x_center - ((e_we - 1) / 2) * dx
            y_bottom = y_center - ((e_sn - 1) / 2) * dy
            x_right = x_center + ((e_we - 1) / 2) * dx
            y_top = y_center + ((e_sn - 1) / 2) * dy
        else:
            # 子网格基于父网格的位置
            parent_domain = domains[parent_id - 1] # todo
            parent_dx = parent_domain['dx']
            parent_dy = parent_domain['dy']
            x_left = parent_domain['x_left'] + (params['i_parent_start'][i] - 1) * parent_dx
            y_bottom = parent_domain['y_bottom'] + (params['j_parent_start'][i] - 1) * parent_dy
            x_right = x_left + (e_we - 1) * dx
            y_top = y_bottom + (e_sn - 1) * dy

        # 计算当前网格边界
        domain = compute_corner(x_left, y_bottom, x_right, y_top, projection, dx, dy)
        domains.append(domain)

        print(f"Domain {i+1}: ")
        print(f"  South_West = ({domain['lon_west_south']:.2f}, {domain['lat_west_south']:.2f})")
        print(f"  North_West = ({domain['lon_west_north']:.2f}, {domain['lat_west_north']:.2f})")
        print(f"  North_East = ({domain['lon_east_north']:.2f}, {domain['lat_east_north']:.2f})")
        print(f"  South_East = ({domain['lon_east_south']:.2f}, {domain['lat_east_south']:.2f})")
        print('-----------------------------------')

    return domains

def add_gridline(ax, labelsize, linewidth=0.5, linestyle='--', color='gray'):
    gl = ax.gridlines(
            draw_labels=True, x_inline=False, y_inline=False,
            linewidth=linewidth, linestyle=linestyle, color=color)
    gl.top_labels = False
    gl.right_labels = False
    gl.rotate_labels = False
    gl.xlabel_style = {'size': labelsize}
    gl.ylabel_style = {'size': labelsize}
    
def add_labels(ax, i, domain, colors,x_offset=0.1, y_offset=0.1):
    label_lon, label_lat = domain['lon_west_north'], domain['lat_west_north']
    ax.text(
        label_lon+x_offset, label_lat-y_offset, 
        f'd0{i+1}', 
        color=colors[i % len(colors)],
        fontsize=16,
        fontweight='bold',
        ha='left', va='top', 
        transform=ccrs.PlateCarree()
    )

def plot_domains(params, domains, coastline=True, show_labels=True, extent=None):
    # 定义 Lambert 投影和等经纬度投影
    proj_lambert = CRS.from_proj4(
        f"+proj=lcc +lat_1={params['truelat1']} +lat_2={params['truelat2']} "
        f"+lat_0={params['ref_lat']} +lon_0={params['ref_lon']} +datum=WGS84"
        )
    proj_plate = CRS.from_proj4("+proj=latlong +datum=WGS84")
    
    # 使用 Transformer 代替 transform
    transformer = Transformer.from_proj(proj_plate, proj_lambert)

    # 初始化 Lambert Conformal 投影
    proj = ccrs.LambertConformal(
        central_longitude=params['ref_lon'],
        central_latitude=params['ref_lat'],
        standard_parallels=(params['truelat1'], params['truelat2'])
    )
    fig, ax = plt.subplots(figsize=(12, 8), subplot_kw={'projection': proj}, dpi=300)

    if coastline:
        ax.coastlines(resolution='10m', color='black', linewidth=1)
        ax.add_feature(cfeature.LAND, facecolor='lightgray')
        ax.add_feature(cfeature.OCEAN, facecolor='lightblue')
        # ax.add_feature(cfeature.LAKES, facecolor='lightblue')
    
    add_gridline(ax, labelsize=16)
    
    # 定义颜色列表，用于区分网格层级
    colors = ['red', 'blue', 'green', 'purple', 'orange']
    
    # 绘制每个网格的范围
    for i, domain in enumerate(domains):
        # 将所有边界点转换到 Lambert 投影的经纬度，并创建多边形
        corners = ['west_south', 'west_north', 'east_north', 'east_south']
        polygon_lambert = Polygon([
            transformer.transform(domain[f'lon_{corner}'], domain[f'lat_{corner}'])
            for corner in corners
        ])

        # 绘制多边形
        ax.add_geometries(
            [polygon_lambert], 
            crs=proj,
            edgecolor=colors[i % len(colors)], 
            facecolor='none', 
            linewidth=2, 
            label=f'Domain {i+1}'
        )
        
        # 添加标签
        if show_labels:
            add_labels(ax, i, domain, colors,x_offset=0.1, y_offset=0.1)

    # 设置绘图范围
    if extent:
        ax.set_extent(extent, crs=ccrs.PlateCarree())
    else:
        lon_min = min(domain['lon_west_south'] for domain in domains)
        lon_max = max(domain['lon_east_north'] for domain in domains)
        lat_min = min(domain['lat_west_south'] for domain in domains)
        lat_max = max(domain['lat_east_north'] for domain in domains)

        lon_margin = (lon_max - lon_min) * 0.05
        lat_margin = (lat_max - lat_min) * 0.05
        ax.set_extent([lon_min - lon_margin, lon_max + lon_margin, 
                    lat_min - lat_margin, lat_max + lat_margin], crs=ccrs.PlateCarree())

    ax.set_title('WRF Nested Domains', fontsize=20)
    
    return fig, ax

if __name__ == '__main__':
    namelist_path = 'your_directory/namelist.wps'
    params = parse_namelist(namelist_path)
    domains = init_domains(params)
    
    # default coastline 10m
    fig, ax = plot_domains(params, domains, coastline=True)
    # with user-defined shapefile
    shapefile_path = 'your_directory/boundary.shp'
    myshapefile = cfeature.ShapelyFeature(Reader(shapefile_path).geometries(),
                ccrs.PlateCarree(), edgecolor='k', facecolor='lightblue')
    fig, ax = plot_domains(params, domains, coastline=True, extent=None)
    ax.add_feature(cfeature.BORDERS, linestyle='-', linewidth=1)
    ax.add_feature(myshapefile)
    
    plt.show()