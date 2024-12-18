import numpy as np
import xarray as xr
import pandas as pd
import os
from glob import glob

def load_obs(csvpath:str, start:str, end:str) -> pd.DataFrame:
    '''
    load observation data from csv files in csvpath
    csvpath: path to the csv files
    start: start date in the format 'YYYY-MM-DD'
    end: end date in the format 'YYYY-MM-DD'
    '''
    
    dates = pd.date_range(start, end, freq='D')
    dfdates = {}
    for date in dates:
        yy = date.strftime('%Y')
        mm = date.strftime('%m')
        dd = date.strftime('%d')
        file = csvpath + '/china_sites_' + yy + mm + dd + '.csv'
        try:
            dfdates[str(yy+mm+dd)] = pd.read_csv(file)
        except:
            print(f'Cannot find file: {file}')
        
    obs = pd.concat(dfdates, ignore_index=True)
    
    return obs

def load_location(locpath:str, extend:list =None) -> pd.DataFrame:
    '''
    load site locations from the latest file in locpath
    locpath: path to the location file
    extend: [lon_min, lon_max, lat_min, lat_max]
    '''
    
    try:
        sitelocations = pd.read_excel(locpath + '/sitelocations_from2022.02.13.xlsx')
    except FileNotFoundError:
        xlsx_files = glob(os.path.join(locpath, '站点列表*.xlsx'))
        csv_files = glob(os.path.join(locpath, '站点列表*.csv'))
        
        if xlsx_files:
            latest_file = max(xlsx_files, key=os.path.getctime)
            sitelocations = pd.read_excel(latest_file)
        elif csv_files:
            latest_file = max(csv_files, key=os.path.getctime)
            sitelocations = pd.read_csv(latest_file)
        else:
            raise FileNotFoundError(f'No suitable file found in {locpath}.')
        
    # 去掉经纬度不是数字的点
    sitelocations = sitelocations[pd.to_numeric(sitelocations['经度'], errors='coerce').notnull()]
    sitelocations = sitelocations[pd.to_numeric(sitelocations['纬度'], errors='coerce').notnull()]
    
    if extend:
        lon_min = extend[0]
        lon_max = extend[1]
        lat_min = extend[2]
        lat_max = extend[3]
        sitelocations = sitelocations[(sitelocations['经度']>=lon_min) & (sitelocations['纬度']>=lat_min) & 
                                      (sitelocations['经度']<=lon_max) & (sitelocations['纬度']<=lat_max)]

    return sitelocations

def process_data(outpath:str, obs:pd.DataFrame, sitelocations:pd.DataFrame, variables: list=None) -> None:
    '''
    process observation data and save to excel files
    outpath: path to save the excel files
    obs: observation data
    sitelocations: site locations
    variables: variables to save, example: ['O3','PM2.5','PM10']
    '''
    
    sitelist = sitelocations['监测点编码'].tolist()
    
    df = obs[['date','hour','type']+sitelist]
    df['datetime'] = pd.to_datetime(df['date'].astype(str) + 'T' + df['hour'].astype(str) + ':00')
    df.set_index('datetime', inplace=True)
    df.drop(columns=['date','hour'], inplace=True)
    
    with pd.ExcelWriter(outpath + '/obs_sites.xlsx') as writer:
        if variables is None:
            variables = df['type'].unique()
        
        dfs = {}
        for var in variables:
            dfs[var] = df.groupby(['type']).get_group(var)
            dfs[var].drop(columns=['type'], inplace=True)
            dfs[var].to_excel(writer, sheet_name=var, index=True)
            
    return None

if __name__ == '__main__':
    csvpath = '/your/directory/obs_files'
    locpath = '/your/directory/site_locations'
    outpath = '/your/directory/outpath'
    start = '2024-04-11'
    end = '2024-04-13'
    extend = [100, 120, 30, 40]
    variables = ['O3','PM2.5','PM10']
    
    obs = load_obs(csvpath, start, end)
    sitelocations = load_location(locpath, extend)
    process_data(outpath, obs, sitelocations, variables)