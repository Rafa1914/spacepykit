import string
from astropy.io import fits
from src.astrophotometry import find_peaks
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from src.astrophotometry import gaia_radecs
from src.astrophotometry.geometry import sparsify
from src.astrophotometry import compute_wcs
from src.astride import Streak
import pandas as pd
import datetime
import os

from astropy.wcs import WCS

def get_wcs_from_fits(filename: string, fov_parameter: float = 1.2) -> WCS:
    '''
    Calculates the World Coordinate System for a FITS image using plate solving.

    Attributes:
    filename: FITS file's name
    fov_parameter: Multiplicator of FOV (Default is 1.2)
    '''

    # Open some FITS image
    image_name = filename
    hdu_list = fits.open(image_name) #Header Data Unit
    img_header = hdu_list[0].header
    img_data = hdu_list[0].data

    # Close the file
    hdu_list.close()

    # Obtaining the 20 brightest stars
    peaks_coordinates = find_peaks(img_data)[0:20]

    # Obtaining image center
    ra, dec = img_header['RA'], img_header['DEC']
    center_header = SkyCoord(ra,dec,unit='deg')

    # Obtaining field of view
    telescope_focal_length = 400 # mm
    pixel_size = 3.76 # um
    pixel_ratio = 206*pixel_size/telescope_focal_length # 0.66
    pixel = pixel_ratio * u.arcsec  # known pixel scale
    shape = img_data.shape
    fov = np.max(shape) * pixel.to(u.deg)

    # Obtained stars catalogued by GAIA
    all_radecs = gaia_radecs(center_header, fov_parameter * fov, circular = False)

    # we only keep stars 0.01 degree apart from each other
    all_radecs = sparsify(all_radecs, 0.01)

    # we only keep the 20 brightest stars from gaia
    wcs = compute_wcs(peaks_coordinates, all_radecs[0:20], tolerance=5)
    return wcs



def get_radec_from_fits(filename:string,output_path:string,wcs:WCS,contour_threshold: float = 3.0):
    '''
    Calculates the RA and DEC angles for a streak in a FITS image, given a WCS.

    Attributes:
    filename: File's name
    output_path: Path of destination file.
    wcs: WCS of the FITS image
    contour_threshold: Threshold to identify the streaks in the image.
    '''

    # Open the FITS image
    image_name = filename
    hdu_list = fits.open(image_name) #Header Data Unit
    img_header = hdu_list[0].header
    
    #Defining initial parameters
    connectivity_angle = 8.0

    # Read a fits image and create a Streak instance.
    streak = Streak(
        filename = filename, 
        contour_threshold = contour_threshold, 
        connectivity_angle = connectivity_angle, 
        )

    # Detect streaks
    streak.detect()

    if len(streak.streaks) > 0:
        # satellite_list = streak.streaks

        # Getting DataFrame data
        data = []
        for satellite in streak.streaks:
            # x_min = satellite['x_min']
            # y_min = satellite['y_min']
            # x_max = satellite['x_max']
            # y_max = satellite['y_max']
            # radec_min = wcs.pixel_to_world(x_min,y_min)
            # radec_max = wcs.pixel_to_world(x_max,y_max)
            # ra_min = radec_min.ra.deg
            # ra_max = radec_max.ra.deg
            # dec_min = radec_min.dec.deg
            # dec_max = radec_max.dec.deg
            # m = satellite['slope']
            # theta = satellite['slope_angle']
            # connectivity = satellite['connectivity']
            # n = satellite['intercept']
            # perimeter = satellite['perimeter']
            # shape_factor = satellite['shape_factor']

            new_row = {
                'Index': satellite['index'],
                'X_min':satellite['x_min'],
                'Y_min':satellite['y_min'],
                'X_max':satellite['x_max'],
                'Y_max':satellite['y_max'],
                'RA_min':wcs.pixel_to_world(satellite['x_min'],satellite['y_min']).ra.deg,
                'RA_max':wcs.pixel_to_world(satellite['x_max'],satellite['y_max']).ra.deg,
                'DEC_min':wcs.pixel_to_world(satellite['x_min'],satellite['y_min']).dec.deg,
                'DEC_max':wcs.pixel_to_world(satellite['x_max'],satellite['y_max']).dec.deg,
                'Coef. Angular':satellite['slope'],
                'Theta':satellite['slope_angle'],
                'Connectivity':satellite['connectivity'],
                'Interception': satellite['intercept'],
                'Perimeter': satellite['perimeter'],
                'Shape Factor': satellite['shape_factor'],
                'Custom Factor': satellite['perimeter']/satellite['shape_factor']
                }            
            # new_row = {
            #     'Index': satellite['index'],
            #     'X_min':x_min,
            #     'Y_min':y_min,
            #     'X_max':x_max,
            #     'Y_max':y_max,
            #     'RA_min':ra_min,
            #     'RA_max':ra_max,
            #     'DEC_min':dec_min,
            #     'DEC_max':dec_max,
            #     'Coef. Angular':m,
            #     'Theta':theta,
            #     'Connectivity':connectivity,
            #     'Interception': n,
            #     'Perimeter': perimeter,
            #     'Shape Factor': shape_factor,
            #     'Custom Factor': perimeter/shape_factor
            #     }            
            data.append(new_row)
        
        # Creating de DataFrame
        # columns = ['Index','X_min','Y_min','X_max','Y_max','RA_min','RA_max','DEC_min','DEC_max','Coef. Angular','Theta','Connectivity','Interception','Perimeter','Shape Factor','Custom Factor']
        columns = list(new_row.keys())

        satellites_df = pd.DataFrame(columns=columns,data = data)

        #Identificação da linha com melhor razão entre comprimento e área
        factor_max = (satellites_df['Custom Factor']).max()
        filter = satellites_df['Custom Factor'] == factor_max
        mvp_streak = satellites_df[filter]

        # ay+bx+c = 0
        b = mvp_streak['Coef. Angular']
        c = mvp_streak['Interception']
        a = -1

        for index, satellite in satellites_df[~filter].iterrows():
            tolerancia = 10.0 # Distância máxima em pixel
            x_center = (satellite['X_max'] + satellite['X_min']) * 0.5
            y_center = (satellite['Y_max'] + satellite['Y_min']) * 0.5
            # Distância
            d = np.abs(a*y_center+b*x_center+c)/np.sqrt(a**2+b**2)
            if d.values[0] < tolerancia:
                if satellite['Index'] > mvp_streak['Index'].values[0]:
                    mvp_streak['X_max'] = satellite['X_max']
                    mvp_streak['Y_max'] = satellite['Y_max']
                    mvp_streak['RA_max'] = satellite['RA_max']
                    mvp_streak['DEC_max'] = satellite['DEC_max']
                    mvp_streak['Perimeter'] = mvp_streak['Perimeter'] + satellite['Perimeter']
                else:
                    mvp_streak['X_min'] = satellite['X_min']
                    mvp_streak['Y_min'] = satellite['Y_min']
                    mvp_streak['RA_min'] = satellite['RA_min']
                    mvp_streak['DEC_min'] = satellite['DEC_min']
                    mvp_streak['Perimeter'] = mvp_streak['Perimeter'] + satellite['Perimeter']

        #Informações de Tempo
        initial_time = datetime.datetime.strptime(img_header['S_EXP'],'%Y-%m-%dT%H:%M:%S.%f')
        interval = img_header['EXPTIME']
        final_time = initial_time + datetime.timedelta(seconds = interval)

        # Resultados da análise
        columns_results = ['Time', 'RA[deg]', 'DEC[deg]']
        initial_obs = {
            'Time':initial_time.strftime('%Y-%m-%dT%H:%M:%S.%f'),
            'RA[deg]': mvp_streak['RA_max'].values[0],
            'DEC[deg]': mvp_streak['DEC_max'].values[0]
        }

        final_obs = {
            'Time':final_time.strftime('%Y-%m-%dT%H:%M:%S.%f'),
            'RA[deg]': mvp_streak['RA_min'].values[0],
            'DEC[deg]': mvp_streak['DEC_min'].values[0]
        }

        # Escrevendo em um arquivo
        with open(output_path,'a') as file:
            file.write(f'ANGLE_1 = {initial_obs["Time"]} {initial_obs["RA[deg]"]}\n')
            file.write(f'ANGLE_2 = {initial_obs["Time"]} {initial_obs["DEC[deg]"]}\n')
            file.write(f'\n')
            file.write(f'ANGLE_1 = {final_obs["Time"]} {final_obs["RA[deg]"]}\n')
            file.write(f'ANGLE_2 = {final_obs["Time"]} {final_obs["DEC[deg]"]}\n')
            file.write(f'\n')



# if __name__ == '__main__':
#     # main('')
#     images_directory = 'Rastreio\\images'
#     output_path = 'Resultados\\Resultados-12.txt'
#     for folder in os.listdir(images_directory):
#         for number in range(1,11):
#             filename = f'{images_directory}\\{folder}\\Capture_{number}.fits'
#             print(f'Analisando imagem {filename}.')
#             if os.path.isfile(filename):
#                 try:
#                     get_radec_from_fits(filename,output_path)
#                     print(f'Imagem {filename} analisada com sucesso.')
#                 except:
#                     print(f'Erro ao processar a imagem: {filename}')
#                     continue
#     print('Análise de imagens finalizadas.')