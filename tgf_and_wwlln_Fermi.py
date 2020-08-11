#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 14:49:35 2019

@author: kate
"""

import matplotlib.pyplot as plt
from PIL import Image
import numpy as np
#from mpl_toolkits.basemap import Basemap
import pandas
import datetime  
import csv
import math


def scalar_multiplication(vector_1, vector_2):
    ans = 0
    for i in range(len(vector_1)):
        ans = ans + vector_1[i] * vector_2[i]
    return ans


# __________returns angle between 2 points in spherical coordinates
def angle_two_points(lon1, lat1, r1, lon2, lat2, r2):
    lat1_rad = lat1 * math.pi / 180
    lon1_rad = lon1 * math.pi / 180
    lat2_rad = lat2 * math.pi / 180
    lon2_rad = lon2 * math.pi / 180
    descartes_1 = np.array([r1 * math.cos(lat1_rad) * math.cos(lon1_rad), r1 * math.cos(lat1_rad) * math.sin(lon1_rad), r1 * math.sin(lat1_rad)])
    descartes_2 = np.array([r2 * math.cos(lat2_rad) * math.cos(lon2_rad), r2 * math.cos(lat2_rad) * math.sin(lon2_rad), r2 * math.sin(lat2_rad)])
    angle = 180 / math.pi * math.acos(scalar_multiplication(descartes_1, descartes_2) / np.sqrt(scalar_multiplication(descartes_1, descartes_1) * scalar_multiplication(descartes_2, descartes_2)))
    return angle


# _________a function to find the angle between the TGF source (correlated lightning source) and the detector
def detector_cloud_angle(lon1, lat1, r1, lon2, lat2, r2):
    # conversion from a spherical coordinate system to a Cartesian one
    lat1_rad = lat1 * math.pi / 180
    lon1_rad = lon1 * math.pi / 180
    lat2_rad = lat2 * math.pi / 180
    lon2_rad = lon2 * math.pi / 180
    descartes_1 = np.array([r1 * math.cos(lat1_rad) * math.cos(lon1_rad), r1 * math.cos(lat1_rad) * math.sin(lon1_rad),
                            r1 * math.sin(lat1_rad)])
    descartes_2 = np.array([r2 * math.cos(lat2_rad) * math.cos(lon2_rad), r2 * math.cos(lat2_rad) * math.sin(lon2_rad),
                            r2 * math.sin(lat2_rad)])
    # calculation of the vector between the cloud and the detector
    vector_between_points = descartes_1 - descartes_2
    # calculation of the vector normal to the cloud surface in the Cartesian system. In spherical coordinates it has the
    # same lon and alt, but r = r2 + some small number
    normal_to_vector_2 = np.array(
        [(r2 + 1) * math.cos(lat2_rad) * math.cos(lon2_rad), (r2 + 1) * math.cos(lat2_rad) * math.sin(lon2_rad),
         (r2 + 1) * math.sin(lat2_rad)])
    # classical way to find an angle between two vectors in the Cartesian system
    angle = 180 / math.pi * math.acos(scalar_multiplication(vector_between_points, normal_to_vector_2) / np.sqrt(
        scalar_multiplication(vector_between_points, vector_between_points) * scalar_multiplication(normal_to_vector_2,
                                                                                                    normal_to_vector_2)))
    return angle


def distance_two_points(lon1, lat1, r1, lon2, lat2, r2):
    lat1_rad = lat1 * math.pi / 180
    lon1_rad = lon1 * math.pi / 180
    lat2_rad = lat2 * math.pi / 180
    lon2_rad = lon2 * math.pi / 180
    delta_x = r1 * math.cos(lat1_rad) * math.cos(lon1_rad) - r2 * math.cos(lat2_rad) * math.cos(lon2_rad)
    delta_y = r1 * math.cos(lat1_rad) * math.sin(lon1_rad) - r2 * math.cos(lat2_rad) * math.sin(lon2_rad)
    delta_z = r1 * math.sin(lat1_rad) - r2 * math.sin(lat2_rad)
    distance = math.sqrt(pow(delta_x, 2) + pow(delta_y, 2) + pow(delta_z, 2))
    return distance  # result in kilometers


def coord_Dekart_from_speric(lon1, lat1, r1):
    lat1_rad = lat1 * math.pi / 180
    lon1_rad = lon1 * math.pi / 180
    x = r1 * math.cos(lat1_rad) * math.cos(lon1_rad)
    y = r1 * math.cos(lat1_rad) * math.sin(lon1_rad)
    z = r1 * math.sin(lat1_rad)
    r = math.sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2))
    horizontal = math.sqrt(pow(x, 2) + pow(y, 2))
    return x, y, z, r, horizontal   # result in kilometers


def create_datetime_array_from_file(data):
    datetime_from_data = []
    for jj in range(0, len(data)):
        jj_year = int(data.iloc[jj, 0][0:4])
        jj_month = int(data.iloc[jj, 0][5:7])
        jj_day = int(data.iloc[jj, 0][8:10])
        jj_hour = int(data.iloc[jj, 1][0:2])
        jj_minute = int(data.iloc[jj, 1][3:5])
        jj_second = int(data.iloc[jj, 1][6:8])
        jj_mcs = int(data.iloc[jj, 1][9:15])
        datetime_from_data.append(datetime.datetime(jj_year, jj_month, jj_day, jj_hour, jj_minute, jj_second, jj_mcs))
    return datetime_from_data


def create_datetime_array_for_Fermi(data):
    datetime_from_data = []
    for jj in range(0, len(data)):
        jj_year = int(data.iloc[jj, 6][1:5])
        jj_month = int(data.iloc[jj, 6][6:8])
        jj_day = int(data.iloc[jj, 6][9:11])
        jj_hour = int(data.iloc[jj, 7][1:3])
        jj_minute = int(data.iloc[jj, 7][4:6])
        jj_second = int(data.iloc[jj, 7][7:9])
        jj_mcs = int(data.iloc[jj, 7][10:16])
        datetime_from_data.append(datetime.datetime(jj_year, jj_month, jj_day, jj_hour, jj_minute, jj_second, jj_mcs))
    return datetime_from_data


def plot_tgf_detection_angle_distribution():
    data = pandas.read_csv('gbm_tgf_catalog_wwlln.csv', index_col=False, header=None)

    r_detector = 6378.1 + 565.0  # km
    r_cloud = 6378.1 + 10.0  # km

    tgf_angles = []
    tgf_distances = []

    for j_tgf in range(1, len(data)):
        tgf_angles.append(detector_cloud_angle(float(data.iloc[j_tgf, 4]), float(data.iloc[j_tgf, 5]), r_detector, float(data.iloc[j_tgf, 7]), float(data.iloc[j_tgf, 8]), r_cloud))
        tgf_distances.append(distance_two_points(float(data.iloc[j_tgf, 4]), float(data.iloc[j_tgf, 5]), r_detector, float(data.iloc[j_tgf, 7]), float(data.iloc[j_tgf, 8]), r_cloud))

    plt.hist(tgf_angles)
    plt.xlabel('TGF observation angle, degrees')
    plt.ylabel('Distribution')
    plt.show()


def main():

    # __________
    # lon1 = 112.1
    # lat1 = 0.49
    # r1 = 6400.0 + 500.0
    #
    # lon2 = 11.642
    # lat2 = 8.2996
    # r2 = 6400.0 + 500.0
    #
    # lon = 111.6418
    # lat = -8.2996
    # r = 6400.0 + 10.0
    #
    # print(detector_cloud_angle(lon1, lat1, r1, lon, lat, r))
    # print(detector_cloud_angle(lon2, lat2, r2, lon, lat, r))
    # __________

    # __________a script to plot tgf angle detection distribution considering detection probability properties
    data = pandas.read_csv('gbm_tgf_catalog_wwlln.csv', index_col=False, header=None)

    r_detector = 6378.1 + 565.0  # km
    r_cloud = 6378.1 + 10.0  # km

    tgf_angles = []
    tgf_distances = []

    for j_tgf in range(1, len(data)):
        tgf_angles.append(detector_cloud_angle(float(data.iloc[j_tgf, 4]), float(data.iloc[j_tgf, 5]), r_detector, float(data.iloc[j_tgf, 7]), float(data.iloc[j_tgf, 8]), r_cloud))
        tgf_distances.append(distance_two_points(float(data.iloc[j_tgf, 4]), float(data.iloc[j_tgf, 5]), r_detector, float(data.iloc[j_tgf, 7]), float(data.iloc[j_tgf, 8]), r_cloud))

    # plt.hist(tgf_angles)
    # plt.xlabel('TGF observation angle, degrees')
    # plt.ylabel('Distribution')
    # plt.show()

    angle_histogram = np.histogram(tgf_angles)
    print(angle_histogram)
    distribution = np.zeros(len(angle_histogram[1]) - 1)
    normalization_factor = 0
    mean_distances = np.zeros(len(angle_histogram[1]) - 1)

    for i in range(len(angle_histogram[1]) - 1):
        for j in range(len(tgf_angles)):
            if (tgf_angles[j] >= angle_histogram[1][i]) and (tgf_angles[j] < angle_histogram[1][i + 1]):
                mean_distances[i] = mean_distances[i] + tgf_distances[j]
        mean_distances[i] = mean_distances[i] / angle_histogram[0][i]
        # distribution[i] = angle_histogram[0][i] / np.sin(math.pi / 180 * angle_histogram[1][i + 1])
        distribution[i] = angle_histogram[0][i] / (np.cos(math.pi / 180 * angle_histogram[1][i]) - np.cos(math.pi / 180 * angle_histogram[1][i + 1])) * (math.pi / 180 * angle_histogram[1][i + 1] - math.pi / 180 * angle_histogram[1][i])
        normalization_factor = normalization_factor + distribution[i]
    plt.plot(angle_histogram[1][0: len(angle_histogram[1]) - 1], distribution / normalization_factor)
    plt.xlabel('TGF observation angle, degrees')
    plt.ylabel('Distribution')
    plt.show()
    # __________

    # time_difference = datetime.timedelta(microseconds=500)
    #
    # # ==============================================================================
    # # data = pandas.read_csv('gbm_tgf_catalog_2016_f.csv', index_col=False, header=None)
    # # f_TGF_info = open('Fermi_TGFs_distances_and_time_diff' + '.csv','w')
    # # f_TGF_with_lightning = open('Fermi_TGFs_2016_with_lightning_' + str(time_difference.microseconds) + 'mcs.csv','w')
    # # f_TGF_no_lightning = open('Fermi_TGFs_2016_no_lightning_'+ str(time_difference.microseconds) + 'mcs.csv','w')
    # # ==============================================================================
    #
    # data = pandas.read_csv('gbm_tgf_catalog_2016.csv', index_col=False, header=None)
    # f_TGF_info = open('Fermi_TGFs_distances_and_time_diff' + '.csv', 'w')
    # f_TGF_with_lightning = open('Fermi_TGFs_2016_with_lightning_' + str(time_difference.microseconds) + 'mcs.csv', 'w')
    # f_TGF_no_lightning = open('Fermi_TGFs_2016_no_lightning_' + str(time_difference.microseconds) + 'mcs.csv', 'w')
    # # Files where will be TGFs altitude, number of counts, duration
    #
    # lon_tgf = []
    # lat_tgf = []
    # alt_tgf = []
    # counts_NaI = []
    # duration_tgf = []
    # Earth_radius = 6371.0
    #
    # # r_l = Earth_radius + 10.0
    # r_l = Earth_radius
    # speed_of_light = 299792.458  # km per s
    #
    # # let's initialize arrays with lat, lon and datetime values for all the TGFs
    # for j_tgf in range(0, len(data)):
    #     lon_tgf.append((data.iloc[j_tgf, 10]))
    #     lat_tgf.append((data.iloc[j_tgf, 11]))
    #     alt_tgf.append((data.iloc[j_tgf, 12]))
    #     counts_NaI.append((data.iloc[j_tgf, 5]))
    #     duration_tgf.append((data.iloc[j_tgf, 8]))
    #
    # time_tgf = create_datetime_array_for_Fermi(data)
    # quantity_of_columns = 1
    #
    # for j_tgf in range(0, len(data)):
    #     f_TGF_info.write(str(time_tgf[j_tgf]) + '   ')
    #
    #     date_and_time = data.iloc[j_tgf, 6][1:5] + data.iloc[j_tgf, 6][6:8] + data.iloc[j_tgf, 6][9:11] + '_'
    #     date_and_time += data.iloc[j_tgf, 7][1:3] + data.iloc[j_tgf, 7][4:6] + data.iloc[j_tgf, 7][7:9]
    #
    # #     data_wwlln = pandas.read_csv('WWLLN_csv/F' + date_and_time + '.csv', index_col=False, header=None)
    # #     time_wwlln = create_datetime_array_from_file(data_wwlln)
    # #
    # #     f_current_lightning = open('dif/F' + date_and_time + '_time_diff' + '.csv', 'w')
    # #     minimal_time_diff = datetime.timedelta(seconds=1)
    # #     minimal_time_distance = 0
    # #     min_time_offset = 0
    # #     min_time_sign = '!'
    # #     sign = '%'
    # #     for j_wwlln in range (0, len(data_wwlln)):
    # #         lat_l = float(data_wwlln.iloc[j_wwlln, 2])
    # #         lon_l = float(data_wwlln.iloc[j_wwlln, 3])
    # # #        time_for_light_transition = datetime.timedelta(0,0,0)
    # #         distance_for_light_transition = distance_two_points(lon_tgf[j_tgf], lat_tgf[j_tgf], alt_tgf[j_tgf] + Earth_radius, lon_l, lat_l, r_l)
    # #         offset = distance_two_points(lon_tgf[j_tgf], lat_tgf[j_tgf], Earth_radius, lon_l, lat_l, Earth_radius)
    # # #        distance_for_light_transition = alt_tgf[j_tgf] - offset
    # #         time_for_light_transition = distance_for_light_transition / speed_of_light  # seconds
    # #         seconds = math.modf(time_for_light_transition)[1]
    # #         microseconds = math.modf(time_for_light_transition)[0] * pow(10, 6)
    # #         time_for_light_transition = datetime.timedelta(seconds=seconds, microseconds=microseconds)
    # #         time_diff = abs(time_tgf[j_tgf] - (time_wwlln[j_wwlln] + time_for_light_transition))
    # #         if time_tgf[j_tgf] > (time_wwlln[j_wwlln] + time_for_light_transition):
    # #             sign = '-'
    # #         else:
    # #             sign = '+'
    # # #        time_diff = abs(time_tgf[j_tgf] - (time_wwlln[j_wwlln]))
    # #         f_current_lightning.write(sign+str(time_diff) + '   ' + str(distance_for_light_transition) + '   WWLLN_lon = ' + str(lon_l) + '   WWLLN_lat = ' + str(lat_l) + '   offset = ' + str(offset) + '\n')
    # #         if time_diff < minimal_time_diff:
    # #             minimal_time_diff = time_diff
    # #             min_time_sign = sign
    # #             minimal_time_distance = distance_for_light_transition
    # #             min_time_offset = offset
    # #             date_stroke = data_wwlln.iloc[j_wwlln, 0]
    # #             time_stroke = data_wwlln.iloc[j_wwlln, 1]
    # #             lat_stroke = data_wwlln.iloc[j_wwlln, 2]
    # #             lon_stroke = data_wwlln.iloc[j_wwlln, 3]
    # #             residual_fit_error = data_wwlln.iloc[j_wwlln, 4]
    # #             number_of_stations = data_wwlln.iloc[j_wwlln, 5]
    # #     f_current_lightning.close()
    # #     f_TGF_info.write('   time_diff = ' + min_time_sign+ str(minimal_time_diff) + '   offset = ' + str(min_time_offset) + '   distance = ' + str(minimal_time_distance) + '\n')
    # #     str_tgf_param = 'date_TGF = ' + str(data.iloc[j_tgf, 6]) + '   time_TGF = ' + str(data.iloc[j_tgf, 7])
    # #     str_tgf_param += '   lat_TGF = ' + str(lat_tgf[j_tgf]) + '   lon_TGF = ' + str(lon_tgf[j_tgf])
    # #     str_tgf_param += '   alt_TGF = ' + str(alt_tgf[j_tgf])+'   counts_TGF = ' + str(counts_NaI[j_tgf])
    # #     str_tgf_param += '   duration = ' + str(duration_tgf[j_tgf])
    # #     str_stroke_param = 'date_stroke = ' + str(date_stroke) + '   time_stroke = ' + str(time_stroke)
    # #     str_stroke_param += '   lat_stroke = ' + str(lat_stroke) + '   lon_stroke = ' + str(lon_stroke)
    # #     str_stroke_param += '   residual_fit_error = ' + str(residual_fit_error) + '   number_of_stations = ' + str(number_of_stations)
    # #     str_comparison_param = 'offset = ' + str(min_time_offset) + '   distance_from_TGF_registration_to_lightning = ' + str(minimal_time_distance)
    # #
    # #     print('min_time_dif = ', min_time_sign, minimal_time_diff, '\n', 'min_time_offset = ', min_time_offset)
    # #     if minimal_time_diff < time_difference:
    # #         f_TGF_with_lightning.write(str_tgf_param + '\n' + str_stroke_param + '\n' + str_comparison_param + '\n' + '\n')
    # #     else:
    # #         f_TGF_no_lightning.write(str_tgf_param + '\n' + str_stroke_param + '\n' + str_comparison_param + '\n' + '\n')
    #
    # f_TGF_info.close()
    # f_TGF_with_lightning.close()
    # f_TGF_no_lightning.close()


if __name__ == '__main__':
    main()
