#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 11:29:25 2022

@author: varun
"""
import read_data
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# READ SNOWPACK REFERENCE DATA
# snowpack_folder = (
#    "/home/varun/WORKING/snowpack_vs_ssa/snowpack_from_trunk/sim_folder/output/"
# )
# fname = snowpack_folder + "WFJ_PSUM_RUN.smet"
# df_snowpack_ref = read_data.read_smet_file(fname)

# READ SNOWPACK FORCING DATA
toolchain_folder = "./input_smet/"
fname = toolchain_folder + "KLO.txt"
df_snowpack_forcing, d_alt = read_data.read_toolchain_file(fname)

df_snowpack_forcing["ISWR_dir"] = df_snowpack_forcing["ISWR"]
df_snowpack_forcing["ISWR_dif"] = df_snowpack_forcing["ISWR"] * 0.0

# df_snowpack_forcing["TA"] = df_snowpack_forcing["TA"] + 273.15
##
const_g = 9.80665
const_lapse_rate = 0.0065
const_gaz_dry_air = 287.058
const_std_press = 101325
const_earth_r0 = 6356766.0
const_std_temp = 288.15
altitude = d_alt

expo = const_g / (const_lapse_rate * const_gaz_dry_air)
P = const_std_press * np.power(
    1.0
    - (
        (const_lapse_rate * const_earth_r0 * altitude)
        / (const_std_temp * (const_earth_r0 + altitude))
    ),
    expo,
)


def calc_qsat(T, P):

    if T < 273.16:
        c2 = 21.88
        c3 = 7.66
    else:
        c2 = 17.27
        c3 = 35.86

    expo = c2 * (T - 273.16) / (T - c3)

    vap_pressure = 273.16 * np.exp(expo)

    vp = 1 / (P - 0.378 * vap_pressure)
    qsat = 0.622 * vp * vap_pressure
    return qsat


def calc_qsat_noah(T, P):

    #
    #! For water vapor (temperature range 0C-100C)
    #
    a0 = 6.11213476
    a1 = 0.444007856
    a2 = 0.143064234e-01
    a3 = 0.264461437e-03
    a4 = 0.305903558e-05
    a5 = 0.196237241e-07
    a6 = 0.892344772e-10
    a7 = -0.373208410e-12
    a8 = 0.209339997e-15

    #
    # For ice (temperature range -75C-0C)
    #
    c0 = 6.11123516
    c1 = 0.503109514
    c2 = 0.188369801e-01
    c3 = 0.420547422e-03
    c4 = 0.614396778e-05
    c5 = 0.602780717e-07
    c6 = 0.387940929e-09
    c7 = 0.149436277e-11
    c8 = 0.262655803e-14

    T_limit = T - 273.15
    if T_limit > 100.0:
        T_limit = 100.0
    if T_limit < -75.0:
        T_limit = -75.0

    td = T_limit
    if td >= 0.0:
        es = a0 + td * (
            a1
            + td
            * (a2 + td * (a3 + td * (a4 + td * (a5 + td * (a6 + td * (a7 + td * a8))))))
        )
    else:
        es = c0 + td * (
            c1
            + td
            * (c2 + td * (c3 + td * (c4 + td * (c5 + td * (c6 + td * (c7 + td * c8))))))
        )

    es = es * 100.0  # pa

    vp = 1.0 / (P - 0.378 * es)
    vp1 = 0.622 * vp

    qs = es * vp1  # kg/kg

    return qs


df_snowpack_forcing["P"] = P
df_snowpack_forcing["QV"] = df_snowpack_forcing.apply(
    lambda row: row["RH"] * calc_qsat_noah(row["TA"], row["P"]), axis=1
)

df_snowpack_forcing = df_snowpack_forcing[
    ["TA", "P", "QV", "VW", "ISWR_dir", "ISWR_dif", "ILWR", "PSUM", "TSG"]
]

# df_snowpack_forcing["TA"] = df_snowpack_forcing["TA"] + 273.15
df_snowpack_forcing["PSUM"] = df_snowpack_forcing["PSUM"] / 30.0


df_snowpack_forcing.to_csv(
    "./ssa_KLO.txt", float_format="%g", header=False, index=False
)
