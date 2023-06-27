#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 11:30:07 2022

@author: varun
"""
import pandas as pd
import numpy as np


def read_smet_file(fname):
    """


    Parameters
    ----------
    fname : string
        It is the input file to be read and stored in a pandas dataframe

    Returns
    -------
    df : pandas dataframe

    """
    a = open(fname, "r")
    Lines = a.readlines()
    count = 0
    for line in Lines:
        count += 1
        if "fields" in line:
            tmp = line.strip()
            col_names = tmp.split(sep="=")[-1].split()
        if "altitude" in line:
            tmp = line.strip()
            altitude = tmp.split(sep="=")[-1].split()
        if "[DATA]" in line:
            break
    a.close()

    df = pd.read_csv(
        fname,
        sep="\s+",
        skiprows=count,
        parse_dates=[0],
        header=None,
        names=col_names,
        low_memory=False,
    )
    altitude = [float(x) for x in altitude][0]
    df = df.set_index(["timestamp"])
    # df.name = "df_" + fname.split("/")[-1].split(".")[0]

    if "TA" in df:
        if df["TA"].max() < 200:
            df["TA"] = df["TA"] + 273.15

    if "HS_mod" in df:
        if df["HS_mod"].max() > 100:
            df["HS_mod"] = df["HS_mod"] / 100.0

    var = "TSS"
    var_cols = [col for col in df.columns if var in col]
    for loc_var in var_cols:
        if df[loc_var].max() < 200:
            df[loc_var] = df[loc_var] + 273.15

    df = df.loc[:, ~df.eq(-999).all()]

    timeresolution = (df.index[1] - df.index[0]).seconds

    return df, altitude, timeresolution


def read_toolchain_file(fname):
    """


    Parameters
    ----------
    fname : string
        It is the input file to be read and stored in a pandas dataframe

    Returns
    -------
    df : pandas dataframe

    """
    a = open(fname, "r")
    Lines = a.readlines()
    count = 0
    for line in Lines:
        count += 1
        if "fields" in line:
            tmp = line.strip()
            col_names = tmp.split(sep="=")[-1].split()
        if "altitude" in line:
            tmp = line.strip()
            altitude = tmp.split(sep="=")[-1].split()
        if "[DATA]" in line:
            break
    a.close()
    col_names = [x.replace(",", "") for x in col_names]

    df = pd.read_csv(
        fname,
        sep="\s+",
        skiprows=count,
        header=None,
        names=col_names,
        low_memory=False,
    )

    time_delta = "30s"

    times = pd.date_range(start="2020-08-01 01:00:00", periods=len(df), freq=time_delta)
    df["timestamp"] = times
    df = df.set_index(["timestamp"])

    # df.name = "df_" + fname.split("/")[-1].split(".")[0]
    print(df.columns)
    df["TA"] = df["t_a"]
    df["RH"] = df["rh"]
    df["VW"] = df["uv"]
    df["ISWR"] = df["sw_d"]
    df["ILWR"] = df["lw_d"]
    df["PSUM"] = df["p_sum"]

    df["TSG"] = 273.15

    if "TA" in df:
        if df["TA"].max() < 200:
            df["TA"] = df["TA"] + 273.15

    if "HS_mod" in df:
        if df["HS_mod"].max() > 100:
            df["HS_mod"] = df["HS_mod"] / 100.0

    var = "TSS"
    var_cols = [col for col in df.columns if var in col]
    for loc_var in var_cols:
        if df[loc_var].max() < 200:
            df[loc_var] = df[loc_var] + 273.15

    df = df.loc[:, ~df.eq(-999).all()]
    df.name = "df_" + fname.split("/")[-1].split(".")[0]
    return df, float(altitude[0])


def read_cosmo_output(fname):

    cosmo_model = fname.split("/")[-1].split("_")[1]
    print(str(cosmo_model))
    if "with" in fname.split("/")[-1]:
        bool_snopo = True
        print("bool_snopo: " + str(bool_snopo))
    else:
        bool_snopo = False
        print("bool_snopo: " + str(bool_snopo))

    if bool_snopo:
        col_names = [
            "u",
            "v",
            "TA",
            "qv",
            "ps",
            "prr_gsp",
            "prs_gsp",
            "prg_gsp",
            "ISWR_dir",
            "ISWR_diff",
            "OSWR",
            "ILWR",
            "OLWR",
            "TSS_mod",
            "HS_mod",
            "T_bottom",
            "top",
            "hn",
            "tch",
        ]
    else:
        col_names = [
            "u",
            "v",
            "TA",
            "qv",
            "ps",
            "prr_gsp",
            "prs_gsp",
            "prg_gsp",
            "SWNET",
            "LWR_net",
            "HS_mod",
            "TSS_mod",
            "T_bottom",
            "Qs",
            "Ql",
        ]

    if int(cosmo_model) == 1:
        time_delta = "10s"

    if int(cosmo_model) == 2:
        time_delta = "20s"

    print(time_delta)
    df = pd.read_csv(fname, sep=",", header=None, names=col_names, low_memory=False)

    times = pd.date_range(start="2020-09-24 12:00:00", periods=len(df), freq=time_delta)
    df["timestamp"] = times
    df = df.set_index(["timestamp"])

    # set wind speed
    df["VW"] = np.sqrt(df["u"] ** 2.0 + df["v"] ** 2.0)
    df.drop(columns=["u", "v"], inplace=True)

    # set temperature
    var = "TSS"
    var_cols = [col for col in df.columns if var in col]
    for loc_var in var_cols:
        df[loc_var] = df[loc_var] - 273.15
    df["TA"] = df["TA"] - 273.15

    # set snow height
    df["HS_mod"] = df["HS_mod"] * 100.0

    if bool_snopo:
        df["SWNET"] = df["ISWR_diff"] + df["ISWR_dir"] - df["OSWR"]
        df["LWR_net"] = df["ILWR"] - df["OLWR"]

    df.name = "df_" + fname.split("/")[-1].split(".")[0]
    return df


def read_ssa_debug(fname, dt):

    col_names = header_str = [
        "for_sn",
        "ISWR_dir",
        "ISWR_diff",
        "OSWR",
        "SWNET",
        "SWNET_SNOW",
        "SW_BOTTOM",
        "top",
        "ILWR",
        "OLWR",
        "Ql",
        "Qs",
        "TSS_mod",
        "HS_MOD",
        "new snow",
        "T_bottom",
        "dz_bottom",
        "dz_top",
        "theta_w",
        "theta_i",
        "rho_sn_top",
        "rho_sn_bottom",
        "merge_p",
        "TA",
        "QV",
    ]
    df = pd.read_csv(fname, sep=",", header=None, names=col_names, low_memory=False)

    time_delta = str(dt) + "s"
    times = pd.date_range(start="2020-09-24 12:00:00", periods=len(df), freq=time_delta)
    df["timestamp"] = times
    df = df.set_index(["timestamp"])
    df.name = "df_" + fname.split("/")[-1].split(".")[0]
    return df


def read_ssa_output(fname, dt):
    col_names = ["n", "top", "HS_mod", "TSS_mod", "SW_ABS", "ILWR", "OLWR", "Ql", "Qs"]
    df = pd.read_csv(fname, sep="\s+", header=None, names=col_names, low_memory=False)
    time_delta = str(dt) + "s"
    times = pd.date_range(start="2020-08-01 01:00:00", periods=len(df), freq=time_delta)
    df["timestamp"] = times
    df = df.set_index(["timestamp"])
    df.name = "df_" + fname.split("/")[-1].split(".")[0]
    return df
