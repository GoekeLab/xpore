import numpy as np
import os


def makedirs(main_dir, sub_dirs=None, opt='depth'):
    if not os.path.exists(main_dir):
        os.makedirs(main_dir)
    filepaths = dict()
    if sub_dirs is not None:
        if opt == 'depth':
            path = main_dir
            for sub_dir in sub_dirs:
                path = os.path.join(path, sub_dir)
                filepaths[sub_dir] = path
                # if not os.path.exists(path):
                try:  # Use try-catch for the case of multiprocessing.
                    os.makedirs(path)
                except:
                    pass

        else:  # opt == 'breadth'
            for sub_dir in sub_dirs:
                path = os.path.join(main_dir, sub_dir)
                filepaths[sub_dir] = path
                # if not os.path.exists(path):
                try:  # Use try-catch for the case of multiprocessing.
                    os.makedirs(path)
                except:
                    pass

    return filepaths


def str_decode(df):
    str_df = df.select_dtypes([np.object])
    str_df = str_df.stack().str.decode('utf-8').unstack()
    for col in str_df:
        df[col] = str_df[col]
    return df


def str_encode(df):
    str_df = df.select_dtypes([np.object])
    str_df = str_df.stack().str.encode('utf-8').unstack()
    for col in str_df:
        df[col] = str_df[col]
    return df
