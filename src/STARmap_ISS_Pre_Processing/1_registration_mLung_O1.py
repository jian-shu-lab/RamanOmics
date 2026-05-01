import os
import sys
import shutil
import multiprocessing as mp
import numpy as np
import tifffile as tiff
import matplotlib.pyplot as plt
import scipy.fftpack as fft
import scipy.ndimage as ndi
from scipy.ndimage import fourier_shift
from scipy.ndimage import fourier_gaussian
from skimage.feature import register_translation

channels = [1, 2, 3, 4]
channel_dict = {1: 1, 2: 2, 3: 3, 4: 4}

def apply_fft_3d(img1, img2):
    img1_query = np.mean(img1, axis = -1).astype(img1.dtype)
    img2_query = np.mean(img2, axis = -1).astype(img2.dtype)

    fft_img1 = np.fft.fftn(img1_query)
    fft_img2 = np.fft.fftn(img2_query)

    cross_correlation = np.fft.ifftn(fft_img1 * np.conj(fft_img2))

    shifts = np.unravel_index(np.argmax(np.abs(cross_correlation)), cross_correlation.shape)
    print(shifts, type(shifts))
    translation = []
    for shift, dim in zip(shifts, img1_query.shape):
        if shift > dim // 2:
            shift -= dim
        translation.append(shift)

    for channel in np.arange(img2.shape[-1]):
        aligned_img2 = fourier_shift(np.fft.fftn(img2[:,:,:,channel]), translation)
        aligned_img2 = np.fft.ifftn(aligned_img2).real.astype(img2.dtype)
        if channel == 0:
            aligned_img2_whole = aligned_img2[:,:,:,None]
        else:
            aligned_img2_whole = np.concatenate((aligned_img2_whole, aligned_img2[:,:,:,None]), axis = 3)

    shift_z, shift_x, shift_y = translation
    shift_x = int(shift_x); shift_y = int(shift_y); shift_z = int(shift_z)
    print(shift_x, shift_y, shift_z)
    if shift_x > 0:
        aligned_img2_whole[:, :shift_x, :, :] = 0
    elif shift_x < 0:
        aligned_img2_whole[:, shift_x:, :, :] = 0
    if shift_y > 0:
        aligned_img2_whole[:, :, :shift_y, :] = 0
    elif shift_y < 0:
        aligned_img2_whole[:, :, shift_y:, :] = 0
    if shift_z > 0:
        aligned_img2_whole[:shift_z, :, :, :] = 0
    elif shift_z < 0:
        aligned_img2_whole[shift_z:, :, :, :] = 0
    return translation, aligned_img2_whole

def concatenate_images(tile):
    root_folder = f'/data/framont/mLung_O1/Deconvoluted/'
    rounds = [1,2,3,4,5,6]

    # Initialize img_round before the loop
    img_round = None
    
    for round_id in rounds:
        print(f'tile {tile}: concatenating round {round_id}')
        for ch in channels:
            img = tiff.imread(f'{root_folder}mLung_O1_tile_{tile}_R{round_id}_ch0{ch}_cmle.tif')
            img = img[:,:,:,None]

            # concatenate channels
            if ch == 1:
                img_round = img
            else:
                img_round = np.concatenate((img_round, img), axis=3)
        
        # concatenate rounds
        if round_id == 1:
            img_tile = img_round[None, :, :, :, :]
        else:
            img_tile = np.concatenate((img_tile, img_round[None, :, :, :, :]), axis=0)
        
    return img_tile

def get_local_shifts(img_tile, tile, ref_round = 1):
    # ref_round: 1
    query_rounds = [2,3,4,5,6]
    img_ref = img_tile[ref_round-1, :, :, :, :] 

    for query_round in query_rounds:
        img_query = img_tile[query_round-1, :, :, :, :]
        shift, img_aligned = apply_fft_3d(img_ref, img_query)
        print(f'tile {tile}: round {query_round} shift: {shift}')

        for ch in channels:
            img_aligned_ch = img_aligned[:,:,:,channel_dict[ch]-1]
            tiff.imwrite(f'/data/framont/mLung_O1/tile_data_registered/tile_{tile}/mLung_O1_R{query_round}_ch0{ch}.tif', img_aligned_ch)
        
def copy_files(tile):
    root_folder = f'/data/framont/mLung_O1/Deconvoluted/'
    gene_tsv_file = '/data/framont/mLung_O1/Deconvoluted/genes.csv'

    # copy round 1
    for ch in channels:
        shutil.copy(f'{root_folder}mLung_O1_tile_{tile}_R1_ch0{ch}_cmle.tif', f'/data/framont/mLung_O1/tile_data_registered/tile_{tile}/mLung_O1_R1_ch0{ch}.tif')
    
    # copy genes
    shutil.copy(gene_tsv_file, f'/data/framont/mLung_O1/tile_data_registered/tile_{tile}/genes.csv')

def run_register_tile(tile):
    print(f'tile {tile}')
    if not os.path.exists('/data/framont/mLung_O1/tile_data_registered'):
        os.makedirs('/data/framont/mLung_O1/tile_data_registered', exist_ok=True)
    if not os.path.exists(f'/data/framont/mLung_O1/tile_data_registered/tile_{tile}'):
        os.makedirs(f'/data/framont/mLung_O1/tile_data_registered/tile_{tile}', exist_ok=True)
    if not os.path.exists('logging'):
        os.makedirs('logging', exist_ok=True)

    img_tile = concatenate_images(tile)
    get_local_shifts(img_tile, tile, ref_round = 1)
    copy_files(tile)


if __name__ == '__main__':
    print('Start')
    tiles = np.arange(16)

    # avail_tiles = [i for i in tiles if os.path.exists(f'tile_data_registered/tile_{i}')]
    # print(len(avail_tiles))
    tiles = [i for i in tiles if not os.path.exists(f'/data/framont/mLung_O1/tile_data_registered/tile_{i}')]
    # for tile in avail_tiles:
    #     print(tile, len(os.listdir(f'tile_data_registered/tile_{tile}')))
    #     # run_register_tile(tile)
    
    n_process = 8
    with mp.Pool(processes=n_process) as pool:
        items = list(tiles)
        results = pool.map(run_register_tile, items)

    print('Done')
