import os
import cv2
import sys
import numpy as np
from tqdm import tqdm
import tifffile as tiff
import multiprocessing as mp
import matplotlib.pyplot as plt
from scipy.ndimage import fourier_shift


# if not os.path.exists('/data/Haochun/ISH/Data/Sample_D4/concat_3_rounds'):
#     os.makedirs('/data/Haochun/ISH/Data/Sample_D4/concat_3_rounds')

def _get_fft_2d(img1, img2,root_folder):
    if img2.shape != img1.shape:
        img2 = cv2.resize(img2, (img1.shape[1], img1.shape[0]))
    
    plt.imshow(img2, cmap = 'gray', vmax = 75)
    plt.axis('off')
    plt.savefig(f'{root_folder}/resized_image.png', dpi = 300, bbox_inches = 'tight')
    print(f'shape of image 1 and image 2: {img1.shape}, {img2.shape}')

    # step 1
    fft_img1 = np.fft.fftn(img1)
    fft_img2 = np.fft.fftn(img2)

    # step 2
    cross_correlation = np.fft.ifftn(fft_img1 * np.conj(fft_img2))

    # step 3
    shifts = np.unravel_index(np.argmax(np.abs(cross_correlation)), cross_correlation.shape)
    print(f'shifts: {shifts}')
    translation = []
    for shift, dim in zip(shifts, img1.shape):
        if shift > dim // 2:
            shift -= dim
        translation.append(shift)
    print(f'translation: {translation}')

    return translation

def _apply_shift(img2, translation,root_folder):
    # step 4
    aligned_img2 = fourier_shift(np.fft.fftn(img2), translation)
    aligned_img2 = np.fft.ifftn(aligned_img2).real.astype(img2.dtype)

    shift_x, shift_y = translation
    shift_x = int(shift_x); shift_y = int(shift_y)
    print(shift_x, shift_y)
    if shift_x > 0:
        aligned_img2[:shift_x, :] = 0
    elif shift_x < 0:
        aligned_img2[shift_x:, :] = 0
    if shift_y > 0:
        aligned_img2[:, :shift_y] = 0
    elif shift_y < 0:
        aligned_img2[:, shift_y:] = 0
    
    plt.imshow(aligned_img2, cmap = 'gray', vmax = 75)
    plt.axis('off')
    plt.savefig(f'{root_folder}/registered_image.png', dpi = 300, bbox_inches = 'tight')

    return aligned_img2

if __name__ == '__main__':
    
    samples = ['E2','F','G','H','I','J','K','L','M','N','O','P']
    #batches = ['G1', 'G2']
    rs = ['R1','R2','R3']
    channels = ['ch01', 'ch02', 'ch03', 'ch04']
    for sample in samples:
            root_folder=f'/data2/haochun/Sennet/ISH/{sample}'
            if not os.path.exists(f'/data2/haochun/Sennet/ISH/{sample}/concat_3_rounds'):
                os.makedirs(f'/data2/haochun/Sennet/ISH/{sample}/concat_3_rounds')
            print(f'Processing {sample}')
        
            dapi_r1 = tiff.imread(f'{root_folder}/dapi_images/R1_dapi_all_v3.tif')
            dapi_r2 = tiff.imread(f'{root_folder}/dapi_images/R2_dapi_all_v3.tif')
            dapi_r3 = tiff.imread(f'{root_folder}/dapi_images/R3_dapi_all_v3.tif')
        # Calculate the shift between round 1 and round 2
            shift_1 = _get_fft_2d(dapi_r1, dapi_r2,root_folder)
            shift_2= _get_fft_2d(dapi_r1,dapi_r3,root_folder)
        # Apply the shift to the dapi image
            dapi_r2_aligned = _apply_shift(dapi_r2, shift_1,root_folder)
            tiff.imsave(f'{root_folder}/dapi_images/R2_dapi_all_v3_aligned.tif', dapi_r2_aligned)
            dapi_r3_aligned = _apply_shift(dapi_r3, shift_2,root_folder)
            tiff.imsave(f'{root_folder}/dapi_images/R3_dapi_all_v3_aligned.tif', dapi_r3_aligned)
        # Apply the same shift to all channels
            for channel in tqdm(channels):
                img_r2 = tiff.imread(f'{root_folder}/spots_images/R2_{channel}_all.tif')
                img_aligned = _apply_shift(img_r2, shift_1,root_folder)
                tiff.imsave(f'{root_folder}/spots_images/R2_{channel}_all_aligned.tif', img_aligned)

            for channel in tqdm(channels):
                img_r3 = tiff.imread(f'{root_folder}/spots_images/R3_{channel}_all.tif')
                img_aligned = _apply_shift(img_r3, shift_2,root_folder)
                tiff.imsave(f'{root_folder}/spots_images/R3_{channel}_all_aligned.tif', img_aligned)
