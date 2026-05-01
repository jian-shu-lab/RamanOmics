import os
import numpy as np
from tqdm import tqdm
import tifffile as tiff
import matplotlib as mpl
import matplotlib.pyplot as plt

# tile -> row
row_dict = {
    16: 0, 15: 0, 14: 0, 13: 0,
    9: 1, 10: 1, 11: 1, 12: 1,
    8: 2, 7: 2, 6: 2, 5: 2,
    1: 3, 2: 3, 3: 3, 4: 3,
}

col_dict = {
    1: 0, 8: 0, 9: 0, 16: 0,
    15: 1, 10: 1, 7: 1, 2: 1,
    14: 2, 11: 2, 6: 2, 3: 2,
    13: 3, 12: 3, 5: 3, 4: 3,
}

dapi_source_file = '/data/framont/mLung_O1/DAPI/mLung_O1_R1_RAW_ch00.tif'
dapi = tiff.imread(dapi_source_file)

### parameters
tiles = np.arange(16)
num_z =  49
num_rows = 4
num_cols = 4
size_x = size_y = 2048
overlap_x = overlap_y = 110
updated_tile_size = size_x - overlap_x * 2

### read data and get individual dapi
img_whole = tiff.imread(dapi_source_file)
if not os.path.exists('dapi_images'):
    os.mkdir('dapi_images')
for tile in tiles:
    print(tile)
    img_tile = img_whole[num_z * tile : num_z * (tile + 1)].max(axis = 0)
    tiff.imwrite(f'dapi_images/dapi_tile_{tile}.tif', img_tile)
    plt.figure(figsize = (2*num_cols, 2*num_rows))
    plt.imshow(img_tile, cmap = 'gray')
    plt.axis('off')
    plt.savefig(f'dapi_images/dapi_tile_{tile}.png', dpi = 300, bbox_inches = 'tight')

### get stitched image for dapi
img_whole = tiff.imread(dapi_source_file)
img_merged = np.zeros(
    (num_rows * updated_tile_size, num_cols * updated_tile_size),
    dtype = np.uint8
)

### merge data tile by tile
for tile in tqdm(tiles):
    img_tile = img_whole[num_z * tile : num_z * (tile + 1)].max(axis = 0)    
    row, col = row_dict[tile+1], col_dict[tile+1]

    img_tile = img_tile[overlap_y:(size_y - overlap_y), :]
    img_tile = img_tile[:, overlap_x:(size_x - overlap_x)]
    img_tile = img_tile[::-1,:]
    img_merged[(row * updated_tile_size):((row+1) * updated_tile_size), (col * updated_tile_size) : ((col+1) * updated_tile_size)] = img_tile
    tiff.imwrite(f'dapi_images/dapi_tile_whole_{tile}_v3.tif', img_tile)

### save
img_merged = img_merged[::-1, :]
tiff.imwrite('dapi_all_v3.tif', img_merged)
print(img_merged.shape)
plt.figure(figsize = (2*num_cols, 2*num_rows))
plt.imshow(img_merged, cmap = 'gray')
plt.axis('off')
plt.savefig('dapi_all_v3.png', dpi = 300, bbox_inches = 'tight')
