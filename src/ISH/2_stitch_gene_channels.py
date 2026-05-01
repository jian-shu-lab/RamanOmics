import os
import numpy as np
from tqdm import tqdm
import tifffile as tiff
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import glob

# if not os.path.exists('/data/Haochun/md3/Haochun/ISh/Ke/D0/M_3/spots_images'):
#     os.mkdir('/data/Haochun/md3/Haochun/ISh/Ke/D0/M_3/spots_images')

def _get_tile_row_col(root, num_tiles):
    tiles = np.arange(num_tiles)
    x = []; y = []
    for tile in tiles:
        x.append(float(root.find(f".//Tile[@FieldX='{tile}']").get('PosX')))
        y.append(float(root.find(f".//Tile[@FieldX='{tile}']").get('PosY')))

    xpos = np.sort(np.unique(x))
    xpos_dict = dict(zip(xpos, np.arange(len(xpos))))

    ypos = np.sort(np.unique(y))
    ypos_dict = dict(zip(ypos, np.arange(len(ypos))))

    row_dict = {}
    col_dict = {}
    for tile in tiles:
        row_dict[tile+1] = ypos_dict[y[tile]]
        col_dict[tile+1] = xpos_dict[x[tile]]
    return row_dict, col_dict

def _stitch_spots_images(sample, r, ch, vmax = 30, overlap_pixels = 80):
    root_folder = f'/data2/haochun/Sennet/ISH/{sample}'

    # if r == 'R1':
    #     root_folder = root_folder
    # elif r == 'R2':
    #     root_folder = root_folder
    # root = ET.parse(f'{root_folder}ORS1234/MetaData/20250509 NEW woundhelaing_24mo_wound site_D0_ORS1234_20250509 NEW woundhelaing_24mo_wound site_D0_ORS1234_Properties.xml').getroot()
    # root = ET.parse('/data/Haochun/md3/Haochun/ISh/Ke/D0/24mo_D0.xml').getroot()
    xml_path=glob.glob(f'{root_folder}/{r}/*.xml')
    print(xml_path)
    root = ET.parse(xml_path[0]).getroot()

    num_tiles = len(root.findall(".//Tile"))

    spots_source_file = glob.glob(f'{root_folder}/{r}/*_{ch}.tif')
    img_whole = tiff.imread(spots_source_file[0])

    row_dict, col_dict = _get_tile_row_col(root, num_tiles=num_tiles)

    ### parameters
    tiles = np.arange(num_tiles)
    num_z = int(img_whole.shape[0] / num_tiles)
    num_rows = max(row_dict.values()) + 1
    num_cols = max(col_dict.values()) + 1
    size_x = size_y = 2048
    overlap_x = overlap_y = int(size_x*0.05)
    updated_tile_size = size_x - overlap_x * 2

    ### get stitched image for spots_images
    img_merged = np.zeros(
        (num_rows * updated_tile_size, num_cols * updated_tile_size),
        dtype = np.uint8
    )

    # append data tile by tile
    for tile in tqdm(tiles):
        img_tile = img_whole[num_z * tile : num_z * (tile + 1)].max(axis = 0)    
        row, col = row_dict[tile+1], col_dict[tile+1]

        img_tile = img_tile[overlap_y:(size_y - overlap_y), :]
        img_tile = img_tile[:, overlap_x:(size_x - overlap_x)]
        img_merged[(row * updated_tile_size):((row+1) * updated_tile_size), (col * updated_tile_size) : ((col+1) * updated_tile_size)] = img_tile

    #img_merged = img_merged[::-1, :]
    # tiff.imwrite(f'{root_folder}/Sample_{sample}/spots_images/{sample}_{r}_{ch}_all.tif', img_merged)
    print(img_merged.shape)
    tiff.imwrite(f'{root_folder}/spots_images/{r}_{ch}_all.tif', img_merged)

    plt.figure(figsize = (2*num_cols, 2*num_rows))
    plt.imshow(img_merged, cmap = 'gray', vmax = vmax)
    plt.axis('off')
    plt.savefig(f'{root_folder}/spots_images/{r}_{ch}_all.png', dpi = 300, bbox_inches = 'tight')

if __name__ == '__main__':
    samples = ['E2','F','G','H','I','J','K','L','M','N','O','P']
    rs = ['R1','R2','R3']
    channels = ['ch01', 'ch02', 'ch03', 'ch04']

    for sample in samples:
        if not os.path.exists(f'/data2/haochun/Sennet/ISH/{sample}/spots_images'):
            os.mkdir(f'/data2/haochun/Sennet/ISH/{sample}/spots_images')
        for r in rs:
            for ch in channels:
                print(f'{sample}_{r}_{ch}')
                _stitch_spots_images(sample, r, ch, vmax = 30, overlap_pixels = 80)
