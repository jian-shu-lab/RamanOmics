import os
import numpy as np
from tqdm import tqdm
import tifffile as tiff
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import glob

# init
# if not os.path.exists('/data/Haochun/ISH/Data/Ke_Final/24mo_D0/dapi_images'):
#     os.mkdir('/data/Haochun/ISH/Data/Ke_Final/24mo_D0/dapi_images')

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

def _stitch_dapi_images(sample, r, vmax = 100, overlap_pixels = 80):
    root_folder = f'/data2/haochun/Sennet/ISH/{sample}/'
    # if r == 'R1':
    #     root_folder = root_folder
    # elif r == 'R2':
    #     root_folder = root_folder
    # elif r=='R3':
    #     root_folder=root_folder
    xml_path=glob.glob(f'{root_folder}/{r}/*.xml')
    print(xml_path)
    root = ET.parse(xml_path[0]).getroot()
    # root = ET.parse(f'{root_folder}ORS1234/MetaData/20250509 NEW woundhelaing_24mo_wound site_D0_ORS1234_20250509 NEW woundhelaing_24mo_wound site_D0_ORS1234_Properties.xml').getroot()
    # root = ET.parse('/data/Haochun/md3/Haochun/ISh/Ke/D0/24mo_D0.xml').getroot()

    num_tiles = len(root.findall(".//Tile"))

    # dapi_source_file = f'{root_folder}/Sample_{sample}/{r}/{sample} {r}_ch00.tif'
    dapi_source_file=glob.glob(f'{root_folder}/{r}/*_ch00.tif')
    # dapi_source_file=f'{root_folder}ORS1234/20250509 NEW woundhelaing_24mo_wound site_D0_ORS1234_20250509 NEW woundhelaing_24mo_wound site_D0_ORS1234_ch00.tif'
    img_whole = tiff.imread(dapi_source_file[0])

    row_dict, col_dict = _get_tile_row_col(root, num_tiles=num_tiles)

    ### parameters
    tiles = np.arange(num_tiles)
    num_z = int(img_whole.shape[0] / num_tiles)
    num_rows = max(row_dict.values()) + 1
    num_cols = max(col_dict.values()) + 1
    size_x = size_y = 2048
    overlap_x = overlap_y = int(size_x*0.05)
    updated_tile_size = size_x - overlap_x * 2

    ### get stitched image for dapi
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
    # tiff.imwrite(f'{root_folder}/Sample_{sample}/dapi_images/{sample}_{r}_dapi_all_v3.tif', img_merged)
    tiff.imwrite(f'{root_folder}/dapi_images/{r}_dapi_all_v3.tif', img_merged)

    plt.figure(figsize = (2*num_cols, 2*num_rows))
    plt.imshow(img_merged, cmap = 'Blues', vmax = 75)

    # draw tile borders
    for row in range(num_rows):
        plt.plot([0, num_cols * updated_tile_size], [row * updated_tile_size, row * updated_tile_size], color = 'black', linewidth = 0.25)
    for col in range(num_cols):
        plt.plot([col * updated_tile_size, col * updated_tile_size], [0, num_rows * updated_tile_size], color = 'black', linewidth = 0.25)
    
    # add tile text
    for tile in tiles:
        row, col = row_dict[tile+1], col_dict[tile+1]
        plt.text(col * updated_tile_size + 10, row * updated_tile_size + 35, str(tile+1), color = 'black', fontsize = 10)

    plt.axis('off')
    # plt.savefig(f'{root_folder}/Sample_{sample}/dapi_images/{sample}_{r}_dapi_all_v3.png', dpi = 300, bbox_inches = 'tight')
    plt.savefig(f'{root_folder}/dapi_images/{r}_dapi_all_v3.png', dpi = 300, bbox_inches = 'tight')

if __name__ == '__main__':
    samples = ['E2','F','G','H','I','J','K','L','M','N','O','P']
    rs = ['R1','R2','R3']

    for sample in samples:
        if not os.path.exists(f'/data2/haochun/Sennet/ISH/{sample}/dapi_images'):
            os.mkdir(f'/data2/haochun/Sennet/ISH/{sample}/dapi_images')
        for r in rs:
            print(f'{sample}_{r}')
            _stitch_dapi_images(sample, r, vmax = 100, overlap_pixels = 50)
