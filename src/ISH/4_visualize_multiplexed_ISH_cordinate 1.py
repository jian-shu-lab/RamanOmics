import os
import numpy as np
import pandas as pd
import tifffile as tiff
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters

def _find_local_maxima(
    image,
    neighborhood_size=5,
    threshold=50,
):
    data_max = filters.maximum_filter(image, neighborhood_size)
    maxima = (image == data_max)
    data_min = filters.minimum_filter(image, neighborhood_size)
    diff = ((data_max - data_min) > threshold)
    maxima[diff == 0] = 0

    labeled, num_objects = ndimage.label(maxima)
    slices = ndimage.find_objects(labeled)
    x, y = [], []

    for dy, dx in slices:
        x_center = (dx.start + dx.stop - 1) / 2
        x.append(x_center)
        y_center = (dy.start + dy.stop - 1) / 2    
        y.append(y_center)

    return num_objects, x, y



def spots_calling(
    image_ch, 
    channel_name,
    local_maxima_neighbors=10, 
    local_maxima_diff_threshold=50,
):
    num_objects, x, y = _find_local_maxima(
        image_ch,
        neighborhood_size=local_maxima_neighbors,
        threshold=local_maxima_diff_threshold,
    )

    return x, y, num_objects

### params
samples = ['E2','F','G','H','I','J','K','L','M','N','O','P']
Rs = ['R1','R2','R3']
channels = ['ch01', 'ch02', 'ch03', 'ch04']
colors = ['#4B4C76', '#687DA5', '#A8CED1', '#91B894', '#CED87F', '#E2CB6B', '#D1A289', '#9E5B52']
for sample in samples:
    root_folder=f'/data2/haochun/Sennet/ISH/{sample}/'
    spots_folder=root_folder+'spots_images/'
    save_folder = root_folder+'figures_v3/'
    scale = 1 / 300
    dapi_amplification = 5
    coordinates_folder = root_folder+'coordinates_v3/'
    # scale bar
    pixel_size = 581.25/2048  # um
    bar_size = 300  # um
    scale_factor = bar_size / pixel_size  # Adjust this value as needed
    scale_x = 200  # X-coordinate for the scale bar

    if not os.path.exists(save_folder):
        os.mkdir(save_folder)
    if not os.path.exists(coordinates_folder):
        os.makedirs(coordinates_folder, exist_ok=True)



    df_genes = pd.read_excel('/data/Haochun/ISH/Data/Ray_v3_ISH2.xlsx', sheet_name='Ray_v3_ISH_Prob2')
    df_genes.set_index('Gene', inplace=True)

    dapi = tiff.imread(f'{root_folder}/dapi_images/R1_dapi_all_v3.tif')
    fig = plt.figure(figsize=(dapi.shape[1] * scale, dapi.shape[0] * scale))
    dapi = dapi * dapi_amplification
    dapi[dapi > 255] = 255
    dapi = dapi.astype(np.uint8)
    plt.imshow(dapi, cmap='gray', vmax=255, alpha=1.0)

    cell_types = []
    for R in Rs:
        for channel in channels:
            if R == 'R1':
                img = tiff.imread(f'{spots_folder}/{R}_{channel}_all.tif')
            elif (R == 'R3')|( R=='R2'):
                img = tiff.imread(f'{spots_folder}{R}_{channel}_all_aligned.tif')

            channel2 = int(channel.split('ch0')[1])
            
            # Debugging prints
            print(f"Checking for R: {R}, channel: {channel2}")
            
            matching_genes = df_genes[(df_genes['Round'] == R) & (df_genes['Channel'] == channel2)]
            
            # Debugging prints
            print(f"Matching genes: {matching_genes}")
            print(matching_genes)
            if matching_genes.empty:
                print(f"No matching gene found for {R}, {channel2}")
                continue

            gene_name = matching_genes.index[0]

            xo, yo = np.where(img > 0)
            x, y, num_objects = spots_calling(
                image_ch=img,
                channel_name=f'{R}_{channel}_{gene_name}',
                local_maxima_neighbors=15, 
                local_maxima_diff_threshold=100
            )
            print(f'{sample}_{R}_{channel}: {gene_name}, {len(x) / 1e6} million spots out of {len(xo) / 1e6} million foreground pixels')
            color = df_genes.loc[gene_name, 'Color']
            cell_type = df_genes.loc[gene_name, 'Cell type']
            cell_types.append(cell_type)
            label = f'{gene_name} ({cell_type})' if not pd.isna(cell_type) else gene_name
            plt.scatter(x, y, color=color, s=1.5, label=label, alpha=0.8)
            # Save coordinates to a .txt file
            with open(f'{coordinates_folder}/{R}_{gene_name}_coordinates_filtered.txt', 'w') as f:
                f.write('x_coordinate\ty_coordinate\n')
                for x_coord, y_coord in zip(x, y):
                    f.write(f'{x_coord}\t{y_coord}\n')

    # sort them by legend labels
    handles, labels = plt.gca().get_legend_handles_labels()
    sorted_indices = np.argsort(cell_types)
    handles = [handles[i] for i in sorted_indices]
    labels = [labels[i] for i in sorted_indices]

    scale_y = np.max(y) - 200
    plt.plot([scale_x, scale_x + scale_factor], [scale_y, scale_y], color='white', linewidth=5)
    plt.annotate(f'{bar_size}µm', (scale_x + scale_factor / 2, scale_y - 75), color='white', ha='center', fontsize=24)

    plt.legend(handles, labels, markerscale=10, fontsize=20)
    plt.axis('off')
    plt.savefig(f'{save_folder}/{sample}.png', dpi=300, bbox_inches='tight')
