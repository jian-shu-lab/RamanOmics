##If use segmentation, need to flip y in segmentation step and also in this step
from PIL import Image
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import matplotlib.cm as cm
import numpy as np

Image.MAX_IMAGE_PIXELS = None  # Disable decompression bomb check
gene_color_map = {
    "Cdkn1a": "#e9392f",
    "Cdkn2a": "#FF8C00",
    "Krt10":  "#4169E1",
    "Sfn":    "#1E90FF",
    "Lor":    "#9370DB",
    "Dmkn":   "#00BFFF",
    "Eepd1":  "#32CD32",
    "Sbsn":   "#00FF7F",
    "Dab2":   "#3CB371",
    "Scn7a":  "#4B0082",
    "Abcg1":  "#FFD700",
    "Serpine1": "#FF1493"
}

dapi_amplification = 5
scale_bar_um = 300
pixel_size = 581.25/2048  # micron/pixel
scale_factor = scale_bar_um / pixel_size
scale_x = 200  # position for scale bar


def plot_spots_on_tile(image_file, spots_file, output_prefix,plot_all_targets=True):
    
    # Load the image
    img = Image.open(image_file)
    img_array = mpimg.imread(image_file)
    print(img_array.shape)
    
    # Load the spots
    spots_df = pd.read_csv(spots_file, delimiter='\t')
    print(spots_df.columns)
    print(spots_df.shape)
    # Get the image height
    img_height = img.height
    print(img_height)
    spots_df['y_flipped'] = img_height - spots_df['y']
    # Get the unique targets
    unique_targets = spots_df['target'].unique()
    print(unique_targets)

    # Loop through each target and create a plot
    for target in unique_targets:
        target_df = spots_df[spots_df['target'] == target]
        
        # Plot the image as is
        plt.figure(figsize=(20, 20))
        plt.imshow(img_array, cmap='gray')

        # Plot the spots for the specific target
        plt.scatter(target_df['x'], target_df['y_flipped'], c=gene_color_map[target], s=5, label=f'Spots - {target}')
        # plt.scatter(target_df['x'], target_df['y'], c=gene_color_map[target], s=2, label=f'Spots - {target}')

        # Add labels and legend
        plt.xlabel('X Coordinate')
        plt.ylabel('Y Coordinate')
        plt.title(f'Spots on Tile - {target}')
        plt.legend(markerscale=5, bbox_to_anchor=(1.05, 1),loc='upper left')
        # Scale bar for gene
        scale_y = np.max(target_df['y_flipped']) - 200 if len(target_df['y_flipped']) > 0 else img.shape[0] - 200
        # scale_y = np.max(target_df['y']) - 200 if len(target_df['y']) > 0 else img.shape[0] - 200

        plt.plot([scale_x, scale_x + scale_factor], [scale_y, scale_y], color='white', linewidth=5)
        plt.annotate(f'{scale_bar_um}µm', (scale_x + scale_factor/2, scale_y - 75),
                            color='white', ha='center', fontsize=50)
        # Save the plot without any inversion
        output_plot = f"{output_prefix}{target}.png"
        plt.savefig(output_plot)

    # Plot all genes with different colors
    if plot_all_targets:
        plt.figure(figsize=(20, 20))
        plt.imshow(img_array, cmap='gray')

        # Assign a distinct color to each target
        cmap = cm.get_cmap('tab20', len(unique_targets))  # Use a colormap with enough colors

        for i, target in enumerate(unique_targets):
            print(target)
            target_df = spots_df[spots_df['target'] == target]
            color = gene_color_map[target]
            print(color)
            plt.scatter(target_df['x'], target_df['y_flipped'], c=color, s=5, label=target)

        plt.xlabel('X Coordinate')
        plt.ylabel('Y Coordinate')
        plt.title('Spots on Tile - All Genes')
        # Scale bar for gene
        scale_y = np.max(spots_df['y_flipped']) - 200 if len(spots_df['y_flipped']) > 0 else img.shape[0] - 200
        plt.plot([scale_x, scale_x + scale_factor], [scale_y, scale_y], color='white', linewidth=5)
        plt.annotate(f'{scale_bar_um}µm', (scale_x + scale_factor/2, scale_y - 75),
                            color='white', ha='center', fontsize=50)
        plt.legend(markerscale=5, bbox_to_anchor=(1.05, 1),loc='upper left')  # move legend outside
        plt.tight_layout()
        output_plot = f"{output_prefix}All_Genes.png"
        plt.savefig(output_plot, dpi=300)
        plt.close()
        


# Usage
import os
# root_path='/data/Haochun/ISH/Data/Ke_Final/'
root_path='/data2/haochun/Sennet/ISH/'
sample='E2'
if not os.path.exists(f'{root_path}/{sample}/segmentation/final_plot_select'):
        os.mkdir(f'{root_path}/{sample}/segmentation/final_plot_select')
plot_spots_on_tile(f'/data2/haochun/Sennet/ISH/{sample}/nuclear_mask.tif', f'/data/Haochun/mgb_server/SenNet/ISH/{sample}/segmentation/new_spots_decoded_all_stitched_globalcoords_v1.txt', f'{root_path}/{sample}/segmentation/final_plot_select/')
