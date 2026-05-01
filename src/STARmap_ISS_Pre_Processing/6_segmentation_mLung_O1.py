dapi_file = 'dapi_all_v3.tif'
transcript_file = 'spots_decoded_all_stitched_globalcoords_v1.txt'



import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from skimage.filters import gaussian
from skimage.filters import threshold_otsu
from skimage.morphology import binary_dilation, disk
from skimage.segmentation import watershed
import scipy.ndimage as ndi
from stardist.models import StarDist2D
from csbdeep.utils import normalize
from skimage.measure import regionprops
import cv2

def read_tif(file):
    img = cv2.imread(file, -1)
    return img

StarDist2D.from_pretrained()
# creates a pretrained model
model = StarDist2D.from_pretrained('2D_versatile_fluo')

img = read_tif(dapi_file)
pi_label, _ = model.predict_instances(normalize(img))

# Get cell locations 
centroids = []
areas = []
for i, region in enumerate(regionprops(pi_label)):
    centroids.append(region.centroid)
    areas.append(region.area)

# Get cell locations 
centroids = []
areas = []
for i, region in enumerate(regionprops(pi_label)):
    centroids.append(region.centroid)
    areas.append(region.area)
centroids = np.array(centroids)
areas = np.array(areas)

# Apply threshold to centroids based on area
# to_keep = areas > 1000 # area_threshold
# centroids = centroids[to_keep, :]

blurred_overlay_seg = gaussian(img.astype(np.float32), 10)
threhold = threshold_otsu(blurred_overlay_seg)
blurred_overlay_seg = blurred_overlay_seg > threhold
blurred_overlay_seg = binary_dilation(blurred_overlay_seg, footprint=disk(10))

print("Assigning markers")
centroids = centroids.astype(int)
markers = np.zeros(blurred_overlay_seg.shape, dtype=np.uint8)
for i in range(centroids.shape[0]):
    x, y = centroids[i, :]
    if x < blurred_overlay_seg.shape[0] and y < blurred_overlay_seg.shape[1]:
        markers[x-1, y-1] = 1
markers = ndi.label(markers)[0]

print("Watershed")
segmentation = watershed(blurred_overlay_seg, markers, mask=blurred_overlay_seg)
labels_line = watershed(blurred_overlay_seg, markers, mask=blurred_overlay_seg, watershed_line=True)
print(f"Labeled {len(np.unique(segmentation)) - 1} cells")


import pandas as pd
import tqdm

transcript = pd.read_csv(transcript_file, sep='\t', index_col=0)
transcript = transcript.dropna(subset=['target'])

max_y, max_x = segmentation.shape
transcript['y'] = transcript['y'].apply(lambda x: max_y-1 if x >= max_y else x)
transcript['x'] = transcript['x'].apply(lambda x: max_x-1 if x >= max_x else x)

cell_ids = []
for i in tqdm.tqdm(range(len(transcript))):
    x = transcript.iloc[i]['x']
    y = transcript.iloc[i]['y']
    cell = segmentation[max_y-1 -y, x]
    cell_ids.append(cell)
transcript['cell_id'] = cell_ids



import numpy as np

cells = np.unique(transcript['cell_id'])
cells = np.delete(cells, 0)
print(f'Number of segmented cells: {len(cells)}')
genes = np.unique(transcript['target'])
print(f'Number of detected genes: {len(genes)}')
expressions = pd.DataFrame(index=cells, columns=genes, data=0)
df_group = transcript.groupby('cell_id')

for cell in tqdm.tqdm(cells):
    cell_transcript = df_group.get_group(cell)
    gene_exp = cell_transcript['target'].value_counts()
    expressions.loc[cell, gene_exp.index] = gene_exp.values
# give index name as cell
expressions.index.name = 'cell'





import tifffile
# save the transcript, expression matrix and segmentation
transcript = transcript[transcript['cell_id'] != 0] # delete the background

transcript.to_csv('new_spots_decoded_all_stitched_globalcoords_v1.txt', sep='\t')
expressions.to_csv('new_cell_gene_matrix.csv')
tifffile.imwrite('new_cell_segmentation.tif', segmentation.astype(np.uint16))
