import os


saved_folder='/data/omicseq/Haochun/js_server/SenNet/ISH/O/segmentation/'
dapi_file='/data/omicseq/Haochun/js_server/SenNet/ISH/O/R1_dapi_all_v3.tif'
transcript_file='/data/omicseq/Haochun/js_server/SenNet/ISH/O/spots_decoded_all_stitched_globalcoords_v1.txt'



if not os.path.exists(saved_folder):
    os.makedirs(saved_folder)
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


import tifffile as tiff

tiff.imwrite(os.path.join(saved_folder, "nuclear_mask.tif"), pi_label.astype(np.uint16)

expanded_label = pi_label


import pandas as pd
import tqdm

transcript = pd.read_csv(transcript_file, sep=',')
print(transcript)
transcript = transcript.dropna(subset=['target'])

transcript['y']=img.shape[0]-transcript['y']
max_y,max_x=expanded_label.shape
transcript['y'] = transcript['y'].apply(lambda x: max_y-1 if x >= max_y else x)
transcript['x'] = transcript['x'].apply(lambda x: max_x-1 if x >= max_x else x)

cell_ids = []
for i in tqdm.tqdm(range(len(transcript))):
    x = int(transcript.iloc[i]['x'])
    y = int(transcript.iloc[i]['y'])
    cell = expanded_label[int(max_y) - 1 - y, x]
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

transcript.to_csv(f'{saved_folder}new_spots_decoded_all_stitched_globalcoords_v1.txt', sep='\t')
expressions.to_csv(f'{saved_folder}new_cell_gene_matrix.csv')

