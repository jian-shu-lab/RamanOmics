import os
import pickle
import logging 
import datetime
import numpy as np
import xarray as xr
from pprint import pprint
import multiprocessing as mp
import multiprocessing as mp
import matplotlib.pyplot as plt
import gzip

from starfish import image
from starfish.types import Axes, Levels
from starfish.image import Filter
from starfish.spots import FindSpots
from starfish.spots import DecodeSpots
from starfish import Experiment
from starfish.spots import FindSpots
from starfish.types import TraceBuildingStrategies


logging.basicConfig(
    filename=f'/data/framont/mLung_O1/starfish_registered/Results/logging/run_pipeline_v1_{datetime.datetime.now().strftime("%Y%m%d-%H%M%S")}.log', 
    level=logging.INFO,
    format='%(asctime)s %(message)s'
)

def filter_white_tophat(imgs, masking_radius):
    wth = Filter.WhiteTophat(masking_radius=masking_radius, is_volume=False)
    return wth.run(imgs)

def filter_clip(imgs, p_min):
    # clip = Filter.Clip(p_min=p_min)
    clip = Filter.ClipPercentileToZero(p_min=p_min, p_max=100, level_method=Levels.SCALE_BY_CHUNK)
    return clip.run(imgs)

def match_histograms_(imgs):
    print("Matching histograms")
    mh_c = Filter.MatchHistograms({Axes.CH})
    scaled_c = mh_c.run(imgs, in_place=False, n_processes=2)
    return scaled_c

def register(imgs, nuclei):
    projection = imgs.reduce([Axes.CH, Axes.ZPLANE], func="max")
    reference_image = nuclei.reduce(dims = [Axes.CH, Axes.ZPLANE], func="max")

    ltt = image.LearnTransform.Translation(
        reference_stack=reference_image,
        axes = Axes.ROUND,
        upsampling=1000,
    )
    transforms = ltt.run(projection)
    warp = image.ApplyTransform.Warp()
    imgs = warp.run(
        stack=imgs,
        transforms_list=transforms,
    )
    return imgs

def equalize_signals(imgs):
    mh = image.Filter.MatchHistograms({Axes.CH, Axes.ROUND})
    scaled = mh.run(imgs, in_place=False, verbose=True, n_processes=2)
    return scaled

def find_spots(scaled, method = 'BlobDetector'):
    ref = scaled.reduce({Axes.CH, Axes.ROUND}, func="max")
    if method == 'BlobDetector':
        bd = FindSpots.BlobDetector(
            min_sigma=0.6,
            max_sigma=10,
            num_sigma=10,
            threshold=0.1,
            exclude_border=0,
            is_volume=True,
            overlap=0.75,
            measurement_type='mean',
        )
        spots = bd.run(scaled, reference_image=ref)
    return spots

def get_codebook():
    file_expr = f'/data/framont/mLung_O1/starfish_registered/tile_data_starmap/tile_0/experiment.json'
    experiment = Experiment.from_json(file_expr)
    codebook = experiment.codebook.copy()
    # channels 1-4: 647, 488, 546, 594 so we use the next parameter
    codebook.values = codebook.sel(c = [3,0,1,2]).values
    return codebook

def decode_spots(spots, method = 'PerRoundMaxChannel'):
    codebook = get_codebook()
    if method == 'PerRoundMaxChannel':
        decoder = DecodeSpots.PerRoundMaxChannel(
            codebook=codebook,
            trace_building_strategy=TraceBuildingStrategies.EXACT_MATCH
        )
    return decoder.run(spots=spots)

def _run_pipeline(
        tile_id, 
        method_findspots = 'BlobDetector',
        method_decodespots = 'PerRoundMaxChannel',
        folder_name = f'/data/framont/mLung_O1/starfish_registered/Results/tile_decoded/'
    ):
    
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
    
    ### load data
    logging.info(f'tile {tile_id}, load data')
    experiment = Experiment.from_json(f'/data/framont/mLung_O1/starfish_registered/tile_data_starmap/tile_{tile_id}/experiment.json')
    fov = experiment.fov()
    registered_imgs = fov.get_image('primary')

    ### match histogram if channels are not balanced
    # logging.info(f'tile {tile_id}, match histogram')
    # registered_imgs = match_histograms_(registered_imgs)

    ### find spots
    logging.info(f'tile {tile_id}, find spots')
    spots = find_spots(registered_imgs, method_findspots)
    with gzip.open(f'{folder_name}tile_{tile_id}_spots.pkl.gz', 'wb') as f:
        pickle.dump(spots, f)
    
    ### decode
    logging.info(f'tile {tile_id}, decode genes')
    decoded = decode_spots(spots, method_decodespots)
    with gzip.open(f'{folder_name}tile_{tile_id}_decoded_v2.pkl.gz','wb') as f:
        pickle.dump(decoded, f)


if __name__ == '__main__':
    tiles_id = [0] # tiles
    tiles_id = [i for i in tiles_id if not os.path.exists(f'/data/framont/mLung_O1/starfish_registered/Results/tile_decoded/tile_{i}_spots.pkl')]

    with mp.Pool(processes = 2) as pool:
        items = list(tiles_id)
        pool.map(_run_pipeline, items)
