# RamanOmics Source Code

This folder contains the major analysis scripts used in the RamanOmics project, including ISH image processing, STARmap/ISS preprocessing and downstream analysis (Lung as example), Raman preprocessing, imputation analysis, scRNA-seq plot codes, and snRNA-seq analysis (Lung as example).

## Directory structure

### `ISH/`

Scripts for multiplexed ISH image processing, including image stitching, round alignment, transcript coordinate visualization, transcript file stitching, segmentation, and segmentation validation.

Main files:

- `1_stitch_nuclei_images.py`  
  Stitches DAPI/nuclei images across imaging tiles and rounds. It reads microscope XML metadata, reconstructs tile row/column positions, crops overlap regions, stitches tiles into a merged DAPI image, and saves both `.tif` and `.png` outputs.

- `2_stitch_gene_channels.py`  
  Stitches gene-channel images across imaging tiles. This is used to reconstruct full-field multiplexed gene-channel images from tiled microscopy outputs.

- `3_round_alignment 2.py`  
  Aligns imaging rounds, likely to correct spatial shifts between different ISH rounds.

- `4_visualize_multiplexed_ISH_cordinate 1.py`  
  Visualizes multiplexed ISH transcript coordinates, likely for checking decoded transcript locations and spatial patterns.

- `5_stitch_transcript_files.ipynb`  
  Notebook for merging/stitching transcript coordinate files from individual tiles into a whole-sample coordinate table.

- `6_segmentation.py`  
  Performs image-based cell/nuclei segmentation for ISH data.

- `7_Verify_final_segmentation.py`  
  Verifies final segmentation results, plot every gene signals after segmentation on the DAPI background.

---

### `STARmap_ISS_Pre_Processing/`

Scripts for preprocessing STARmap/ISS data before downstream biological analysis.

Main files:

- `1_registration_mLung_O1.py`  
  Performs image registration for mouse lung STARmap/ISS data.

- `2_create_starfish_format_data_mLung_O1.py`  
  Converts registered imaging data into Starfish-compatible format.

- `3_run_starfish_pipeline_mLung_O1.py`  
  Runs the Starfish decoding pipeline to identify spatial transcript spots.

- `4_stitch_decoded_spots_mLung_O1.py`  
  Stitches decoded spots from different tiles into a whole-sample coordinate space.

- `5_stitch_nuclei_mLung_O1.py`  
  Stitches nuclei images (DAPI) across tiles.

- `6_segmentation_mLung_O1.py`  
  Performs segmentation on the stitched mouse lung nuclei (DAPI) data.

---

### `STARmap_ISS_Downstream_Analysis_mLung/`

R scripts for downstream STARmap/ISS analysis of old and young mouse lung samples.

Main files:

- `1_1_mLung_OLD_STARMAP.R`  
  Processes and analyzes old mouse lung STARmap data.

- `1_2_mLung_OLD_STARMAP_LABEL_TRANFER.R`  
  Performs label transfer for old mouse lung STARmap data from scRNA.

- `2_1_mLung_YOUNG_STARMAP.R`  
  Processes and analyzes young mouse lung STARmap data.

- `2_2_mLung_YOUNG_STARMAP_LABEL_TRANFER.R`  
  Performs label transfer for young mouse lung STARmap data from scRNA.

---

### `snRNA_seq_mLung/`

R scripts for mouse lung snRNA-seq preprocessing, integration, annotation, differential expression, GO analysis, and cell-cell communication analysis.

Main files:

- `1_ambient_rna_removal.r`  
  Removes ambient RNA contamination from snRNA-seq data.

- `2_doublet_removal.r`  
  Detects and removes doublets.

- `3_low_qual_filtering.r`  
  Filters low-quality nuclei/cells.

- `4_integration.r`  
  Integrates snRNA-seq datasets across samples/batches.

- `5_scRNA_seq_LUNG_ANNOTATION_Label_Transfer.R`  
  Performs lung cell-type annotation by label transfer using data from GSE124872.

- `6_Senn_DEGs.R`  
  Identifies differentially expressed genes related to senescence analysis.

- `7_mLung_DEGs_Aging.R`  
  Performs aging-associated DEG analysis in mouse lung.

- `8_GO_Loop.R`  
  Runs Gene Ontology enrichment analysis for DEG results.

- `9_CellChat.r`  
  Performs CellChat-based cell-cell communication analysis.

---

### `scRNA/`

R scripts for scRNA-seq downstream analysis, including differential expression, distribution comparison, and CellChat visualization.

Main files:

- `DE_analysis_Skin_Lung_Old_vs_Young_subtype_Plot.R`  
  Performs and visualizes old-vs-young differential expression analysis for skin and lung cell subtypes.

- `EMD.R`  
  Computes Earth Mover's Distance for Lung and Skin.

- `Young_vs_Old_scRNA_Lung_Skin_cellChat_Plot.R`  
  Generates CellChat plots comparing young and old scRNA-seq datasets from lung and skin.

---

### `Imputation_mLung/`

Notebooks for mouse lung imputation analysis.

Main files:

- `mLung_O1_cell.ipynb`  
  Main imputation notebook for mouse lung O1 cell-level analysis with tangram.

- `mLung_O1_cell_LOOV_out.ipynb`  
  Leave-one-out validation or output notebook for mouse lung O1 imputation analysis.

---

### `Raman/`

Notebook for Raman spectral preprocessing.

Main file:

- `preprocess_0.ipynb`  
  Preprocesses Raman spectroscopy data before downstream integration or analysis.

---


