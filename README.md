# thesisbsc_supermarkets
Codebase for my BSc Thesis in Urban Economics &amp; Industrial Organization 





### How does it work?

This codebase is a sampler+analyzer that uses OSM data from the Italian Northwest to find:
  - Small Shops (bakeries, butcheries, greengroceries and convenience stores), Points of Interest (POIs);
  - Supermarkets;
  - Funeral Services and Graveyards, Points of Control (POCs, not yet implemented in econometric analysis);
  - Public Transport (all kinds).

within the Metropolitan City of Milan to answer the following research question:
  - #### How does proximity to a supermarket influence the probability of finding a smaller retail shop that offers a subset of the services a supermarket already offers?

Sampling is carried out through a quadtree recursive procedure (rationale in the updated methodological note):
  - Designate an Area of Interest (AOI) by latitude, longitude and radius;
  - Split the AOI in four quadrants;
  - If in one of the quadrants there is at least two POIs/POCs, split the quadrant in four;
  - Repeat splitting until each quadrant has either one or zero POIs/POCs.

After this procedure is completed, two files are created and stored in the "data" folder:
   - A .gpkg file is created and can be opened in QGIS (open source software) for visualizing the parsing and splitting results;
   - A .csv file is created and it contains:
     - tile_id: numerical ID for the tile obtained through quadtree;
     - depth: number of subdivisions necessary to isolate that tile;
     - n_pois: number of POIs in the tile (1 max, residual from previous experiment, will be fixed later);
     - has_poi: binary variable that states if in a tile there is a POI;
     - n_poc: number of POCs in the tile (1 max, residual from previous experiment, will be fixed later);
     - has_poc: binary variable that states if in a tile there is a POC;
     - center_lon: longitude of the tile's centroid;
     - center_lat: latitude of the tile's centroid;
     - tile_side_m: length in meters of the tile's side;
     - tile_area_m2: area in squared meters of the tile;
     - dist_poly_to_pt_m: distance in meters from the tile's boundary to the closest public transport stop (=0 if it is inside the tile);
     - dist_poly_to_supermarket_m: distance in meters from the tile's boundary to the closest supermarket stop (=0 if it is inside the tile);
     - dist_poly_to_duomo_m: distance in meters from the tile's boundary to the closest Duomo (proxy for city center) stop (=0 if it is inside the tile);
     - dist_poly_to_pt_outside_m: distance in meters from the tile's boundary to the closest public transport stop outside the tile;
     - dist_poly_to_supermarket_outside_m: distance in meters from the tile's boundary to the closest supermarket outside the tile;
     - poi2super_min_m: distance in meters from the tile's POI (if present) to the nearest supermarket;
     - poi2pt_min_m: distance in meters from the tile's POI (if present) to the nearest public transport stop;
     - poi2duomo_min_m: distance in meters from the tile's POI (if present) to Duomo;
     - pt_cnt_200m: count of public transport stops within 200m of the tile's boundary;
     - pt_cnt_400m: count of public transport stops within 400m of the tile's boundary;
     - super_cnt_300m: count of supermarkets within 300m of the tile's boundary;
     - super_cnt_600m: count of supermarkets within 600m of the tile's boundary

After the sampling is carried out, a fullanalysis.r file begins an analysis in the following way:
   - Allows the user to load any .csv file from the sampling at choice;
   - Creates an output folder in the data folder;
   - If polygon-to-POI/Duomo/supermarket/PT = 0, swap to outside-of-tile measurements;
   - Prefers POI-supermarket/PT/Duomo distances, if available, or use polygon distances;
   - Transforms distances into logarithms (reduces heavy-tails, betters handling);
   - Creates three models:
      - Full logistic regression model: has_poi ~ lpt + lsuper + lduomo + squares + interaction + density and area controls;
      - Simpler logistic regression model: has_poi ~ lsuper + lduomo + interactions;
      - Simpler LPM model: has_poi ~ dist_pt_100m + dist_super_100m + dist_duomo_km + interactions + squares;

   - Conducts robust inference twice by:
      - Computing heteroskedasticity-robust SEs on the entire sample;
      - Computing heteroskedasticity-robust SEs on 1km radius clusters per-tile in the entire sample.

   - Compiles a neat summary .txt file found in the data folder.

### How to make it work?

To make this codebase run on your system do the following:

   - Install the following:
      - R 4.2 or above;
      - Python 3.10 or above;

   - On Windows:
      - Double click run_all.bat;
         - 
      - When prompted, search for the latest .csv file in the data folder;
      - At the end, you can see the .txt file with the models' results;
        
   - On Mac/Linux
      - 

 

