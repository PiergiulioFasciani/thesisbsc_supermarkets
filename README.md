# thesisbsc_supermarkets
Codebase for my BSc Thesis in Urban Economics &amp; Industrial Organization 





### How does it work?

This codebase is a sampler+analyzer that uses OSM data from the Italian Northwest to find:
  - Small Shops (bakeries, butcheries, greengroceries and convenience stores), Points of Interest (POIs);
  - Supermarkets;
  - Funeral Services and Graveyards, Points of Control (POCs, not yet implemented);
  - Public Transport (all kinds).

within the Metropolitan City of Milan to answer the following research question:
  - #### How does proximity to a supermarket influence the probability of finding a smaller retail shop that offers a subset of the services a supermarket already offers?

Sampling is carried out through a quadtree recursive procedure (rationale in the updated methodological note):
  - Designate an Area of Interest (AOI);
  - Split the AOI in four quadrants;
  - If in one of the quadrants there is at least one POI/POC, split the quadrant in four;
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
     - pt_cnt_200m: count of public transport stops within 200m of the tile's boundary;
     - pt_cnt_400m: count of public transport stops within 400m of the tile's boundary;
     - super_cnt_300m: count of supermarkets within 300m of the tile's boundary;
     - super_cnt_600m: count of supermarkets within 600m of the tile's boundary
    which is used for econometric analysis through the fullanalysis.r file.

### How to make it work?


