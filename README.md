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
     - has_poi: 0 or 1 if the cell has a POI?

### How to make it work?


