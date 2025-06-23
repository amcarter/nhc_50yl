// Load watershed asset
var watershed = ee.FeatureCollection('users/vlahm13/NHC/nhc_wbd');

// Define date range
var startDate = '1960-01-01';
var endDate = ee.Date(Date.now()).format('YYYY-MM-dd');

// Load PRISM monthly data
var prism = ee.ImageCollection('OREGONSTATE/PRISM/AN81m')
              .filterDate(startDate, endDate)
              .filterBounds(watershed)
              .select(['ppt', 'tmean']);

// Reduce to median value over the watershed
var stats = prism.map(function(img) {
  var date = img.date().format('YYYY-MM-dd');
  var reduced = img.reduceRegion({
    reducer: ee.Reducer.median(),
    geometry: watershed.geometry(),
    scale: 4000,
    maxPixels: 1e13
  });
  return ee.Feature(null, {
    'date': date,
    'ppt': reduced.get('ppt'),
    'tmean': reduced.get('tmean')
  });
});

// Export to CSV
Export.table.toDrive({
  collection: stats,
  description: 'PRISM_PPT_TMEAN_Median',
  fileFormat: 'CSV'
});
