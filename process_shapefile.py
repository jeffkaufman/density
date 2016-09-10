# Usage: python process_shapefile.py tabblock2010_??_pophu.shp tl_2010_??_place10.zip

# Derived from makedots.py at https://github.com/meetar/dotmap, which either was written
# by Brendan Martin-Anderson, or Peter Richardson.  Since that's under the MIT license, so
# is this.

import sys
import os
from random import uniform
import ogr
from osgeo import osr
from collections import defaultdict
from math import radians, cos, sin, asin, sqrt
import subprocess
from collections import defaultdict

# If you set this to higher numbers you'll get more accurate results
# but they'll take longer.
OVERSAMPLING = 1

MINDIST = 402.336  # 1/4 mile in meters
ROUGHDIST = 0.009  # distance in latlng space that's > 1/4 mi everywhere
def latlngkey(lat, lng):
  # this rounds to a box 0.01 x 0.01 in latlng space, which will always be at least 1/4 mile in each direction
  return (int(lat * 100), int (lng*100))

# http://stackoverflow.com/questions/4913349/
def distance(lat1, lon1, lat2, lon2):
  """
  Calculate the great circle distance between two points 
  on the earth (specified in decimal degrees)
  """
  # convert decimal degrees to radians 
  lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])
  
    # haversine formula 
  dlon = lon2 - lon1 
  dlat = lat2 - lat1 
  a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
  c = 2 * asin(sqrt(a)) 
  r = 6371 # Radius of earth in kilometers. Use 3956 for miles
  return c * r * 1000

def lookup_field_index(lyr, name):
  feat_defn = lyr.GetLayerDefn()
  field_defns = [feat_defn.GetFieldDefn(i)
                 for i in range(feat_defn.GetFieldCount())]

  for i, defn in enumerate(field_defns):
    if defn.GetName() == name:
      return i
  raise Exception("getting %s" % name)

def get_bbox(geom):
  ll = float("inf")
  bb = float("inf")
  rr = float("-inf")
  tt = float("-inf")

  ch = geom.ConvexHull()
  if not ch:
    return None

  bd = ch.GetBoundary()
  if not bd:
    return None

  pts = bd.GetPoints()
  if not pts:
    return None
  
  for x,y in pts:
    ll = min(ll, x)
    rr = max(rr, x)
    bb = min(bb, y)
    tt = max(tt, y)
    
  return (ll, bb, rr, tt)

def make_ogr_point(x, y):
  return ogr.Geometry(wkt="POINT(%f %f)" % (x, y))
  
def generate_points(tabblock_file, points_file):
  ds = ogr.Open(tabblock_file)
  if not ds:
    raise Exception("tabblock %s open failed" % tabblock_file)

  lyr = ds.GetLayerByIndex(0)
  lyr.ResetReading()

  pop_field = lookup_field_index(lyr, "POP10")
  
  with open(points_file, "w") as outf:
    n_features = len(lyr)
    for j, feat in enumerate( lyr ):
      if j % 10000 == 0:
        print " %s/%s (%0.2f%%)" % (
          j+1, n_features, 100*((j+1)/float(n_features)))
      elif j % 100 == 0:
        sys.stdout.write(".")
        sys.stdout.flush()
        outf.flush()

      pop = feat.GetField(pop_field)
      if pop < 1:
        continue

      geom = feat.GetGeometryRef()
      if geom is None:
        continue

      bbox = get_bbox(geom)
      if not bbox:
        continue
      ll,bb,rr,tt = bbox

      # generate a sample within the geometry for every person
      for i in range(pop * OVERSAMPLING):
        while True:
          samplepoint = make_ogr_point(uniform(ll, rr), uniform(bb, tt))
          if geom.Intersects(samplepoint):
            # lat, lng
            outf.write("%s %s\n" % (samplepoint.GetY(), samplepoint.GetX()))
            outf.flush()
            break

def sort_points(points_file, points_sorted_file):
  subprocess.call(['sort', '-n' , points_file], stdout=open(points_sorted_file, 'w'))

# includes self
def adjacent(key):
  tr = []
  key_lat, key_lng = key
  for i in [-1, 0, 1]:
    for j in [-1, 0, 1]:
      n_key_lat = key_lat + i
      n_key_lng = key_lng + j
      n_key = (n_key_lat, n_key_lng)
      tr.append(n_key)
  return tr  

def compute_density(points_file, points_density_file):
  print "loading points..."
  points = defaultdict(list)
  with open(points_file) as inf:
    for line in inf:
      lat, lng = line.strip().split()
      lat, lng = float(lat), float(lng)
      points[latlngkey(lat, lng)].append((lat, lng))
    
  print "computing density..."
  with open(points_density_file, "w") as outf:
    n_cells = len(points)
    for i, key in enumerate(points):
      print "%s %s/%s (%0.2f%%)" % (
        key, i+1, n_cells, 100*((i+1)/float(n_cells)))
      for adj_key in adjacent(key):
        pass

      for j, (lat1, lng1) in enumerate(points[key]):
        neighbors = 0
        
        for adj_key in adjacent(key):
          if adj_key in points:
            for k, (lat2, lng2) in enumerate(points[adj_key]):
              if adj_key != key or j != k:
                if abs(lat1-lat2) < ROUGHDIST and abs(lng1-lng2) < ROUGHDIST:
                  if distance(lat1, lng1, lat2, lng2) < MINDIST: 
                    neighbors += 1
        outf.write("%s %s %s\n" % (lat1, lng1, float(neighbors)/OVERSAMPLING))  
        outf.flush()
        sys.stdout.write(".")
        sys.stdout.flush()

def compute_town_density(points_density_file, tl_file, points_town_density_file):
  ds = ogr.Open(tl_file)
  if not ds:
    raise Exception("tl open failed: %s" % tl_file)

  lyr = ds.GetLayerByIndex(0)
  lyr.ResetReading()

  name_field = lookup_field_index(lyr, "NAME10")

  points_density = []
  with open(points_density_file) as inf:
    for line in inf:
      lat, lng, neighbors = line.strip().split()
      points_density.append((float(lat), float(lng), float(neighbors)))

  with open(points_town_density_file, "w") as outf:
    n_cities = len(lyr)
    for i, feat in enumerate(lyr):
      geom = feat.GetGeometryRef()
      assert geom
      if geom.IsEmpty():
        continue
      
      name = feat.GetField(name_field)

      print "%s (%s/%s)" % (name, i, n_cities)
      
      for i, (lat, lng, neighbors) in enumerate(points_density):
        if geom.Intersects(make_ogr_point(lng, lat)):
          outf.write("%s %s %s %s\n" % (lat, lng, neighbors, name))
          if i % 100 == 0:
            outf.flush()

def compute_town_averages(points_town_density_file, town_density_file):
  with open(points_town_density_file) as inf:
    with open(town_density_file, "w") as outf:
      prev_name = None
      t_pop = 0
      t_n = 0
      for line in inf:
        lat, lng, neighbors, name = line.strip().split(" ", 4)
        if name != prev_name:
          if t_pop != 0:
            outf.write("%s %s %s %s\n" % (t_pop, t_n, t_n/t_pop, name))
          prev_name = name
          t_pop = 0
          t_n = 0
        t_pop += 1
        t_n += float(neighbors)
      outf.write("%s %s %s %s\n" % (t_pop, t_n, t_n/t_pop, name))

# To get the input files, download:
#    http://www2.census.gov/geo/tiger/TIGER2010/PLACE/2010/tl_2010_${STATE}_place10.zip
#    ftp://ftp2.census.gov/geo/tiger/TIGER2010BLKPOPHU/tabblock2010_${STATE}_pophu.zip
# where $STATE is the number of the state you care about, from 01 to 56.  For example, MA
# is state 25.
def start(tabblock_file, tl_file):
  points_file = "%s.points" % tabblock_file
  points_sorted_file = "%s.points_sorted" % tabblock_file
  points_density_file = "%s.points_density" % tabblock_file
  points_town_density_file = "%s.points_town_density" % tabblock_file
  town_density_file = "%s.town_density" % tabblock_file

  if not os.path.exists(points_file):
    generate_points(tabblock_file, points_file)
  if not os.path.exists(points_sorted_file):
    sort_points(points_file, points_sorted_file)
  if not os.path.exists(points_density_file):
    compute_density(points_file, points_density_file)
  if not os.path.exists(points_town_density_file):
    compute_town_density(points_density_file, tl_file, points_town_density_file)
  if not os.path.exists(town_density_file):
    compute_town_averages(points_town_density_file, town_density_file)

if __name__ == "__main__":
  start(*sys.argv[1:])
