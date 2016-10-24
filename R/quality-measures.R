### quality-measures.R
# version 0.1
# (c) 2016 Gregory Breard, University of Rhode Island
#
# This file contains a set of functions used for 
# evaluating the quality of self-organizing maps (SOMs).
### License
# This program is free software; you can redistribute it 
# and/or modify it under the terms of the GNU General 
# Public License as published by the Free Software 
# Foundation.
#
# This program is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the 
# implied warranty of MERCHANTABILITY or FITNESS FOR A 
# PARTICULAR PURPOSE. See the GNU General Public License
# for more details.
#
# A copy of the GNU General Public License is available 
# at: http://www.r-project.org/Licenses/
###

### get.distances - returns the distances between the 
#                   data, neurons, and across both.
#
# parameters:
# - map is an object returned by map.build
# return
# - list containing three distance matrices:
#   dist.data - distance between all data points
#   dist.neurons - distance between all map neurons
#   dist.cross - distance between the data points and the
#               map neurons
#   dist.proj - distance between the projected data points
#
get.distances <- function(map) {
  if (class(map) != "map")
    stop("get.distances: not a map object")
  
  # Get the data set
  data.df <- data.frame(map$data)
  
  # Get the neurons
  neurons.df <- data.frame(map$neurons)
  
  # Merge the data
  colnames(neurons.df) <- colnames(data.df)
  all <- rbind(data.df, neurons.df)
  
  # Calculate the distances
  d <- as.matrix(dist(all))
  
  # Pull out the distances between data and neurons
  n <- dim(data.df)[1]
  m <- dim(neurons.df)[1]
  dist.data <- d[1:n, 1:n]
  dist.neurons <- d[(n + 1):(n + m), (n + 1):(n + m)]
  dist.cross <- d[1:n, (n + 1):(n + m)]
  
  # Get the projected points and distances
  projection <- map$neurons[map$visual,]
  dist.proj <- as.matrix(dist(projection))
  
  # Return distances
  list(dist.data = dist.data,
       dist.neurons = dist.neurons,
       dist.cross = dist.cross,
       dist.proj = dist.proj)
} # end get.distances

### get.distances - returns the ratio of neurons that
#               are best matching units for a data point.
#
# parameters:
# - dist.cross is the distance matrix between the data
#   points and the map neurons
# return
# - list containing a value and a vector:
#   ratio - the ratio of neurons that are a BMU for a 
#           data point.
#   neurons - number of data points mapped to each neuron
#
get.bmu.ratio <- function(dist.cross) {
  # Initialize
  n <- dim(dist.cross)[1]
  m <- dim(dist.cross)[2]
  neurons <- rep(0, m)
  
  # Check each data point
  for (i in 1:n) {
    between <- dist.cross[i, ]
    bmu.idx = which.min(between)
    neurons[bmu.idx] = neurons[bmu.idx] + 1;
  } # end for
    
  # Calculate the error
  bmus <- 0;
  for (i in 1:m)
    if (neurons[i] > 0)
      bmus <- bmus + 1
  ratio <- bmus / m
  
  # Return list
  list(ratio = ratio, neurons = neurons)
} # end get.bmu.ratio


### get.map.diff - returns the difference between two
#                  maps, i.e. the average distance between
#                  the closest pairs of neurons.
#
# parameters:
# - map1 first map
# - map2 second map
# return
# - list containing a value:
#   map.diff - the difference in the maps
#
get.map.diff <- function(map1, map2) {
  # Make sure the maps have the same dimensions
  if (map1$xdim != map2$xdim || map1$ydim != map2$ydim)
    stop("maps have different dimesions")
  
  # Get the neurons from the maps
  neurons.1.df <- data.frame(map1$neurons)
  neurons.2.df <- data.frame(map2$neurons)
  
  # Make sure the neurons have the same number of dims
  if (dim(neurons.1.df)[2] != dim(neurons.2.df)[2])
    stop("map neurons have different dimesions")
  
  # Get the number of neurons (for extracting dist)
  n <- dim(neurons.1.df)[1]
  
  # Merge the data
  colnames(neurons.1.df) <- colnames(neurons.2.df)
  all <- rbind(neurons.1.df, neurons.2.df)
  
  # Calculate the distances
  d <- as.matrix(dist(all))
  
  # Pull out the distances between neurons in two maps
  dist.cross <- d[1:n, (n + 1):(n + n)]
  
  # Use quantization error to find the distance
  # between the nodes in the two maps
  m.diff <- get.quant.err(dist.cross)$val
  
  # Return list
  list(val = m.diff)
} # end get.map.diff
