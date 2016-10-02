### quality-measures.R
# version 0.1
# (c) 2016 Gregory Breard, University of Rhode Island
#
# This file constitues a set of functions used for evaluating the
# quality of self-organizing maps (SOMs).
### License
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
###

### get.distances - returns the distances between the data points, neurons, and across both.
#
# parameters:
# - map is an object returned by map.build
# return
# - list containing three distance matrices:
#   dist.data - distance between all data points
#   dist.neurons - distance between all map neurons
#   dist.cross - distance between the data points and the map neurons
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
