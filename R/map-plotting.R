### map-plotting.R
# version 0.1
# (c) 2016 Gregory Breard, University of Rhode Island
#
# This file contains a set of functions used for 
# plottiing self-organizing maps (SOMs).
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

# load packages
library(scatterplot3d)
library(rgl)
library(RColorBrewer)

### map.plot2d - plots the map in 2-dimensions with data
#                point and neurons colored by label
#
# parameters:
# - map is an object returned by map.build
# - x.idx the index of the column to use as x-axis
# - y.idx the index of the column to use as y-axis
# - show.dead logical, if true dead nodes are highlighted.
# - highlight a list of indices of nodes to highlight
#             (overrides show.dead)
# - use.rgl logical, true if the visualization should use
#           the rgl package.
# return:
# - nothing
#
map.plot2d <- function(map, x.idx = 1, y.idx = 2, 
                       show.dead = F, highlight = NULL, 
                       use.rgl = F) {
  n <- dim(map$data)[1]
  
  # need to check if the data has labels
  if (is.null(map$labels)) {
    # get the node labels
    node.labels <- rep(1, dim(map$neurons)[1])
    pal = "white"
    
    # plot data points
    if (use.rgl)
      plot3d(map$data[, x.idx], map$data[, y.idx], 
             rep(0, n),
             col = c("grey"), pch = rep(20, n), 
             cex = rep(0.5, n))
    else
      plot(map$data[, x.idx], map$data[, y.idx], 
           col = c("grey"), pch = 20, cex = 0.5)
  } else {
    # get the node labels
    node.labels <- as.numeric(as.factor(
                        map.label.accuracy(map)$neurons))
    if (max(node.labels, na.rm = T) < 10)
      pal = c(brewer.pal(max(node.labels, na.rm = T), 
                         "Set1"), "white")
    else
      pal = c(rainbow(max(node.labels, na.rm = T)), 
                                                "white")
    na.lab <- max(node.labels, na.rm = T) + 1
    node.labels[which(is.na(node.labels))] = na.lab
    
    # plot data points
    if (use.rgl)
      plot3d(map$data[, x.idx], map$data[, y.idx], 
             rep(0, n), 
             col = pal[as.factor(
                        as.data.frame(map$labels)[,1])],
             pch = rep(20, n), cex = rep(0.5, n))
    else
      plot(map$data[, x.idx], map$data[, y.idx], 
           col = pal[as.factor(
                        as.data.frame(map$labels)[,1])], 
           pch = 20, cex = 0.5)
  } # end if
  
  # plot edges before nodeso we can see the colors
  for (i in 1:dim(map$neurons)[1]) {
    # get the position in the map
    coord <- get.map.coord(i, map$xdim, map$ydim)
    
    # get the neighboring nodes
    hood <- get.hood(coord$x, coord$y, map$xdim, map$ydim)
    
    # add segments between the nodes
    for (j in 1:length(hood)) {
      # convert back to index
      point <- hood[[j]]
      l <- get.map.idx(point[1], point[2], map$xdim)
      if (use.rgl)
        segments3d(c(map$neurons[i, x.idx], 
                     map$neurons[l, x.idx]), 
                   c(map$neurons[i, y.idx], 
                     map$neurons[l, y.idx]), c(0,0))
      else
        segments(map$neurons[i, x.idx], 
                 map$neurons[i, y.idx], 
                 map$neurons[l, x.idx], 
                 map$neurons[l, y.idx])
    } # end for (j)
  } # end for (i)
  
  # check if we should highlight dead nodes
  if (is.null(highlight) && show.dead)
    highlight <- which(!(1:(map$xdim * map$ydim) 
                              %in% unique(map$visual)))
  
  # plot map nodes
  if (is.null(highlight)) {
    if (use.rgl)
      points3d(map$neurons[, x.idx], map$neurons[, y.idx],
               rep(0, n), col = "black", 
               bg = pal[node.labels], 
               pch = rep(21, n), cex = rep(20, n))
    else
      points(map$neurons[, x.idx], map$neurons[, y.idx], 
             col = "black", bg = pal[node.labels], pch = 21)
  } else {
    all.n <- 1:dim(map$neurons)[1]
    h <- all.n %in% highlight
    not.h <- !h
    if (use.rgl) {
      points3d(map$neurons[not.h, x.idx], 
               map$neurons[not.h, y.idx], 
               rep(0, n), col = "black", 
               bg = pal[node.labels[not.h]], 
               pch = rep(21, n))
      points3d(map$neurons[h, x.idx], 
               map$neurons[h, y.idx], rep(0, n), 
               col = "darkgoldenrod1", 
               bg = pal[node.labels[h]], pch = rep(23, n))
    } else {
      points(map$neurons[not.h, x.idx], 
             map$neurons[not.h, y.idx], 
             col = "black", 
             bg = pal[node.labels[not.h]], 
             pch = rep(21, n))
      points(map$neurons[h, x.idx], 
             map$neurons[h, y.idx], 
             col = "darkgoldenrod1", 
             bg = pal[node.labels[h]], pch = rep(23, n))
    } # end if
  } # end if
} # end map.plot2d

### map.plot3d - plots the map in 3-dimensions with data
#                point and neurons colored by label
#
# parameters:
# - map is an object returned by map.build
# - x.idx the index of the column to use as x-axis
# - y.idx the index of the column to use as y-axis
# - z.idx the index of the column to use as z-axis
# - show.dead logical, if true dead nodes are highlighted.
# - highlight a list of indices of nodes to highlight
#             (overrides show.dead)
# - use.rgl logical, true if the visualization should use
#           the rgl package.
# return:
# - nothing
#
map.plot3d <- function(map, x.idx = 1, y.idx = 2, 
                       z.idx = 3, show.data = T, 
                       show.map = T, show.dead = F, 
                       highlight = NULL, use.rgl = F) {
  # reset the device before drawing
  if (use.rgl)
    clear3d()
  
  # need to check if the data has labels
  if (is.null(map$labels)) {
    # get the node labels
    node.labels <- rep(1, dim(map$neurons)[1])
    pal = "white"
    
    # check if we should draw the data points
    if (show.data) {
      # plot data points
      if (use.rgl)
        plot3d(map$data[, x.idx], map$data[, y.idx], 
                      map$data[, z.idx], col = c("grey"))
      else
        s <- scatterplot3d(map$data[, x.idx], 
                           map$data[, y.idx], 
                           map$data[, z.idx], 
                           color = "grey", 
                           pch = 20, cex.symbols = 0.5)
    } # end if  
  } else {
    # get the node labels
    node.labels <- as.numeric(
              as.factor(map.label.accuracy(map)$neurons))
    if (max(node.labels, na.rm = T) < 10)
      pal = c(brewer.pal(max(node.labels, na.rm = T),
                      "Set1"), "white") # get max with NAs
    else
      pal = c(rainbow(max(node.labels, na.rm = T)), 
              "white")
    node.labels[which(is.na(node.labels))] = 
                          max(node.labels, na.rm = T) + 1
    
    # check if we should draw the data points
    if (show.data) {
      # plot data points
      if (use.rgl) {
        plot3d(map$data[, x.idx], map$data[, y.idx], 
               map$data[, z.idx], 
               col = pal[as.factor(
                    as.data.frame(map$labels)[,1])])
      } else
        s <- scatterplot3d(map$data[, x.idx], 
                           map$data[, y.idx], 
                           map$data[, z.idx], 
                           color = pal[as.factor(
                        as.data.frame(map$labels)[,1])],
                           pch = 20, cex.symbols = 0.5)
    } # end if
  } # end if
  
  # check if we should draw the map
  if (show.map) {
    # use grey instead of white for unlabeled with rgl
    if (use.rgl) pal[length(pal)] = "grey"
   
    # plot edges before nodeso we can see the colors
    for (i in 1:dim(map$neurons)[1]) {
      # get the position in the map
      coord <- get.map.coord(i, map$xdim, map$ydim)
      
      # get the neighboring nodes
      hood <- get.hood(coord$x, coord$y, map$xdim, 
                       map$ydim)
      
      # add segments between the nodes
      for (j in 1:length(hood)) {
        # convert back to index
        point <- hood[[j]]
        l <- get.map.idx(point[1], point[2], map$xdim)
        if (use.rgl)
          segments3d(c(map$neurons[i, x.idx], 
                       map$neurons[l, x.idx]), 
                     c(map$neurons[i, y.idx], 
                       map$neurons[l, y.idx]), 
                     c(map$neurons[i, z.idx], 
                       map$neurons[l, z.idx]))
        else {
          p1 <- s$xyz.convert(map$neurons[i, x.idx], 
                              map$neurons[i, y.idx], 
                              map$neurons[i, z.idx])
          p2 <- s$xyz.convert(map$neurons[l, x.idx], 
                              map$neurons[l, y.idx], 
                              map$neurons[l, z.idx])
          segments(p1$x, p1$y, p2$x, p2$y)
        } # end if
      } # end for (j)
    } # end for (i)
      
    # check if we should highlight dead nodes
    if (is.null(highlight) && show.dead)
      highlight <- which(!(1:(map$xdim * map$ydim) 
                                %in% unique(map$visual)))
    
    # plot map nodes
    if (is.null(highlight)) {
      if (use.rgl)
        points3d(map$neurons[, x.idx], 
                 map$neurons[, y.idx], 
                 map$neurons[, z.idx], 
                 col = pal[node.labels], size = 10)
      else {
        nodes.2d <- s$xyz.convert(map$neurons[, x.idx], 
                                  map$neurons[, y.idx], 
                                  map$neurons[, z.idx])
        points(nodes.2d$x, nodes.2d$y, col = "black", 
               bg = pal[node.labels], pch = 21)
      }
    } else {
      all.n <- 1:dim(map$neurons)[1]
      h <- all.n %in% highlight
      not.h <- !h
      if (use.rgl) {
        points3d(map$neurons[not.h, x.idx], 
                 map$neurons[not.h, y.idx], 
                 map$neurons[not.h, z.idx], 
                 col = pal[node.labels[not.h]], size = 2)
        points3d(map$neurons[h, x.idx], 
                 map$neurons[h, y.idx], 
                 map$neurons[not.h, z.idx],
                 col = "darkgoldenrod1", size = 2)
      } else {
        nodes.2d <- s$xyz.convert(
                        map$neurons[not.h, x.idx], 
                        map$neurons[not.h, y.idx], 
                        map$neurons[not.h, z.idx])
        nodes.h.2d <- s$xyz.convert(
                        map$neurons[h, x.idx], 
                        map$neurons[h, y.idx], 
                        map$neurons[h, z.idx])
        points(nodes.2d$x, nodes.2d$y, col = "black", 
               bg = pal[node.labels[not.h]], pch = 21)
        points(nodes.h.2d$x, nodes.h.2d$y, 
               col = "darkgoldenrod1", 
               bg = pal[node.labels[h]], pch = 23)
      } # end if
    } # end if
  } #end if
} # end map.plot3d

### map.label.accuracy - returns the labeling of the 
#         neurons and other measures
#
# parameters:
# - map is an object returned by map.build
# return:
# - list containing:
#   neurons - majority label mapped to each neuron
#   acc - labeling accuracy of the map
#   bmu - ratio of BMUs
#   fit - combined accuracy and ratio value
#
map.label.accuracy <- function(map) {
  # get the neuron labels
  neuron.l <- rep(NA, dim(map$neurons)[1])
  proj <- map$visual
  labels <- as.vector(as.data.frame(map$labels)[,1])
  
  # assign majority labels to neurons
  for (i in 1:length(neuron.l)) {
    if (length(which(proj == i)) > 0) {
      labels.n <- labels[which(proj == i)]
      neuron.l[i] <- names(which.max(table(labels.n)))
    } # end if
  } # end forlength(unique(proj)) / length(neuron.l)
  
  # check how many of the labels of the neurons match 
  # the input and how many neurons are mapped to
  acc <- (length(which(neuron.l[proj] == labels)) 
                                      / length(proj))
  bmu <- length(unique(proj)) / length(neuron.l)
  list(neurons = neuron.l, acc = acc, bmu = bmu, 
       fit = (acc * bmu))
} # end map.label.accuracy

# ----------------- Helper Functions ----------------- #

### get.hood - get the indices of the neighboring nodes
#
# parameters:
# - x is the x index
# - y is the y index
# - xdim is the width of the map
# - ydim is the height of the map
# return:
# - list containing:
#   t - top neighbor coordinates
#   r - right neighbor coordinates
#   b - bottom neighbor coordinates
#   l - left neighbor coordinates
#
get.hood <- function(x, y, xdim, ydim) {
  hood <- list()
  
  # get the neighbors
  t <- c(x, y + 1)
  r <- c(x + 1, y)
  b <- c(x, y - 1)
  l <- c(x - 1, y)
  
  # check if they are valid
  if (t[2] <= ydim) {
    hood[["t"]] = t
  }
  if (r[1] <= xdim) {
    hood[["r"]] = r
  }
  if (b[2] > 0) {
    hood[["b"]] = b
  }
  if (l[1] > 0) {
    hood[["l"]] = l
  }
  
  hood
} # end get.hood

### get.map.coord - gets the x, y position in the map
#
# parameters:
# - i is the index of the neuron
# - xdim is the width of the map
# - ydim is the height of the map
# return:
# - list containing:
#   x - the x index
#   y - the y index
#
get.map.coord <- function(i, xdim, ydim) {
  # get the position in the map
  x <- i %% xdim
  if (x == 0) x <- xdim
  y <- ceiling(i / xdim)
  list(x = x, y = y)
} # end get.map.coord

### get.map.idx - gets the index in the map
#
# parameters:
# - x is the x index
# - xdim is the width of the map
# - ydim is the height of the map
# return:
# - the index
#
get.map.idx <- function(x, y, xdim) {
  # get the index in the map
  ((y - 1) * xdim) + x
} # end get.map.idx

