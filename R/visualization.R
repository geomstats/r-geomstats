requireNamespace(plotly, quietly = TRUE)
requireNamespace(reticulate, quietly = TRUE)
requireNamespace(rlist)

use_python('/usr/local/bin/python3', required = T)
gs <- import("geomstats")
h2 <- gs$hyperbolic_space$HyperbolicSpace(dimension=2L)

PoincareDisk <- setRefClass("PoincareDisk",
                        fields = c("center", "points"),
                        methods = list(

                          initialize = function(center=c(0., 0., 0.)){
                            .self$center <- array(c(0., 0., 0.))
                            .self$points <- data.frame(NULL)
                          },

                          ConvertToPoincareCoordinates = function(points){
                            poincare.coords <- (points[, 2:dim(points)[2]]) / (1 + points[, 1])
                            return(poincare.coords)
                          },

                          AddPoints = function(points){
                            poincare_coords <- ConvertToPoincareCoordinates(points)
                            rbind(.self$points, poincare.coords)
                          },

                          GetCircle = function(center, radius, n.points=100){
                            t <- seq(0, 1, 1 / n.points)
                            x <- center[1] + radius * cos(2 * pi * t)
                            y <- center[2] + radius * sin(2 * pi * t)
                            circle <- data.frame(x=x, y=y)
                            return(circle)
                          },

                          GetGeodesicBall = function(center, radius, n.points=100){
                            angles = seq(0, 2 * pi, 2 * pi / n.points)
                            n.angles = length(angles)
                            vectors = array(0, c(n.angles, 3))

                            for (k in 1:n.angles){
                              angle= angles[k]
                              vectors[k,] = c(
                                0.,
                                radius * cos(angle),
                                radius * sin(angle))
                              }

                            tangent_vectors = h2$projection_to_tangent_space(
                              vectors, base_point=center)
                            ball <- h2$metric$exp(tangent_vectors, base_point=center)

                            return(ball)
                          },

                          Draw = function(p){
                             circle <- GetCircle(center=c(0, 0), radius=1)
                             p <- add_trace(p, x=circle$x, y=circle$y,
                                          name = 'Boundary of Poincare Disk',
                                          hoverinfo='skip',
                                          type = 'scatter',
                                          mode = 'lines',
                                          line = list(width=2, color='black'))

                             p <- layout(p,
                               title = "Poincare Disk",
                               xaxis = list(
                                 domain = c(-1, 1),
                                 title = "",
                                 zeroline = FALSE,
                                 showline = FALSE,
                                 showticklabels = FALSE,
                                 showgrid = FALSE
                               ),
                               yaxis = list(
                                 scaleanchor = "x",
                                 domain = c(-1, 1),
                                 title = "",
                                 zeroline = FALSE,
                                 showline = FALSE,
                                 showticklabels = FALSE,
                                 showgrid = FALSE
                               ))
                             return(p)
                          },

                          DrawPoints = function(p, data, col='rgba(255, 182, 193, 1)', size=10){
                            mean <- h2$metric$mean(data)
                            variance <- h2$metric$variance(data, base_point = mean)
                            std.dev <- sqrt(variance)

                            data.poincare <- ConvertToPoincareCoordinates(data)
                            mean.poincare <- ConvertToPoincareCoordinates(mean)
                            data.df <- data.frame(data.poincare)
                            mean.df <- data.frame(mean.poincare)

                            p <- Draw(p)

                            p <- add_trace(p, x=data.df$X1, y=data.df$X2,
                                        name = 'Data',
                                        type = 'scatter',
                                        mode = 'markers+line',
                                        marker = list(
                                          size = size,
                                          color = col))

                            p <- add_trace(p, x=mean.df$mean.poincare[1], y=mean.df$mean.poincare[2],
                                           name = 'Mean',
                                           type = 'scatter',
                                           mode = 'markers',
                                           marker = list(
                                             size = 10,
                                             color = 'rgba(0, 0, 93, .9)'))

                            p <- DrawGeodesicBall(p,
                                        center=mean,
                                        radius=std.dev,
                                        name='Geodesic ball centered at the mean and of radius the standard deviation.')

                            p <- layout(p,
                                        xaxis = list(
                                          hoverformat = '.2f'),
                                        yaxis = list(
                                          hoverformat = '.2f')
                                        )

                            return(p)
                          },

                          DrawGeodesicBall = function(p, center, radius, name='Geodesic Ball', col='rgba(0, 0, 93, .9)', size=2){
                            ball <- GetGeodesicBall(center=center, radius=radius)
                            ball <- ConvertToPoincareCoordinates(ball)
                            ball.df <- data.frame(ball)
                            p <- add_trace(p, x=ball.df$X1, y=ball.df$X2,
                                           name = name,
                                           hoverinfo='skip',
                                           type = 'scatter',
                                           mode = 'markers+line',
                                           marker = list(
                                             size = size,
                                             color = col))
                            return(p)
                          })
)

