library(scatterplot3d)
library(rgl)
library(matlib)
library(plotly)

trials <- c("data/ur5testresult_fullspeed_payload1.6lb_1.csv",
            "data/ur5testresult_fullspeed_payload1.6lb_2.csv",
            "data/ur5testresult_fullspeed_payload1.6lb_3.csv",
            "data/ur5testresult_fullspeed_payload4.5lb_1.csv",
            "data/ur5testresult_fullspeed_payload4.5lb_2.csv",
            "data/ur5testresult_fullspeed_payload4.5lb_3.csv",
            "data/ur5testresult_halfspeed_payload1.6lb_1.csv",
            "data/ur5testresult_halfspeed_payload1.6lb_2.csv",
            "data/ur5testresult_halfspeed_payload1.6lb_3.csv",
            "data/ur5testresult_halfspeed_payload4.5lb_1.csv",
            "data/ur5testresult_halfspeed_payload4.5lb_2.csv",
            "data/ur5testresult_halfspeed_payload4.5lb_3.csv",
            "data/ur5testresult_coldstart_halfspeed_payload4.5lb_1.csv",
            "data/ur5testresult_coldstart_halfspeed_payload4.5lb_2.csv",
            "data/ur5testresult_coldstart_fullspeed_payload4.5lb_3.csv",
            "data/ur5testresult_coldstart_fullspeed_payload4.5lb_1.csv",
            "data/ur5testresult_coldstart_fullspeed_payload4.5lb_2.csv",
            "data/ur5testresult_coldstart_fullspeed_payload4.5lb_3.csv")

df <- data.frame(time = seq(0, 47.992, length = 6000))

for (name in trials) {
  ur5 <- read.csv(file = name, header = FALSE)
  ur5.headers <- scan("data/ur5_header.csv", sep = ',', what = "", quiet = TRUE)
  colnames(ur5) <- ur5.headers
  x <- ur5$`ROBOT_CARTESIAN_COORD_TOOL (x)`
  y <- ur5$`ROBOT_CARTESIAN_COORD_TOOL (y)`
  z <- ur5$`ROBOT_CARTESIAN_COORD_TOOL (z)`
  rx <- ur5$`ROBOT_CARTESIAN_COORD_TOOL (rx)`
  ry <- ur5$`ROBOT_CARTESIAN_COORD_TOOL (ry)`
  rz <- ur5$`ROBOT_CARTESIAN_COORD_TOOL (rz)`

  x <- gsub("[^0-9\\.]", "", x) # cleans non-numeric characters
  x <- as.numeric(x)

  rz <- gsub("[^0-9\\.]", "", rz) # cleans non-numeric characters
  rz <- as.numeric(rz)

  rx <- rx[1:6000]
  ry <- ry[1:6000]
  rz <- rz[1:6000]

  df$placeholder_name = rx
  names(df)[names(df) == "placeholder_name"] <- paste(name, "rx")

  df$placeholder_name = ry
  names(df)[names(df) == "placeholder_name"] <- paste(name, "ry")

  df$placeholder_name = rz
  names(df)[names(df) == "placeholder_name"] <- paste(name, "rz")

  x <- x[1:6000]
  y <- y[1:6000]
  z <- z[1:6000]

  df$placeholder_name = x
  names(df)[names(df) == "placeholder_name"] <- paste(name, "x")

  df$placeholder_name = y
  names(df)[names(df) == "placeholder_name"] <- paste(name, "y")

  df$placeholder_name = z
  names(df)[names(df) == "placeholder_name"] <- paste(name, "z")

}

scatterplot3d(x, y, z)

r1.6 <- cbind(df$`data/ur5testresult_fullspeed_payload1.6lb_3.csv rx`,
              df$`data/ur5testresult_fullspeed_payload1.6lb_3.csv ry`,
              df$`data/ur5testresult_fullspeed_payload1.6lb_3.csv rz`)


r4.5 <- cbind(df$`data/ur5testresult_fullspeed_payload4.5lb_3.csv rx`,
              df$`data/ur5testresult_fullspeed_payload4.5lb_3.csv ry`,
              df$`data/ur5testresult_fullspeed_payload4.5lb_3.csv rz`)

so3 <- SpecialOrthogonalGroup(n = 3)
group.log <- so3$GroupLog(point = r4.5, base.point = r1.6)

jacobian.r1.6 <- so3$JacobianTranslation(r1.6)
for (n in 1:6000){
  group.log[n,] <- group.log[n,] %*% solve(jacobian.r1.6[n,,])
}

plot_ly(
  type = "cone",
  x = c(1, 0, 0),
  y = c(0, 1, 0),
  z = c(0, 0, 1),
  u = c(-sqrt(mean((group.log[,1]) ^ 2)), 0, 0),
  v = c(0, -sqrt(mean((group.log[,2]) ^ 2)), 0),
  w = c(0, 0, -sqrt(mean((group.log[,3]) ^ 2))),
  sizemode = "absolute",
  cmin = 0,
  cmax = max(sqrt(mean((group.log[,1]) ^ 2)),
             sqrt(mean((group.log[,2]) ^ 2)),
             sqrt(mean((group.log[,3]) ^ 2))),

  anchor = "base")%>%
  add_trace(type = 'scatter3d', x = c(0,1), y = c(0,0), z = c(0,0), mode = 'lines',
            opacity = 1, line = list(width = 6, reverscale = FALSE), showlegend = FALSE)%>%
  add_trace(type = 'scatter3d', x = c(0,0), y = c(0,1), z = c(0,0), mode = 'lines',
            opacity = 1, line = list(width = 6, reverscale = FALSE), showlegend = FALSE)%>%
  add_trace(type = 'scatter3d', x = c(0,0), y = c(0,0), z = c(0,1), mode = 'lines',
            opacity = 1, line = list(width = 6, reverscale = FALSE), showlegend = FALSE)%>%
  layout(
    title = "Average Group Log in Rotation Vector: 4.5lb vs 1.6lb load",
    scene = list(
      xaxis = list(
        title = "Y",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE
      ),
      yaxis = list(
        title = "X",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE
      ),
      zaxis = list(
        title = "Z",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE
      )
    ))

rcold <- cbind(df$`data/ur5testresult_coldstart_fullspeed_payload4.5lb_1.csv rx`,
                df$`data/ur5testresult_coldstart_fullspeed_payload4.5lb_1.csv ry`,
                df$`data/ur5testresult_coldstart_fullspeed_payload4.5lb_1.csv rz`)

rwarm <- cbind(df$`data/ur5testresult_fullspeed_payload4.5lb_1.csv rx`,
                df$`data/ur5testresult_fullspeed_payload4.5lb_1.csv ry`,
                df$`data/ur5testresult_fullspeed_payload4.5lb_1.csv rz`)

so3 <- SpecialOrthogonalGroup(n = 3)
group.log <- so3$GroupLog(point = rcold, base.point = rwarm)

jacobian.rwarm <- so3$JacobianTranslation(rwarm)
for (n in 1:6000){
  group.log[n,] <- group.log[n,] %*% solve(jacobian.rwarm[n,,])
}

plot_ly(
  type = "cone",
  x = c(1, 0, 0),
  y = c(0, 1, 0),
  z = c(0, 0, 1),
  u = c(-sqrt(mean((group.log[,1]) ^ 2)), 0, 0),
  v = c(0, -sqrt(mean((group.log[,2]) ^ 2)), 0),
  w = c(0, 0, -sqrt(mean((group.log[,3]) ^ 2))),
  sizemode = "absolute",
  cmin = 0,
  cmax = max(sqrt(mean((group.log[,1]) ^ 2)),
             sqrt(mean((group.log[,2]) ^ 2)),
             sqrt(mean((group.log[,3]) ^ 2))),

  anchor = "base")%>%
  add_trace(type = 'scatter3d', x = c(0,1), y = c(0,0), z = c(0,0), mode = 'lines',
            opacity = 1, line = list(width = 6, reverscale = FALSE), showlegend = FALSE)%>%
  add_trace(type = 'scatter3d', x = c(0,0), y = c(0,1), z = c(0,0), mode = 'lines',
            opacity = 1, line = list(width = 6, reverscale = FALSE), showlegend = FALSE)%>%
  add_trace(type = 'scatter3d', x = c(0,0), y = c(0,0), z = c(0,1), mode = 'lines',
            opacity = 1, line = list(width = 6, reverscale = FALSE), showlegend = FALSE)%>%
  layout(
    title = "Average Group Log in Rotation Vector: Warm vs Cold Start",
    scene = list(
      xaxis = list(
        title = "Y",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE
      ),
      yaxis = list(
        title = "X",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE
      ),
      zaxis = list(
        title = "Z",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE
      )
    ))
